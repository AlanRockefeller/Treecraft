from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QComboBox, QPushButton, QProgressBar, QMessageBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
import logging

from treecraft.utils.external_tools import (
    get_available_alignment_methods,
    find_external_tool,
    run_alignment_tool
)

logger = logging.getLogger("treecraft.gui.alignment_dialog")

class AlignmentThread(QThread):
    """Thread for running sequence alignment"""
    progress_update = pyqtSignal(int)
    alignment_complete = pyqtSignal(object)
    error_occurred = pyqtSignal(str)
    log_message = pyqtSignal(str)  # New signal for thread-safe logging
    
    def __init__(self, sequences, method_id):
        super().__init__()
        self.sequences = sequences
        self.method_id = method_id  # This is now the internal ID (halign, muscle, simple)
        
        # Use the log_message signal for thread-safe logging
        self.log_message.emit(f"Creating alignment thread with method: {method_id}")
        
        # Don't create or manipulate any GUI objects here or in run()
    
    def run(self):
        try:
            from Bio import AlignIO
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            import tempfile
            import os
            
            # Create temporary files
            input_fd, input_path = tempfile.mkstemp(suffix=".fasta")
            output_fd, output_path = tempfile.mkstemp(suffix=".aln")
            os.close(input_fd)
            os.close(output_fd)
            
            # Write sequences to input file with their full descriptions
            with open(input_path, "w") as f:
                for i, record in enumerate(self.sequences):
                    # Use full description to preserve names
                    # Also include a custom prefix so we can match them back later
                    description = record.description if record.description else record.id
                    
                    # Store the sequence ID at the start to ensure it's preserved for matching
                    f.write(f">{record.id} {description}\n{record.seq}\n")
                    
                    # Update progress for sequence writing
                    self.progress_update.emit(int((i+1) / len(self.sequences) * 30))
            
            self.log_message.emit(f"Wrote {len(self.sequences)} sequences to {input_path}")
            
            # Run alignment based on method
            if self.method_id in ["halign", "muscle", "mafft"]:
                try:
                    # Use the external alignment tool
                    self.progress_update.emit(40)
                    success = run_alignment_tool(self.method_id, input_path, output_path, self)
                    
                    if success:
                        self.progress_update.emit(80)
                        self.log_message.emit(f"Successfully ran {self.method_id}, reading alignment from {output_path}")
                        
                        # Read the aligned sequences
                        alignment = list(AlignIO.read(output_path, "fasta"))
                        
                        # Preserve full descriptions and create new records
                        # This ensures we maintain the original sequence descriptions
                        preserved_alignment = []
                        
                        # Create mapping from original seq IDs to their full descriptions
                        original_descriptions = {}
                        original_id_map = {}
                        
                        for record in self.sequences:
                            # Store both ID and the full description
                            original_descriptions[record.id] = record.description
                            
                            # Also map from ID without spaces to original ID
                            # (halign may truncate at spaces or special chars)
                            base_id = record.id.split()[0]
                            original_id_map[base_id] = record.id
                            
                            # Also try mapping using the first word of the description
                            if record.description:
                                first_word = record.description.split()[0]
                                original_id_map[first_word] = record.id
                                
                            # Also map full description for more precise matching
                            if record.description:
                                original_id_map[record.description] = record.id
                        
                        self.log_message.emit(f"Original sequence IDs: {list(original_descriptions.keys())}")
                        self.log_message.emit(f"Aligned sequence IDs: {[rec.id for rec in alignment]}")
                        
                        for aligned_record in alignment:
                            # Try to match the ID to original sequences
                            seq_id = aligned_record.id
                            seq_desc = aligned_record.description if aligned_record.description else seq_id
                            base_id = seq_id.split()[0]
                            
                            # Find the original ID - try direct match first, then base ID
                            original_id = None
                            if seq_id in original_descriptions:
                                original_id = seq_id
                            elif seq_desc in original_id_map:
                                original_id = original_id_map[seq_desc]
                            elif base_id in original_id_map:
                                original_id = original_id_map[base_id]
                            else:
                                # Search for a partial match if direct match fails
                                for orig_id in original_descriptions:
                                    if orig_id.startswith(seq_id) or seq_id.startswith(orig_id):
                                        original_id = orig_id
                                        break
                                
                                # If still no match, try more fuzzy matching with descriptions
                                if not original_id:
                                    for desc in original_descriptions.values():
                                        if desc and seq_desc and (desc in seq_desc or seq_desc in desc):
                                            # Find the ID that corresponds to this description
                                            for oid, odesc in original_descriptions.items():
                                                if odesc == desc:
                                                    original_id = oid
                                                    break
                                            break
                            
                            # If we found a match, use the original description
                            if original_id:
                                description = original_descriptions[original_id]
                            else:
                                # Fallback to using the current description
                                description = aligned_record.description or seq_id
                                self.log_message.emit(f"Could not find original ID for '{seq_id}' - using current description")
                            
                            # Create a new record with the original description but aligned sequence
                            from Bio.SeqRecord import SeqRecord
                            
                            # Check for duplicate IDs in description (e.g., "HQ604103.1_ HQ604103.1_")
                            if original_id and description.startswith(original_id) and description.count(original_id) > 1:
                                # Remove duplicated IDs from description
                                parts = description.split()
                                seen = set()
                                unique_parts = []
                                for part in parts:
                                    if part not in seen:
                                        seen.add(part)
                                        unique_parts.append(part)
                                cleaned_description = " ".join(unique_parts)
                                self.log_message.emit(f"Fixed duplicate ID in description: {description} -> {cleaned_description}")
                                description = cleaned_description
                            
                            new_record = SeqRecord(
                                aligned_record.seq,
                                # Use full description as ID to preserve spaces in tree
                                id=description,
                                description=description
                            )
                            preserved_alignment.append(new_record)
                        
                        self.progress_update.emit(100)
                        self.alignment_complete.emit(preserved_alignment)
                    else:
                        self.log_message.emit(f"Failed to run {self.method_id}, falling back to simple alignment")
                        self._fallback_alignment()
                except Exception as e:
                    self.log_message.emit(f"Error running {self.method_id}: {str(e)}")
                    # Fallback to simple alignment if external tool fails
                    self._fallback_alignment()
            else:
                # Use simple alignment
                self.log_message.emit("Using simple alignment")
                self._fallback_alignment()
                
        except Exception as e:
            self.log_message.emit(f"Alignment thread error: {str(e)}")
            self.error_occurred.emit(str(e))
        finally:
            # Clean up temp files
            try:
                if os.path.exists(input_path):
                    os.remove(input_path)
                if os.path.exists(output_path):
                    os.remove(output_path)
            except Exception as e:
                self.log_message.emit(f"Error cleaning up temp files: {str(e)}")
                pass


    def _fallback_alignment(self):
        """Simple alignment method as fallback"""
        try:
            # Find the longest sequence
            max_len = max(len(record.seq) for record in self.sequences)

            # Create aligned sequences by padding with gaps
            aligned_sequences = []
            for i, record in enumerate(self.sequences):
                # Pad sequence to maximum length
                seq_str = str(record.seq)
                aligned_seq = seq_str + "-" * (max_len - len(seq_str))

                # Create new SeqRecord with aligned sequence - preserve full description
                from Bio.SeqRecord import SeqRecord
                from Bio.Seq import Seq
                
                # Get the description, avoiding duplicates
                if record.description and record.description != record.id:
                    # Check for duplicate IDs in the description
                    if record.description.startswith(record.id) and record.description.count(record.id) > 1:
                        # Remove duplicated IDs from description
                        parts = record.description.split()
                        seen = set()
                        unique_parts = []
                        for part in parts:
                            if part not in seen:
                                seen.add(part)
                                unique_parts.append(part)
                        description = " ".join(unique_parts)
                        self.log_message.emit(f"Fixed duplicate ID in description: {record.description} -> {description}")
                    else:
                        description = record.description
                else:
                    description = record.id
                
                aligned_record = SeqRecord(
                    Seq(aligned_seq),
                    id=description,  # Use cleaned description as ID to preserve spaces
                    description=description  # Keep cleaned description
                )
                aligned_sequences.append(aligned_record)

                # Update progress
                self.progress_update.emit(int((i+1) / len(self.sequences) * 80 + 20))

            self.progress_update.emit(100)
            self.alignment_complete.emit(aligned_sequences)
        except Exception as e:
            self.error_occurred.emit(f"Fallback alignment failed: {str(e)}")
    

class SequenceAlignmentDialog(QDialog):
    """Dialog for aligning sequences"""
    
    def __init__(self, parent=None, sequences=None):
        super().__init__(parent)
        self.setWindowTitle("Align Sequences")
        self.setMinimumWidth(450)
        self.setMinimumHeight(240)
        
        self.sequences = sequences or []
        self.aligned_sequences = None
        
        # Get available alignment methods
        self.alignment_methods = get_available_alignment_methods()
        self.method_map = {display_name: method_id for method_id, display_name in self.alignment_methods}
        
        # Create dialog layout and components
        self.create_layout()
        
        # Inherit dark mode setting from parent if available
        self.dark_mode = False
        if parent and hasattr(parent, 'dark_mode'):
            self.dark_mode = parent.dark_mode
            self.apply_theme()
        
        # Show information about missing tools if needed
        self.check_alignment_tools()
    
    def apply_theme(self):
        """Apply the current theme to the dialog"""
        pass
    
    def check_alignment_tools(self):
        """Check if alignment tools are installed and show warnings if needed"""
        # Check which tools are installed
        halign_installed = find_external_tool("halign") is not None
        muscle_installed = find_external_tool("muscle") is not None
        mafft_installed = find_external_tool("mafft") is not None
        
        # If halign is selected but not installed, show installation info
        selected_text = self.method_combo.currentText()
        if "halign" in selected_text and not halign_installed and "install" in selected_text:
            QMessageBox.information(
                self,
                "halign Installation",
                "The halign alignment program is recommended for best results but was not found "
                "in your PATH.\n\n"
                "To install halign:\n"
                "1. Download it from its source or use your package manager\n"
                "2. Make sure it's in your PATH environment variable\n"
                "3. The program should be executable with the 'halign' command\n\n"
                "Until halign is installed, TreeCraft will use alternative alignment methods."
            )
        
        # If MAFFT is selected but not installed, show installation info
        elif "MAFFT" in selected_text and not mafft_installed and "not found" in selected_text:
            QMessageBox.information(
                self,
                "MAFFT Installation",
                "MAFFT alignment program was not found in your PATH.\n\n"
                "To install MAFFT:\n"
                "1. Download it from https://mafft.cbrc.jp/alignment/software/ or use your package manager\n"
                "2. Make sure it's in your PATH environment variable\n"
                "3. The program should be executable with the 'mafft' command\n\n"
                "Until MAFFT is installed, TreeCraft will use alternative alignment methods."
            )
        
        # If no external tools are installed, warn the user
        if not halign_installed and not muscle_installed and not mafft_installed:
            QMessageBox.warning(
                self,
                "No Alignment Tools Found",
                "No external alignment tools (halign, MUSCLE, MAFFT) were found.\n\n"
                "TreeCraft will use a simple alignment algorithm, but results will be basic.\n"
                "Please consider installing halign, MUSCLE, or MAFFT for better alignments."
            )
    
    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Info label
        info_label = QLabel(f"Aligning {len(self.sequences)} sequences.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # Method selection
        method_layout = QHBoxLayout()
        method_label = QLabel("Alignment Method:")
        self.method_combo = QComboBox()
        
        # Add methods from available list
        for method_id, display_name in self.alignment_methods:
            self.method_combo.addItem(display_name)
        
        # Connect selection change signal
        self.method_combo.currentTextChanged.connect(self.on_method_changed)
        
        method_layout.addWidget(method_label)
        method_layout.addWidget(self.method_combo)
        layout.addLayout(method_layout)
        
        # Method description label
        self.method_desc_label = QLabel()
        self.method_desc_label.setWordWrap(True)
        layout.addWidget(self.method_desc_label)
        
        # Set initial description
        self.update_method_description(self.method_combo.currentText())
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        # Don't set initial value to avoid showing 0%
        self.progress_bar.setTextVisible(False)  # Initially hide text until alignment starts
        layout.addWidget(self.progress_bar)
        
        # Status label
        self.status_label = QLabel("Ready to align")
        layout.addWidget(self.status_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        
        self.align_button = QPushButton("Align")
        self.align_button.clicked.connect(self.start_alignment)
        self.align_button.setDefault(True)
        
        button_layout.addWidget(self.cancel_button)
        button_layout.addWidget(self.align_button)
        layout.addLayout(button_layout)
    
    def on_method_changed(self, method_text):
        """Handle alignment method selection change"""
        self.update_method_description(method_text)
    
    def update_method_description(self, method_text):
        """Update the description text based on selected method"""
        if "halign" in method_text:
            if "install" in method_text:
                self.method_desc_label.setText(
                    "halign is recommended for optimal alignments but is not installed. "
                    "Select a different method or install halign."
                )
            else:
                self.method_desc_label.setText(
                    "halign provides high-quality protein and DNA sequence alignments."
                )
        elif "MUSCLE" in method_text:
            if "not found" in method_text:
                self.method_desc_label.setText(
                    "MUSCLE alignment tool is not installed. "
                    "Please select a different method or install MUSCLE."
                )
            else:
                self.method_desc_label.setText(
                    "MUSCLE is a widely-used alignment algorithm with good accuracy."
                )
        elif "MAFFT" in method_text:
            if "not found" in method_text:
                self.method_desc_label.setText(
                    "MAFFT alignment tool is not installed. "
                    "Please select a different method or install MAFFT."
                )
            else:
                self.method_desc_label.setText(
                    "MAFFT is a fast, accurate multiple sequence alignment program "
                    "with --auto parameter for automatic algorithm selection."
                )
        else:
            self.method_desc_label.setText(
                "Simple alignment provides basic sequence alignment using gap insertion."
            )
    
    def start_alignment(self):
        """Start the alignment process"""
        # Disable controls
        self.align_button.setEnabled(False)
        self.method_combo.setEnabled(False)
        
        # Get selected method (internal ID from display name)
        selected_text = self.method_combo.currentText()
        
        # Extract method ID from the method map
        method_id = None
        for display, method in self.method_map.items():
            if display == selected_text:
                method_id = method
                break
        
        # Fallback to simple if no method was found
        if not method_id:
            method_id = "simple"
        
        logger.info(f"Starting alignment with method_id: {method_id}")
        
        # Update status
        method_name = method_id.upper()
        self.status_label.setText(f"Aligning sequences using {method_name}...")
        
        # Enable progress bar text now that alignment is starting
        self.progress_bar.setTextVisible(True)
        
        # Create and start the alignment thread
        self.align_thread = AlignmentThread(self.sequences, method_id)
        self.align_thread.progress_update.connect(self.update_progress)
        self.align_thread.alignment_complete.connect(self.alignment_complete)
        self.align_thread.error_occurred.connect(self.alignment_error)
        self.align_thread.log_message.connect(self.log_thread_message)
        self.align_thread.start()
    
    def update_progress(self, value):
        """Update the progress bar"""
        self.progress_bar.setValue(value)
    
    def alignment_complete(self, alignment):
        """Handle completed alignment"""
        self.aligned_sequences = alignment
        self.status_label.setText("Alignment complete!")
        
        # Auto-accept after successful alignment
        self.accept()
    
    def alignment_error(self, error_message):
        """Handle alignment error"""
        self.status_label.setText(f"Error: {error_message}")
        
        # Re-enable controls
        self.align_button.setEnabled(True)
        self.method_combo.setEnabled(True)
    
    def log_thread_message(self, message):
        """Safely log messages from the alignment thread"""
        # This runs in the main thread, so it's safe to call logger methods
        logger.info(message)
        
    def get_aligned_sequences(self):
        """Return the aligned sequences"""
        return self.aligned_sequences
