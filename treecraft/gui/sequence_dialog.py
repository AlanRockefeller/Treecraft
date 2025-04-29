from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QTextEdit, QPushButton, QMessageBox
)
from PyQt6.QtGui import QPalette, QColor
from Bio import SeqIO
from io import StringIO

class AddSequenceDialog(QDialog):
    """Dialog for adding new sequences"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Sequence")
        self.setMinimumWidth(500)
        self.setMinimumHeight(300)
        
        self.sequence_id = ""
        self.sequence = ""
        
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Sequence")
        self.setMinimumWidth(500)
        self.setMinimumHeight(300)

        # Make the dialog non-modal
        self.setModal(False)

        self.sequence_id = ""
        self.sequence = ""

        # Inherit dark mode setting from parent if available
        self.dark_mode = False
        if parent and hasattr(parent, 'dark_mode'):
            self.dark_mode = parent.dark_mode
            self.apply_theme()

        self.create_layout()
    
    def apply_theme(self):
        """Apply the current theme to the dialog"""
        if self.dark_mode:
            # Set placeholder text color for dark mode
            palette = self.palette()
            placeholder_color = QColor(180, 180, 180)  # Lighter gray for dark mode
            palette.setColor(QPalette.ColorRole.PlaceholderText, placeholder_color)
            self.setPalette(palette)
    
    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Sequence ID
        id_layout = QHBoxLayout()
        id_label = QLabel("Sequence ID:")
        self.id_edit = QLineEdit()
        id_layout.addWidget(id_label)
        id_layout.addWidget(self.id_edit)
        layout.addLayout(id_layout)
        
        # Sequence
        seq_label = QLabel("Sequence (FASTA/FASTQ format accepted, including .txt files):")
        layout.addWidget(seq_label)
        
        self.seq_edit = QTextEdit()
        self.seq_edit.setPlaceholderText("Enter sequence or paste FASTA/FASTQ content...\n\nFASTA Example:\n>Sequence1\nATGCATGCATGC...\n\nFASTQ Example:\n@Sequence1\nATGCATGCATGC...\n+\nIIIIIIIIIIII...")
        layout.addWidget(self.seq_edit)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        
        self.add_button = QPushButton("Add")
        self.add_button.clicked.connect(self.validate_and_accept)
        self.add_button.setDefault(True)
        
        button_layout.addWidget(self.cancel_button)
        button_layout.addWidget(self.add_button)
        layout.addLayout(button_layout)
    
    def validate_and_accept(self):
        """Validate the sequence and accept the dialog"""
        text = self.seq_edit.toPlainText().strip()
        
        # Check if it's a FASTA format
        if text.startswith(">"):
            self.handle_fasta_input(text)
        # Check if it's a FASTQ format
        elif text.startswith("@"):
            self.handle_fastq_input(text)
        else:
            self.handle_raw_sequence_input(text)
    
    def handle_fasta_input(self, text):
        """Process FASTA format input"""
        # Check for malformed FASTA with embedded headers
        has_embedded_headers = self.check_for_embedded_headers(text)
        invalid_characters = self.check_for_invalid_characters(text)
        
        # Handle embedded headers first
        if has_embedded_headers:
            reply = QMessageBox.question(
                self,
                "Malformed FASTA",
                "The pasted FASTA content contains '>' characters in sequence lines, "
                "which indicates embedded sequence headers.\n\n"
                "This usually happens when FASTA files are concatenated incorrectly.\n\n"
                "Would you like to attempt to fix this?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel
            )
            
            if reply == QMessageBox.StandardButton.Cancel:
                return
            elif reply == QMessageBox.StandardButton.Yes:
                # Attempt to fix the FASTA
                try:
                    records = self.fix_malformed_fasta_text(text)
                    if not records:
                        QMessageBox.warning(self, "Fix Failed", "No valid sequences found after fixing.")
                        return
                        
                    # Use the first sequence
                    self.sequence_id = records[0].id
                    self.sequence = str(records[0].seq)
                    
                    # Update the ID field
                    self.id_edit.setText(self.sequence_id)
                    
                    # Inform user about all sequences
                    if len(records) > 1:
                        QMessageBox.information(
                            self,
                            "Multiple Sequences",
                            f"Fixed FASTA content contains {len(records)} sequences. "
                            f"Only the first sequence has been imported. "
                            f"Use 'File > Import FASTA' to import all sequences."
                        )
                    
                    self.accept()
                    return
                except Exception as e:
                    QMessageBox.warning(self, "Fix Failed", f"Could not fix malformed FASTA: {str(e)}")
                    return
        
        # Handle invalid characters
        if invalid_characters:
            # Build a message about the invalid characters
            lines_msg = "\n".join([f"Line {line_num}: {line[:40]}... (Invalid chars: '{invalid_chars}')" 
                                for line_num, line, invalid_chars in invalid_characters[:5]])
            if len(invalid_characters) > 5:
                lines_msg += f"\n(and {len(invalid_characters) - 5} more lines with invalid characters)"
                
            # Notify user about invalid characters
            reply = QMessageBox.warning(
                self,
                "Invalid Characters in FASTA Sequence",
                f"The FASTA content contains invalid characters in sequence lines.\n\n"
                f"Valid characters are: A, T, G, C, U, R, Y, K, M, S, W, B, D, H, V, N, and -\n\n"
                f"Problematic locations:\n{lines_msg}\n\n"
                f"Would you like to continue anyway? Invalid characters will be replaced with 'N'.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.No:
                return
        
        # Parse FASTA
        try:
            fasta_io = StringIO(text)
            records = list(SeqIO.parse(fasta_io, "fasta"))
            
            if not records:
                QMessageBox.warning(self, "Invalid FASTA", "No valid sequences found in FASTA input.")
                return
            
            # Use the first sequence
            self.sequence_id = records[0].id
            self.sequence = str(records[0].seq)
            
            # Update the ID field
            self.id_edit.setText(self.sequence_id)
            self.accept()
            
        except Exception as e:
            QMessageBox.warning(self, "Invalid FASTA", f"Could not parse FASTA: {str(e)}")
            return
    
    def check_for_embedded_headers(self, text):
        """
        Check if FASTA text has embedded sequence headers (> in sequence lines)
        
        Args:
            text: FASTA format text
            
        Returns:
            Boolean indicating if embedded headers were found
        """
        lines = text.strip().split('\n')
        in_header = False
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                in_header = True
                continue
                
            if in_header:
                in_header = False
                
            # Check if a sequence line contains '>' character (embedded header)
            if not line.startswith('>') and '>' in line:
                return True
                
        return False
        
    def check_for_invalid_characters(self, text):
        """
        Check if FASTA text has invalid characters in sequence lines
        
        Args:
            text: FASTA format text
            
        Returns:
            List of tuples (line_number, line_content, invalid_chars) for problematic lines,
            or empty list if no issues found
        """
        lines = text.strip().split('\n')
        in_header = False
        line_num = 0
        problematic_lines = []
        
        # Define valid DNA/protein characters
        valid_chars = set("ATGCURYKMSWBDHVNatgcurykmswbdhvn-")
        
        for line in lines:
            line_num += 1
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # If this is a header line
            if line.startswith('>'):
                in_header = True
                continue
                
            # If this is a sequence line
            if in_header:
                in_header = False
            
            # Only check sequence lines (not headers)
            if not line.startswith('>'):
                # Check for invalid characters
                invalid_chars = set()
                for char in line:
                    if char not in valid_chars:
                        invalid_chars.add(char)
                
                if invalid_chars:
                    problematic_lines.append((line_num, line, ''.join(sorted(invalid_chars))))
        
        return problematic_lines
        
    def fix_malformed_fasta_text(self, text):
        """
        Fix malformed FASTA text with embedded headers
        
        Args:
            text: Malformed FASTA text
            
        Returns:
            List of Bio.SeqRecord objects
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        fixed_records = []
        current_header = None
        current_seq = ""
        
        lines = text.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # If this is a header line
            if line.startswith('>'):
                # Save previous record if exists
                if current_header and current_seq:
                    fixed_records.append(
                        SeqRecord(Seq(current_seq), id=current_header, description=current_header)
                    )
                
                # Start new record
                current_header = line[1:].strip()
                current_seq = ""
            # If this is a sequence line that might contain an embedded header
            elif '>' in line:
                # Split the line at the embedded header
                parts = line.split('>')
                
                # Add the first part to the current sequence
                current_seq += parts[0]
                
                # Save the current record
                if current_header and current_seq:
                    fixed_records.append(
                        SeqRecord(Seq(current_seq), id=current_header, description=current_header)
                    )
                
                # Process the remaining parts as new records
                for i, part in enumerate(parts[1:]):
                    if not part:
                        continue
                        
                    # Split into header and any remaining sequence on this line
                    header_parts = part.split(None, 1)
                    
                    if len(header_parts) == 1:
                        # Just a header with no sequence on this line
                        new_header = header_parts[0].strip()
                        new_seq = ""
                    else:
                        # Header and some sequence
                        new_header = header_parts[0].strip()
                        new_seq = header_parts[1].strip()
                    
                    # If this is the last part and it has no sequence yet,
                    # set up as current record for next lines
                    if i == len(parts) - 2 and not new_seq:
                        current_header = new_header
                        current_seq = ""
                    else:
                        # Otherwise add as a complete record
                        fixed_records.append(
                            SeqRecord(Seq(new_seq), id=new_header, description=new_header)
                        )
                        # Reset if this was the last part
                        if i == len(parts) - 2:
                            current_header = None
                            current_seq = ""
            # Regular sequence line
            else:
                current_seq += line
        
        # Add the last record if there is one
        if current_header and current_seq:
            fixed_records.append(
                SeqRecord(Seq(current_seq), id=current_header, description=current_header)
            )
        
        return fixed_records
        
    def handle_fastq_input(self, text):
        """Process FASTQ format input"""
        # Parse FASTQ
        try:
            fastq_io = StringIO(text)
            records = list(SeqIO.parse(fastq_io, "fastq"))
            
            if not records:
                QMessageBox.warning(self, "Invalid FASTQ", "No valid sequences found in FASTQ input.")
                return
            
            # Use the first sequence
            self.sequence_id = records[0].id
            self.sequence = str(records[0].seq)
            
            # Update the ID field
            self.id_edit.setText(self.sequence_id)
            self.accept()
            
        except Exception as e:
            QMessageBox.warning(self, "Invalid FASTQ", f"Could not parse FASTQ: {str(e)}")
            return
    
    def handle_raw_sequence_input(self, text):
        """Process raw sequence input"""
        self.sequence_id = self.id_edit.text().strip()
        self.sequence = ''.join(text.split())  # Remove whitespace
        
        if not self.sequence_id:
            QMessageBox.warning(self, "Missing ID", "Please provide a sequence ID.")
            return
        
        if not self.sequence:
            QMessageBox.warning(self, "Missing Sequence", "Please provide a sequence.")
            return
        
        # Validate sequence characters (basic check)
        valid_chars = set("ATGCURYKMSWBDHVN-")
        if not all(c.upper() in valid_chars for c in self.sequence):
            reply = QMessageBox.question(
                self, 
                "Invalid Sequence", 
                "Sequence contains non-standard characters. Add anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.No:
                return
        
        self.accept()
    
    def get_sequence_data(self):
        """Return the sequence data"""
        return self.sequence_id, self.sequence
