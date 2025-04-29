import os
from Bio import Phylo, SeqIO
from PyQt6.QtWidgets import QFileDialog, QMessageBox
from PyQt6.QtSvg import QSvgGenerator
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPainter

def open_fasta_file(parent=None):
    """Open a FASTA file and load sequences (*.fasta, *.fa, *.fas, *.fna, *.txt)"""
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open FASTA File", "", "FASTA Files (*.fasta *.fa *.fas *.fna *.txt);;All Files (*)"
    )
    
    if not file_path:
        return None, 0
        
    try:
        # Check for malformed FASTA files
        validation_results = validate_fasta_file(file_path)
        
        # Handle files with embedded headers
        if validation_results['embedded_headers']:
            problematic_lines = validation_results['embedded_headers']
            # Build a message about the problematic lines
            lines_msg = "\n".join([f"Line {line_num}: {line[:50]}..." for line_num, line in problematic_lines[:5]])
            if len(problematic_lines) > 5:
                lines_msg += f"\n(and {len(problematic_lines) - 5} more problematic lines)"
                
            # Ask user if they want to attempt to fix the file
            reply = QMessageBox.question(
                parent,
                "Malformed FASTA File",
                f"The file contains embedded sequence headers (> character in sequence lines).\n\n"
                f"This usually happens when FASTA files are concatenated incorrectly.\n\n"
                f"Problematic locations:\n{lines_msg}\n\n"
                f"Would you like to attempt to fix the file?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel
            )
            
            if reply == QMessageBox.StandardButton.Cancel:
                return None, 0
            elif reply == QMessageBox.StandardButton.Yes:
                # Try to fix the file
                try:
                    fixed_sequences = fix_malformed_fasta(file_path)
                    count = len(fixed_sequences)
                    QMessageBox.information(
                        parent, 
                        "FASTA File Fixed",
                        f"Successfully fixed and loaded {count} sequences."
                    )
                    return fixed_sequences, count
                except Exception as e:
                    QMessageBox.critical(
                        parent, 
                        "Error", 
                        f"Failed to fix FASTA file: {str(e)}"
                    )
                    return None, 0
                    
        # Handle files with invalid characters
        if validation_results['invalid_chars']:
            invalid_lines = validation_results['invalid_chars']
            # Build a message about the invalid characters
            lines_msg = "\n".join([f"Line {line_num}: {line[:40]}... (Invalid chars: '{invalid_chars}')" 
                                 for line_num, line, invalid_chars in invalid_lines[:5]])
            if len(invalid_lines) > 5:
                lines_msg += f"\n(and {len(invalid_lines) - 5} more lines with invalid characters)"
                
            # Notify user about invalid characters
            reply = QMessageBox.warning(
                parent,
                "Invalid Characters in FASTA File",
                f"The file contains invalid characters in sequence lines.\n\n"
                f"Valid characters are: A, T, G, C, U, R, Y, K, M, S, W, B, D, H, V, N, and -\n\n"
                f"Problematic locations:\n{lines_msg}\n\n"
                f"Would you like to load the file anyway? Invalid characters will be replaced with 'N'.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.No:
                return None, 0
        
        # Load sequences from FASTA
        sequences = list(SeqIO.parse(file_path, "fasta"))
        return sequences, len(sequences)
    except Exception as e:
        QMessageBox.critical(parent, "Error", f"Failed to load FASTA file: {str(e)}")
        return None, 0

def validate_fasta_file(file_path):
    """
    Check for malformed FASTA files with embedded sequence headers or invalid characters.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Dictionary with two lists:
            - 'embedded_headers': List of tuples (line_number, line_content) for lines with embedded headers
            - 'invalid_chars': List of tuples (line_number, line_content, invalid_chars) for lines with invalid characters
    """
    results = {
        'embedded_headers': [],
        'invalid_chars': []
    }
    
    # Define valid DNA/protein characters
    valid_chars = set("ATGCURYKMSWBDHVNatgcurykmswbdhvn-")
    
    with open(file_path, 'r') as f:
        in_header = False
        line_num = 0
        
        for line in f:
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
                # Check for embedded headers
                if '>' in line:
                    results['embedded_headers'].append((line_num, line))
                
                # Check for invalid characters
                invalid_chars = set()
                for char in line:
                    if char not in valid_chars:
                        invalid_chars.add(char)
                
                if invalid_chars:
                    results['invalid_chars'].append((line_num, line, ''.join(sorted(invalid_chars))))
    
    return results

def fix_malformed_fasta(file_path):
    """
    Fix a malformed FASTA file by splitting lines that contain embedded headers.
    
    Args:
        file_path: Path to the malformed FASTA file
        
    Returns:
        List of Bio.SeqRecord objects from the fixed file
    """
    from io import StringIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    fixed_records = []
    current_header = None
    current_seq = ""
    
    with open(file_path, 'r') as f:
        for line in f:
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
                        
                    # Split into header and any remaining sequence
                    header_parts = part.split('\n', 1)
                    
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
        
def open_fastq_file(parent=None):
    """Open a FASTQ file and load sequences"""
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open FASTQ File", "", "FASTQ Files (*.fastq *.fq);;All Files (*)"
    )
    
    if not file_path:
        return None, 0
        
    try:
        # Load sequences from FASTQ
        sequences = list(SeqIO.parse(file_path, "fastq"))
        return sequences, len(sequences)
    except Exception as e:
        QMessageBox.critical(parent, "Error", f"Failed to load FASTQ file: {str(e)}")
        return None, 0

def open_tree_file(parent=None):
    """Open a tree file (Newick, etc.)"""
    file_path, _ = QFileDialog.getOpenFileName(
        parent, "Open Tree File", "", "Newick Files (*.nwk *.newick);;Nexus Files (*.nex);;PhyloXML Files (*.xml);;All Files (*)"
    )
    
    if not file_path:
        return None
    
    # Determine format based on file extension
    ext = os.path.splitext(file_path)[1].lower()
    if ext in (".nex", ".nexus"):
        format_name = "nexus"
    elif ext in (".xml", ".phyloxml"):
        format_name = "phyloxml"
    else:
        format_name = "newick"  # Default
        
    try:
        # Read the tree file
        tree = Phylo.read(file_path, format_name)
        
        # Process the tree to handle quoted names
        for clade in tree.get_terminals() + list(tree.get_nonterminals()):
            if clade.name:
                # Remove quotes from quoted names for display
                if clade.name.startswith("'") and clade.name.endswith("'"):
                    clade.name = clade.name[1:-1]
                # Also support double quotes
                elif clade.name.startswith('"') and clade.name.endswith('"'):
                    clade.name = clade.name[1:-1]
                    
        return tree
    except Exception as e:
        QMessageBox.critical(parent, "Error", f"Failed to load tree file: {str(e)}")
        return None


def save_tree_file(parent=None, tree=None, file_path=None, format_name=None):
    """Save the current tree to a file"""
    if not tree:
        if parent:
            QMessageBox.warning(parent, "Warning", "No tree to save")
        return False
    
    if not file_path:
        # Define supported formats
        formats = "Newick (*.nwk);;Nexus (*.nex);;PhyloXML (*.xml);;All Files (*)"
        
        # Open save dialog
        file_path, selected_filter = QFileDialog.getSaveFileName(
            parent, "Save Tree", "", formats
        )
        
        if not file_path:
            return False
        
        # Determine format based on selected filter
        if not format_name:
            if "Nexus" in selected_filter:
                format_name = "nexus"
            elif "PhyloXML" in selected_filter:
                format_name = "phyloxml"
            else:
                format_name = "newick"  # Default
    
    # If format is still not specified, try to determine from file extension
    if not format_name:
        ext = os.path.splitext(file_path)[1].lower()
        if ext in (".nex", ".nexus"):
            format_name = "nexus"
        elif ext in (".xml", ".phyloxml"):
            format_name = "phyloxml"
        else:
            format_name = "newick"  # Default
    
    try:
        from Bio import Phylo
        
        # Ensure all node names with spaces or parentheses are properly quoted for Newick format
        if format_name == "newick":
            # Process all nodes to ensure names with spaces or parentheses are properly quoted
            for clade in tree.get_terminals() + list(tree.get_nonterminals()):
                if clade.name and (
                    " " in clade.name or 
                    "(" in clade.name or 
                    ")" in clade.name
                ) and not (clade.name.startswith("'") and clade.name.endswith("'")):
                    # Quote the name to preserve spaces and parentheses in Newick format
                    clade.name = f"'{clade.name}'"
                    
        # Write the tree to file
        Phylo.write(tree, file_path, format_name)
        return True
    except Exception as e:
        if parent:
            QMessageBox.critical(parent, "Error", f"Failed to save tree: {str(e)}")
        return False



def export_image(parent=None, tree_canvas=None):
    """Export the tree visualization as an image"""
    if not tree_canvas or not tree_canvas.tree:
        QMessageBox.warning(parent, "Warning", "No tree to export")
        return False
    
    file_path, selected_filter = QFileDialog.getSaveFileName(
        parent, "Export Image", "", 
        "PNG Image (*.png);;JPEG Image (*.jpg);;SVG Image (*.svg)"
    )
    
    if not file_path:
        return False
        
    try:
        # Create a pixmap of the tree canvas
        pixmap = tree_canvas.grab()
        
        # Determine format based on file extension
        if file_path.lower().endswith('.png'):
            pixmap.save(file_path, 'PNG')
        elif file_path.lower().endswith('.jpg') or file_path.lower().endswith('.jpeg'):
            pixmap.save(file_path, 'JPG', quality=90)
        elif file_path.lower().endswith('.svg'):
            # For SVG we need a different approach
            generator = QSvgGenerator()
            generator.setFileName(file_path)
            generator.setSize(tree_canvas.size())
            generator.setViewBox(tree_canvas.rect())
            
            painter = QPainter()
            painter.begin(generator)
            tree_canvas.render(painter)
            painter.end()
        
        return True
    except Exception as e:
        QMessageBox.critical(parent, "Error", f"Failed to export image: {str(e)}")
        return False
