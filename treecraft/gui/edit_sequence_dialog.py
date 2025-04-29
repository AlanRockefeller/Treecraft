from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QTextEdit, QPushButton, QMessageBox
)
from PyQt6.QtCore import Qt

class EditSequenceDialog(QDialog):
    """Dialog for editing a sequence"""
    
    def __init__(self, parent=None, sequence_id="", sequence="", description=""):
        super().__init__(parent)
        self.setWindowTitle("Edit Sequence")
        self.setMinimumWidth(500)
        self.setMinimumHeight(300)
        
        self.original_id = sequence_id
        self.sequence_id = sequence_id
        self.sequence = sequence
        self.description = description
        
        self.create_layout()
        
        # Inherit dark mode setting from parent if available
        self.dark_mode = False
        if parent and hasattr(parent, 'dark_mode'):
            self.dark_mode = parent.dark_mode
            self.apply_theme()
    
    def apply_theme(self):
        """Apply the current theme to the dialog"""
        pass
    

    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Sequence ID
        id_layout = QHBoxLayout()
        id_label = QLabel("Sequence ID (internal identifier):")
        self.id_edit = QLineEdit(self.sequence_id)
        id_layout.addWidget(id_label)
        id_layout.addWidget(self.id_edit)
        layout.addLayout(id_layout)
        
        # Description
        desc_layout = QHBoxLayout()
        desc_label = QLabel("Display Name (shown in list):")
        self.desc_edit = QLineEdit(self.description)
        desc_layout.addWidget(desc_label)
        desc_layout.addWidget(self.desc_edit)
        layout.addLayout(desc_layout)
        
        # Add a note explaining the difference
        note = QLabel("Note: The Display Name is what appears in the sequence list.")
        note.setStyleSheet("color: blue;")
        layout.addWidget(note)
        
        # Sequence
        seq_label = QLabel("Sequence:")
        layout.addWidget(seq_label)
        
        self.seq_edit = QTextEdit()
        self.seq_edit.setText(self.sequence)
        layout.addWidget(self.seq_edit)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        
        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.validate_and_accept)
        self.save_button.setDefault(True)
        
        button_layout.addWidget(self.cancel_button)
        button_layout.addWidget(self.save_button)
        layout.addLayout(button_layout)


    def validate_and_accept(self):
        """Validate the edited sequence and accept the dialog"""
        # Get values from the UI controls
        new_id = self.id_edit.text().strip()
        new_description = self.desc_edit.text().strip()
        new_sequence = self.seq_edit.toPlainText().strip()
        
        print(f"Debug - validate_and_accept: ID={new_id}, Desc={new_description}")
        
        # Automatic update: If description was not changed but ID was, update the description to match the ID
        if new_description == self.description and new_id != self.original_id:
            new_description = new_id
            self.desc_edit.setText(new_description)
            print(f"Auto-updating description to match new ID: {new_description}")
        
        # Store the values
        self.sequence_id = new_id
        self.description = new_description  # Store the updated description
        self.sequence = ''.join(new_sequence.split())  # Remove whitespace
        
        # Validation
        if not self.sequence_id:
            QMessageBox.warning(self, "Missing ID", "Please provide a sequence ID.")
            return
        
        if not self.sequence:
            QMessageBox.warning(self, "Missing Sequence", "Please provide a sequence.")
            return
        
        # Check sequence characters
        valid_chars = set("ATGCURYKMSWBDHVN-")
        if not all(c.upper() in valid_chars for c in self.sequence):
            reply = QMessageBox.question(
                self, 
                "Invalid Sequence", 
                "Sequence contains non-standard characters. Save anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.No:
                return
        
        self.accept()
        

    def get_sequence_data(self):
        """Return the edited sequence data"""
        # Make sure we're returning the current text from the description field
        # This is crucial - it should return the EDITED description
        logger.debug(f"get_sequence_data returning: {self.original_id}, {self.sequence_id}, len={len(self.sequence)}, '{self.desc_edit.text()}'")
        return self.original_id, self.sequence_id, self.sequence, self.desc_edit.text()
