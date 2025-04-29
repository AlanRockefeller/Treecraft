# treecraft/gui/tree_dialog.py
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QComboBox, QPushButton, QGroupBox, QRadioButton,
    QSpinBox, QFormLayout, QCheckBox
)
from PyQt6.QtCore import Qt

class TreeBuildingDialog(QDialog):
    """Dialog for selecting tree-building parameters"""
    
    def __init__(self, parent=None, sequence_count=0):
        super().__init__(parent)
        self.setWindowTitle("Build Phylogenetic Tree")
        self.setMinimumWidth(400)
        
        self.sequence_count = sequence_count
        self.create_layout()
        
        # Inherit dark mode setting from parent if available
        self.dark_mode = False
        if parent and hasattr(parent, 'dark_mode'):
            self.dark_mode = parent.dark_mode
            self.apply_theme()
    
    def apply_theme(self):
        """Apply the current theme to the dialog"""
        # Theme-specific styling could be added here
        pass
    


    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Method selection
        method_group = QGroupBox("Tree Building Method")
        method_layout = QVBoxLayout()
        
        self.upgma_radio = QRadioButton("UPGMA (Fast, ultrametric tree)")
        self.nj_radio = QRadioButton("Neighbor Joining (Fast, non-ultrametric tree)")
        self.ml_radio = QRadioButton("Maximum Likelihood (Slow, more accurate)")
        
        # Select UPGMA by default
        self.upgma_radio.setChecked(True)
        
        # Disable ML if too few sequences
        if self.sequence_count < 4:
            self.ml_radio.setEnabled(False)
            self.ml_radio.setToolTip("Requires at least 4 sequences")
        
        method_layout.addWidget(self.upgma_radio)
        method_layout.addWidget(self.nj_radio)
        method_layout.addWidget(self.ml_radio)
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)
        
        # Distance model selection
        distance_group = QGroupBox("Distance Model")
        distance_layout = QFormLayout()
        
        self.distance_model = QComboBox()
        self.distance_model.addItems(["Identity", "BLOSUM62", "PAM250"])
        distance_layout.addRow("Distance Model:", self.distance_model)
        
        distance_group.setLayout(distance_layout)
        layout.addWidget(distance_group)
        
        # Add bootstrap options
        bootstrap_group = QGroupBox("Bootstrap Analysis")
        bootstrap_layout = QVBoxLayout()
        
        self.bootstrap_check = QCheckBox("Perform bootstrap analysis")
        bootstrap_layout.addWidget(self.bootstrap_check)
        
        replicates_layout = QHBoxLayout()
        replicates_layout.addWidget(QLabel("Number of replicates:"))
        self.bootstrap_replicates = QSpinBox()
        self.bootstrap_replicates.setRange(10, 1000)
        self.bootstrap_replicates.setValue(100)
        self.bootstrap_replicates.setEnabled(False)
        replicates_layout.addWidget(self.bootstrap_replicates)
        
        bootstrap_layout.addLayout(replicates_layout)
        
        # Connect bootstrap checkbox to enable/disable replicates spinbox
        self.bootstrap_check.toggled.connect(self.bootstrap_replicates.setEnabled)
        
        bootstrap_group.setLayout(bootstrap_layout)
        layout.addWidget(bootstrap_group)
        
        # Advanced options
        advanced_group = QGroupBox("Advanced Options")
        advanced_layout = QFormLayout()
        
        advanced_group.setLayout(advanced_layout)
        layout.addWidget(advanced_group)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        
        self.build_button = QPushButton("Build Tree")
        self.build_button.clicked.connect(self.accept)
        self.build_button.setDefault(True)
        
        button_layout.addWidget(self.cancel_button)
        button_layout.addWidget(self.build_button)
        layout.addLayout(button_layout)

    def get_parameters(self):
        """Get the selected tree-building parameters"""
        method = "upgma"
        if self.nj_radio.isChecked():
            method = "nj"
        elif self.ml_radio.isChecked():
            method = "ml"
        
        distance_model = self.distance_model.currentText()
        
        bootstrap = self.bootstrap_check.isChecked()
        bootstrap_replicates = self.bootstrap_replicates.value() if bootstrap else 0
        
        return {
            "method": method,
            "distance_model": distance_model,
            "bootstrap": bootstrap,
            "bootstrap_replicates": bootstrap_replicates
        }
