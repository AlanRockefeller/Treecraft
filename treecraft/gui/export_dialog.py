from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QComboBox, QCheckBox, QGroupBox, QFormLayout,
    QPushButton, QFileDialog
)
from PyQt6.QtCore import Qt

class ExportDialog(QDialog):
    """Dialog for exporting trees to various formats"""
    
    def __init__(self, parent=None, tree=None):
        super().__init__(parent)
        self.tree = tree
        self.setWindowTitle("Export Tree")
        self.setMinimumWidth(350)
        
        # Initialize layout
        self.create_layout()
        
    def create_layout(self):
        """Create dialog layout with export options"""
        layout = QVBoxLayout(self)
        
        # Format selection
        format_group = QGroupBox("Export Format")
        format_layout = QVBoxLayout()
        
        self.format_combo = QComboBox()
        self.format_combo.addItems(["Newick (.nwk)", "Nexus (.nex)", "PhyloXML (.xml)", 
                                    "PNG Image (.png)", "JPEG Image (.jpg)", "SVG Image (.svg)"])
        format_layout.addWidget(self.format_combo)
        format_group.setLayout(format_layout)
        layout.addWidget(format_group)
        
        # Connect format change to handle different options
        self.format_combo.currentIndexChanged.connect(self.on_format_changed)
        
        # Image export options info
        self.image_options_info = QLabel("Configure image options when exporting")
        self.image_options_info.setVisible(False)
        layout.addWidget(self.image_options_info)
        
        # Buttons
        button_layout = QHBoxLayout()
        
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.cancel_button)
        
        self.export_button = QPushButton("Export")
        self.export_button.clicked.connect(self.accept)
        self.export_button.setDefault(True)
        button_layout.addWidget(self.export_button)
        
        layout.addLayout(button_layout)
        
    def on_format_changed(self, index):
        """Handle format selection change"""
        format_text = self.format_combo.currentText()
        # Show image options info for image formats
        is_image = any(ext in format_text for ext in [".png", ".jpg", ".svg"])
        self.image_options_info.setVisible(is_image)
        
    def get_selected_format(self):
        """Get the selected export format"""
        format_text = self.format_combo.currentText()
        # Extract extension from the format string
        extension = format_text.split("(")[1].split(")")[0].strip(".")
        return extension

class ImageExportOptionsDialog(QDialog):
    """Dialog for configuring image export options"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Image Export Options")
        self.setMinimumWidth(350)
        
        # Initialize default options
        self.background_color = "current"  # current, black, or white
        self.high_resolution = False
        self.ultra_resolution = False
        
        self.create_layout()
    
    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Background color options
        bg_group = QGroupBox("Background Color")
        bg_layout = QVBoxLayout()
        
        self.bg_current_radio = QCheckBox("Use current theme")
        self.bg_current_radio.setChecked(True)
        bg_layout.addWidget(self.bg_current_radio)
        
        self.bg_black_radio = QCheckBox("Force black background with white text")
        bg_layout.addWidget(self.bg_black_radio)
        
        self.bg_white_radio = QCheckBox("Force white background with black text")
        bg_layout.addWidget(self.bg_white_radio)
        
        # Connect radio buttons to ensure only one is checked
        self.bg_current_radio.toggled.connect(self.bg_radio_toggled)
        self.bg_black_radio.toggled.connect(self.bg_radio_toggled)
        self.bg_white_radio.toggled.connect(self.bg_radio_toggled)
        
        bg_group.setLayout(bg_layout)
        layout.addWidget(bg_group)
        
        # Resolution options
        res_group = QGroupBox("Resolution Options")
        res_layout = QVBoxLayout()
        
        self.high_res_check = QCheckBox("High Resolution (2x)")
        res_layout.addWidget(self.high_res_check)
        
        self.ultra_res_check = QCheckBox("Ultra Resolution (4x)")
        res_layout.addWidget(self.ultra_res_check)
        
        # Make resolution options mutually exclusive
        self.high_res_check.toggled.connect(self.res_check_toggled)
        self.ultra_res_check.toggled.connect(self.res_check_toggled)
        
        res_group.setLayout(res_layout)
        layout.addWidget(res_group)
        
        # Info label
        info_label = QLabel("Note: SVG exports are always high quality regardless of resolution settings.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # Buttons
        button_layout = QHBoxLayout()
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        ok_button.setDefault(True)
        
        button_layout.addWidget(cancel_button)
        button_layout.addWidget(ok_button)
        layout.addLayout(button_layout)
    
    def bg_radio_toggled(self, checked):
        """Handle background radio button toggling"""
        if checked:
            sender = self.sender()
            if sender == self.bg_current_radio:
                self.bg_black_radio.setChecked(False)
                self.bg_white_radio.setChecked(False)
                self.background_color = "current"
            elif sender == self.bg_black_radio:
                self.bg_current_radio.setChecked(False)
                self.bg_white_radio.setChecked(False)
                self.background_color = "black"
            elif sender == self.bg_white_radio:
                self.bg_current_radio.setChecked(False)
                self.bg_black_radio.setChecked(False)
                self.background_color = "white"
    
    def res_check_toggled(self, checked):
        """Handle resolution checkbox toggling"""
        if checked:
            sender = self.sender()
            if sender == self.high_res_check:
                self.ultra_res_check.setChecked(False)
                self.high_resolution = True
                self.ultra_resolution = False
            elif sender == self.ultra_res_check:
                self.high_res_check.setChecked(False)
                self.high_resolution = False
                self.ultra_resolution = True
        else:
            # If unchecked and it was the active one, reset both flags
            sender = self.sender()
            if (sender == self.high_res_check and self.high_resolution) or \
               (sender == self.ultra_res_check and self.ultra_resolution):
                self.high_resolution = False
                self.ultra_resolution = False
    
    def get_options(self):
        """Get the selected export options"""
        return {
            "background_color": self.background_color,
            "high_resolution": self.high_resolution,
            "ultra_resolution": self.ultra_resolution
        }
