from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
    QComboBox, QPushButton, QGroupBox, QRadioButton,
    QSpinBox, QFormLayout, QCheckBox, QApplication,
    QTabWidget, QWidget, QLineEdit, QMessageBox,
    QProgressDialog, QTextEdit
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QTextCharFormat, QColor
import logging
import os
import tempfile
import subprocess
from treecraft.utils.external_tools import find_external_tool

class RaxMLDialog(QDialog):
    """Dialog for RAxML tree-building parameters"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("RAxML Tree Building")
        self.setMinimumWidth(500)
        self.setMinimumHeight(400)
        
        # Initialize state variables
        self.find_model = False
        self.model_tool = "modeltest-ng"
        self.trim_tool = "trimal"
        
        # Check for availability of external tools
        logger = logging.getLogger("treecraft.raxml_dialog")
        
        # Check for all the needed external tools with more detailed logging
        logger.info("Checking for external tool availability...")
        
        # Check for Gblocks with full details
        logger.info("Checking for Gblocks...")
        gblocks_path = find_external_tool("gblocks")
        logger.info(f"Gblocks detection result: {gblocks_path}")
        
        # Check for other tools
        trimal_path = find_external_tool("trimal")
        logger.info(f"trimAl detection result: {trimal_path}")
        
        modeltest_path = find_external_tool("modeltest-ng")
        logger.info(f"ModelTest-NG detection result: {modeltest_path}")
        
        iqtree_path = find_external_tool("iqtree")
        logger.info(f"IQ-TREE detection result: {iqtree_path}")
        
        # Store availability flags
        self.has_trimal = trimal_path is not None
        self.has_gblocks = gblocks_path is not None
        self.has_modeltest = modeltest_path is not None
        self.has_iqtree = iqtree_path is not None
        
        # Store actual paths for diagnostic purposes
        self.gblocks_path = gblocks_path
        self.trimal_path = trimal_path
        
        # Log tool availability for debugging
        logger.info(f"Tool availability summary:")
        logger.info(f"  trimAl: {self.has_trimal} (Path: {trimal_path})")
        logger.info(f"  Gblocks: {self.has_gblocks} (Path: {gblocks_path})")
        logger.info(f"  ModelTest-NG: {self.has_modeltest} (Path: {modeltest_path})")
        logger.info(f"  IQ-TREE: {self.has_iqtree} (Path: {iqtree_path})")
        
        # Inherit dark mode setting from parent if available
        self.dark_mode = False
        if parent and hasattr(parent, 'dark_mode'):
            self.dark_mode = parent.dark_mode
            self.apply_theme()
        
        self.create_layout()
        
        # Make a final verification check for Gblocks, but don't use the -h flag which can hang
        if self.has_gblocks:
            logger.info("Performing a verification check for Gblocks using the diagnosed path")
            # If we have a gblocks_path, verify it exists and is executable - don't try to run it
            if self.gblocks_path:
                import os
                if os.path.isfile(self.gblocks_path) and os.access(self.gblocks_path, os.X_OK):
                    logger.info(f"Verified Gblocks at {self.gblocks_path} exists and is executable")
                    # This is sufficient verification - don't try to run the command which may hang
                    self.has_gblocks = True
                else:
                    logger.warning(f"Gblocks at {self.gblocks_path} is not a valid executable")
                    self.has_gblocks = False
        
        # Apply initial styles to indicate default tools
        self.trimal_button.setStyleSheet("font-weight: bold;" if self.trim_tool == "trimal" else "font-weight: normal;")
        self.gblocks_button.setStyleSheet("font-weight: bold;" if self.trim_tool == "gblocks" else "font-weight: normal;")
        
        # Set proper tooltips for all tools with more detailed information
        self.trimal_button.setToolTip("trimAl: " + ("Available" if self.has_trimal else "Not found in PATH") + 
                                  (f" ({self.trimal_path})" if self.has_trimal and hasattr(self, 'trimal_path') else ""))
        
        self.gblocks_button.setToolTip("Gblocks: " + ("Available" if self.has_gblocks else "Not found in PATH") +
                                  (f" ({self.gblocks_path})" if self.has_gblocks and hasattr(self, 'gblocks_path') else ""))
        
        self.findmodel_button.setToolTip("ModelTest-NG: " + ("Available" if self.has_modeltest else "Not found in PATH"))
        self.iqtree_button.setToolTip("IQ-TREE: " + ("Available" if self.has_iqtree else "Not found in PATH"))
        
        # Disable model finding buttons if not available
        if not self.has_modeltest:
            self.findmodel_button.setEnabled(False)
            
        if not self.has_iqtree:
            self.iqtree_button.setEnabled(False)
            
        # If neither trim tool is available, disable the trim checkbox
        if not self.has_trimal and not self.has_gblocks:
            self.trim_check.setEnabled(False)
            self.trim_check.setToolTip("No trimming tools available (install trimAl or Gblocks)")
            logger.warning("Both trimming tools are unavailable, disabling trim checkbox")
            
        # Select appropriate default trim tool based on availability
        if not self.has_trimal and self.has_gblocks:
            self.trim_tool = "gblocks"
            self.trimal_button.setStyleSheet("font-weight: normal;")
            self.gblocks_button.setStyleSheet("font-weight: bold;")
            logger.info("trimAl unavailable but Gblocks is available - selecting Gblocks as default")
        
        # Log the final tool selection
        logger.info(f"Final trim tool selection: {self.trim_tool}")
            
        # Now initialize button states - this will properly enable/disable based on checkbox and availability
        self.toggle_trim_tools(self.trim_check.isChecked())
    
    def apply_theme(self):
        """Apply the current theme to the dialog"""
        # Theme-specific styling could be added here
        pass
    
    def create_layout(self):
        """Create the dialog layout"""
        layout = QVBoxLayout(self)
        
        # Create tabs for different parameter groups
        tabs = QTabWidget()
        
        # Basic tab
        basic_tab = QWidget()
        basic_layout = QVBoxLayout(basic_tab)
        
        # Model selection
        model_group = QGroupBox("Substitution Model")
        model_layout = QVBoxLayout()
        
        self.model_combobox = QComboBox()
        self.model_combobox.addItems([
            "GTR+GAMMA", "GTR+GAMMA+I", "GTR+CAT", "GTR+CAT+I",
            "HKY+GAMMA", "HKY+GAMMA+I", "JC+GAMMA", "WAG+GAMMA",
            "LG+GAMMA", "BLOSUM62"
        ])
        
        model_auto_layout = QHBoxLayout()
        self.findmodel_button = QPushButton("Find Best Model with ModelTest-NG")
        self.iqtree_button = QPushButton("Find Model with IQ-TREE (slow)")
        model_auto_layout.addWidget(self.findmodel_button)
        model_auto_layout.addWidget(self.iqtree_button)
        
        # Add tooltip for the IQ-TREE button
        self.iqtree_button.setToolTip("Runs IQ-TREE immediately to find the best evolutionary model for your data")
        
        # Connect model finder buttons
        self.findmodel_button.clicked.connect(self.on_modeltest_clicked)
        self.iqtree_button.clicked.connect(self.on_iqtree_clicked)
        
        model_layout.addWidget(QLabel("Select Evolutionary Model:"))
        model_layout.addWidget(self.model_combobox)
        model_layout.addLayout(model_auto_layout)
        model_group.setLayout(model_layout)
        basic_layout.addWidget(model_group)
        
        # Bootstrap options
        bootstrap_group = QGroupBox("Bootstrap Analysis")
        bootstrap_layout = QVBoxLayout()
        
        self.bootstrap_check = QCheckBox("Perform bootstrap analysis")
        bootstrap_layout.addWidget(self.bootstrap_check)
        
        bootstrap_replicates_layout = QHBoxLayout()
        bootstrap_replicates_layout.addWidget(QLabel("Number of replicates:"))
        self.bootstrap_replicates = QSpinBox()
        self.bootstrap_replicates.setRange(10, 1000)
        self.bootstrap_replicates.setValue(100)
        self.bootstrap_replicates.setEnabled(False)
        bootstrap_replicates_layout.addWidget(self.bootstrap_replicates)
        bootstrap_layout.addLayout(bootstrap_replicates_layout)
        
        # Connect bootstrap checkbox to enable/disable replicates spinbox
        self.bootstrap_check.toggled.connect(self.bootstrap_replicates.setEnabled)
        
        bootstrap_group.setLayout(bootstrap_layout)
        basic_layout.addWidget(bootstrap_group)
        
        # Alignment cleaning tab
        alignment_tab = QWidget()
        alignment_layout = QVBoxLayout(alignment_tab)
        
        # Alignment cleaning options
        trim_group = QGroupBox("Trim Poorly Aligned Regions")
        trim_layout = QVBoxLayout()
        
        self.trim_check = QCheckBox("Trim alignment before tree building")
        trim_layout.addWidget(self.trim_check)
        
        # Trim method options
        trim_method_layout = QHBoxLayout()
        self.trimal_button = QPushButton("Trim with trimAl")
        self.gblocks_button = QPushButton("Trim with Gblocks")
        trim_method_layout.addWidget(self.trimal_button)
        trim_method_layout.addWidget(self.gblocks_button)
        trim_layout.addLayout(trim_method_layout)
        
        # Connect trim buttons to flag trimming
        self.trim_check.toggled.connect(self.toggle_trim_tools)
        self.trimal_button.clicked.connect(self.on_trimal_clicked)
        self.gblocks_button.clicked.connect(self.on_gblocks_clicked)
        
        # Add right-click context menus for advanced options
        self.gblocks_button.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.gblocks_button.customContextMenuRequested.connect(self.show_gblocks_context_menu)
        
        # Trim parameters
        trim_params_layout = QFormLayout()
        self.gap_threshold = QSpinBox()
        self.gap_threshold.setRange(0, 100)
        self.gap_threshold.setValue(50)
        self.gap_threshold.setSuffix("%")
        trim_params_layout.addRow("Gap threshold:", self.gap_threshold)
        
        self.min_conservation = QSpinBox()
        self.min_conservation.setRange(0, 100)
        self.min_conservation.setValue(60)
        self.min_conservation.setSuffix("%")
        trim_params_layout.addRow("Minimum conservation:", self.min_conservation)
        
        trim_layout.addLayout(trim_params_layout)
        trim_group.setLayout(trim_layout)
        alignment_layout.addWidget(trim_group)
        
        # Additional DNA barcode options
        barcode_group = QGroupBox("DNA Barcode Analysis Options")
        barcode_layout = QVBoxLayout()
        
        self.partition_check = QCheckBox("Use separate models for codon positions")
        barcode_layout.addWidget(self.partition_check)
        
        self.ascertainment_check = QCheckBox("Apply ascertainment bias correction")
        barcode_layout.addWidget(self.ascertainment_check)
        
        self.outgroup_layout = QHBoxLayout()
        self.outgroup_layout.addWidget(QLabel("Outgroup sequence:"))
        self.outgroup_edit = QLineEdit()
        self.outgroup_edit.setPlaceholderText("Enter sequence name for outgroup")
        self.outgroup_layout.addWidget(self.outgroup_edit)
        barcode_layout.addLayout(self.outgroup_layout)
        
        barcode_group.setLayout(barcode_layout)
        alignment_layout.addWidget(barcode_group)
        
        # Advanced tab
        advanced_tab = QWidget()
        advanced_layout = QVBoxLayout(advanced_tab)
        
        # RAxML binary selection
        raxml_group = QGroupBox("RAxML Binary")
        raxml_layout = QVBoxLayout()
        
        self.raxml_auto_radio = QRadioButton("Auto-detect RAxML binary")
        self.raxml_auto_radio.setChecked(True)
        raxml_layout.addWidget(self.raxml_auto_radio)
        
        self.raxml_ng_radio = QRadioButton("Use raxml-ng")
        raxml_layout.addWidget(self.raxml_ng_radio)
        
        self.raxml_pthreads_radio = QRadioButton("Use raxmlHPC-PTHREADS-AVX")
        raxml_layout.addWidget(self.raxml_pthreads_radio)
        
        self.raxml_std_radio = QRadioButton("Use standard raxml binary")
        raxml_layout.addWidget(self.raxml_std_radio)
        
        raxml_group.setLayout(raxml_layout)
        advanced_layout.addWidget(raxml_group)
        
        # Additional parameters
        params_group = QGroupBox("Additional Parameters")
        params_layout = QFormLayout()
        
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 64)
        self.threads_spinbox.setValue(8)  # Updated default to 8 threads
        params_layout.addRow("Number of threads:", self.threads_spinbox)
        
        self.seed_spinbox = QSpinBox()
        self.seed_spinbox.setRange(1, 99999)
        self.seed_spinbox.setValue(12345)
        params_layout.addRow("Random seed:", self.seed_spinbox)
        
        self.additional_args = QLineEdit()
        self.additional_args.setPlaceholderText("Additional arguments for RAxML")
        params_layout.addRow("Additional arguments:", self.additional_args)
        
        params_group.setLayout(params_layout)
        advanced_layout.addWidget(params_group)
        
        # Add tabs to tab widget
        tabs.addTab(basic_tab, "Basic Settings")
        tabs.addTab(alignment_tab, "Alignment & Barcode")
        tabs.addTab(advanced_tab, "Advanced")
        
        layout.addWidget(tabs)
        
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
    
    def toggle_trim_tools(self, checked):
        """Enable or disable trim tool buttons based on checkbox state"""
        logger = logging.getLogger("treecraft.raxml_dialog")
        
        # Log detailed information about the toggle operation
        logger.info(f"Toggle trim tools called with checked={checked}")
        logger.info(f"Current tool availability: trimal={self.has_trimal}, gblocks={self.has_gblocks}")
        
        if hasattr(self, 'gblocks_path'):
            logger.info(f"Gblocks path: {self.gblocks_path}")
        if hasattr(self, 'trimal_path'):
            logger.info(f"trimAl path: {self.trimal_path}")
        
        # First handle trimal button - enable if checked AND trimal is available
        trimal_enabled = checked and self.has_trimal
        self.trimal_button.setEnabled(trimal_enabled)
        
        # Handle gblocks button - enable if checked AND gblocks is available
        gblocks_enabled = checked and self.has_gblocks
        self.gblocks_button.setEnabled(gblocks_enabled)
        
        # Update the button appearance to provide visual feedback
        self.trimal_button.setStyleSheet(
            "font-weight: bold;" if self.trim_tool == "trimal" else "font-weight: normal;"
        )
        self.gblocks_button.setStyleSheet(
            "font-weight: bold;" if self.trim_tool == "gblocks" else "font-weight: normal;"
        )
        
        # Log detailed information about the result
        logger.info(f"Trim tools toggle result:")
        logger.info(f"  Trim checkbox checked: {checked}")
        logger.info(f"  trimAl button enabled: {trimal_enabled}")
        logger.info(f"  Gblocks button enabled: {gblocks_enabled}")
        logger.info(f"  Current trim tool: {self.trim_tool}")
        
        # Parameters should be enabled if checkbox is checked
        self.gap_threshold.setEnabled(checked)
        self.min_conservation.setEnabled(checked)
        
    def on_trimal_clicked(self):
        """Handle trimAl button click"""
        # Log that the trimAl button was clicked
        logger = logging.getLogger("treecraft.raxml_dialog")
        logger.info(f"trimAl button clicked - has_trimal: {self.has_trimal}")
        
        # Only proceed if trimAl is available
        if not self.has_trimal:
            QMessageBox.warning(self, "trimAl Not Found", 
                              "trimAl was not found in your PATH. Please install trimAl and make sure it's in your PATH.")
            return
            
        # Store in a class variable
        self.trim_tool = "trimal"
        
        # Make sure the trim checkbox is checked
        self.trim_check.setChecked(True)
        
        # Update button appearance to indicate which is selected
        self.trimal_button.setStyleSheet("font-weight: bold;")
        self.gblocks_button.setStyleSheet("font-weight: normal;")
        
        # Show a message about what will happen
        QMessageBox.information(self, "Trim Alignment", 
                              "trimAl will be used to trim poorly aligned regions from the alignment.")
    
    def on_gblocks_clicked(self):
        """Handle Gblocks button click"""
        # Get logger
        logger = logging.getLogger("treecraft.raxml_dialog")
        
        # Log detailed information about the Gblocks button click
        logger.info("Gblocks button clicked - starting detailed diagnostics:")
        logger.info(f"  has_gblocks flag: {self.has_gblocks}")
        if hasattr(self, 'gblocks_path'):
            logger.info(f"  gblocks_path: {self.gblocks_path}")
        
        # Double-check Gblocks availability to ensure it's up-to-date
        from treecraft.utils.external_tools import find_external_tool
        fresh_gblocks_path = find_external_tool("gblocks")
        logger.info(f"  Re-checking Gblocks availability: {fresh_gblocks_path}")
        
        # Update our flag based on the fresh check
        self.has_gblocks = fresh_gblocks_path is not None
        self.gblocks_path = fresh_gblocks_path
        logger.info(f"  Updated has_gblocks flag: {self.has_gblocks}")
        
        # Only proceed if Gblocks is available
        if not self.has_gblocks:
            # Offer the user a chance to locate Gblocks manually
            from PyQt6.QtWidgets import QFileDialog
            
            error_message = "Gblocks was not found in your PATH. Please install Gblocks and make sure it's in your PATH."
            
            # Ask if they want to locate it manually
            reply = QMessageBox.question(
                self, 
                "Gblocks Not Found", 
                error_message + "\n\nWould you like to locate the Gblocks executable manually?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # Open file dialog to select Gblocks executable
                file_path, _ = QFileDialog.getOpenFileName(
                    self, 
                    "Locate Gblocks Executable",
                    "",
                    "All Files (*)"
                )
                
                if file_path:
                    logger.info(f"  User selected custom Gblocks path: {file_path}")
                    
                    # Verify if the selected file is executable
                    import os
                    if os.path.isfile(file_path) and os.access(file_path, os.X_OK):
                        # Success! Use this path
                        self.gblocks_path = file_path
                        self.has_gblocks = True
                        logger.info(f"  Custom Gblocks path verified: {file_path}")
                        
                        # Update tooltip
                        self.gblocks_button.setToolTip(f"Gblocks: Available (Custom Path: {file_path})")
                    else:
                        # Not a valid executable
                        error_message = f"The selected file '{file_path}' is not an executable.\nPlease select a valid Gblocks executable."
                        logger.error(f"  Error: {error_message}")
                        QMessageBox.warning(self, "Invalid Executable", error_message)
                        return
                else:
                    # User canceled the file dialog
                    logger.info("  User canceled manual Gblocks selection")
                    return
            else:
                # User declined to locate manually
                logger.error(f"  Error: {error_message}")
                return
            
        # If we're here, Gblocks is available
        logger.info("  Gblocks is available, proceeding with selection")
        
        # Store in a class variable
        self.trim_tool = "gblocks"
        logger.info(f"  Setting trim_tool to: {self.trim_tool}")
        
        # Make sure the trim checkbox is checked
        self.trim_check.setChecked(True)
        
        # Update button appearance to indicate which is selected
        self.trimal_button.setStyleSheet("font-weight: normal;")
        self.gblocks_button.setStyleSheet("font-weight: bold;")
        logger.info("  Updated button styles to highlight Gblocks selection")
        
        # Show a message about what will happen
        success_message = "Gblocks will be used to trim poorly aligned regions from the alignment."
        logger.info(f"  Success: {success_message}")
        QMessageBox.information(self, "Trim Alignment", success_message)
        
        # Make sure the button is actually enabled
        self.gblocks_button.setEnabled(True)
        logger.info("  Ensured Gblocks button is enabled")
    
    def on_modeltest_clicked(self):
        """Handle ModelTest-NG button click to find best model immediately"""
        # Store in class variables
        self.find_model = True
        self.model_tool = "modeltest-ng"
        
        # Update button appearance
        self.findmodel_button.setStyleSheet("font-weight: bold;")
        self.iqtree_button.setStyleSheet("font-weight: normal;")
        
        # Get logger for detailed information
        logger = logging.getLogger("treecraft.raxml_dialog")
        
        # Find ModelTest-NG path
        modeltest_path = find_external_tool("modeltest-ng")
        if not modeltest_path:
            QMessageBox.warning(self, "ModelTest-NG Not Found", 
                              "ModelTest-NG was not found in your PATH. Please install ModelTest-NG and make sure it's in your PATH.")
            return
            
        # Check for alignment
        alignment = self.parent().current_alignment
        if not alignment or len(alignment) < 3:
            QMessageBox.warning(self, "No Alignment", 
                               "No alignment available or too few sequences. Please align at least 3 sequences first.")
            return
        
        # Create temporary file for alignment
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
            input_file = tmp.name
            
        # Write alignment to the temporary file
        logger.info(f"Writing alignment to temporary file: {input_file}")
        try:
            # Create the file with alignment data
            with open(input_file, 'w') as f:
                for record in alignment:
                    f.write(f">{record.id}\n{record.seq}\n")
                    
            # Create a custom progress dialog with a label for the command
            from PyQt6.QtWidgets import QVBoxLayout, QDialog, QProgressBar, QPushButton, QLabel
            
            progress_dialog = QDialog(self)
            progress_dialog.setWindowTitle("Finding Best Model")
            progress_dialog.setModal(True)
            progress_dialog.resize(700, 200)
            
            # Create layout
            progress_layout = QVBoxLayout(progress_dialog)
            
            # Add status label
            status_label = QLabel("Running ModelTest-NG to find best model...")
            progress_layout.addWidget(status_label)
            
            # Add command label
            actual_cmd = f"{modeltest_path} -i {input_file} -p 4 -d nt"
            cmd_label = QLabel(f"Command: {actual_cmd}")
            cmd_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            progress_layout.addWidget(cmd_label)
            
            # Add progress bar
            progress_bar = QProgressBar()
            progress_bar.setRange(0, 100)
            progress_bar.setValue(0)
            progress_layout.addWidget(progress_bar)
            
            # Add Cancel button
            button_layout = QHBoxLayout()
            cancel_button = QPushButton("Cancel")
            cancel_button.clicked.connect(progress_dialog.reject)
            button_layout.addWidget(cancel_button)
            progress_layout.addLayout(button_layout)
            
            # Create custom progress object that can be passed to external_tools.py
            class Progress:
                def __init__(self, dialog, progress_bar, status_label):
                    self.dialog = dialog
                    self.progress_bar = progress_bar
                    self.status_label = status_label
                    self._cancelled = False
                    
                def setValue(self, value):
                    self.progress_bar.setValue(value)
                    QApplication.processEvents()
                    
                def wasCanceled(self):
                    return self._cancelled
                    
                def setLabelText(self, text):
                    self.status_label.setText(text)
                    QApplication.processEvents()
                    
            # Create progress object
            progress = Progress(progress_dialog, progress_bar, status_label)
            
            # Connect cancel button
            cancel_button.clicked.connect(lambda: setattr(progress, '_cancelled', True))
            
            # Show the dialog
            progress_dialog.show()
            QApplication.processEvents()
            
            # Run ModelTest-NG
            from treecraft.utils.external_tools import run_modeltest
            model_results = run_modeltest(input_file, None, progress)
            
            # Close the progress dialog
            progress_dialog.accept()
            
            # Process the results
            if model_results and isinstance(model_results, dict) and model_results.get("best"):
                # Get best model (which should be from AICc)
                best_model = model_results["best"]
                logger.info(f"Found best model with ModelTest-NG: {best_model}")
                
                # Add all models found to the dropdown if not already there
                for criterion in ["AICc", "AIC", "BIC"]:
                    model_value = model_results.get(criterion)
                    if model_value:
                        # Check if already in dropdown
                        model_in_dropdown = False
                        for i in range(self.model_combobox.count()):
                            if self.model_combobox.itemText(i) == model_value:
                                model_in_dropdown = True
                                if model_value == best_model:
                                    # If this is our best model (AICc), select it
                                    self.model_combobox.setCurrentIndex(i)
                                break
                        
                        # Add if not found
                        if not model_in_dropdown:
                            self.model_combobox.addItem(model_value)
                            logger.info(f"Added {criterion} model {model_value} to dropdown")
                            
                # Make sure the AICc model (best_model) is selected
                self.model_combobox.setCurrentText(best_model)
                
                # Create a results dialog with the information
                result_dialog = QDialog(self)
                result_dialog.setWindowTitle("ModelTest-NG Results")
                result_dialog.resize(600, 400)
                
                result_layout = QVBoxLayout(result_dialog)
                
                # Add a header
                header_label = QLabel("<h3>ModelTest-NG Results</h3>")
                header_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                result_layout.addWidget(header_label)
                
                # Create a table showing all three criteria results
                result_table = """
                <style>
                table { border-collapse: collapse; width: 100%; }
                th, td { padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }
                th { background-color: #f2f2f2; }
                .selected { background-color: #e6ffe6; font-weight: bold; }
                </style>
                <table>
                  <tr>
                    <th>Criterion</th>
                    <th>Best Model</th>
                    <th>Description</th>
                  </tr>
                  <tr class='%s'>
                    <td>AICc</td>
                    <td>%s</td>
                    <td>Corrected Akaike Information Criterion</td>
                  </tr>
                  <tr class='%s'>
                    <td>AIC</td>
                    <td>%s</td>
                    <td>Akaike Information Criterion</td>
                  </tr>
                  <tr class='%s'>
                    <td>BIC</td>
                    <td>%s</td>
                    <td>Bayesian Information Criterion</td>
                  </tr>
                </table>
                <p><b>Selected model: %s</b> (based on AICc)</p>
                """ % (
                    "selected" if best_model == model_results.get("AICc") else "",
                    model_results.get("AICc", "Not found"),
                    "selected" if best_model == model_results.get("AIC") else "",
                    model_results.get("AIC", "Not found"),
                    "selected" if best_model == model_results.get("BIC") else "",
                    model_results.get("BIC", "Not found"),
                    best_model
                )
                
                results_view = QLabel(result_table)
                result_layout.addWidget(results_view)
                
                # Add explanation
                explanation = QLabel(
                    f"<p>ModelTest-NG analyzed your alignment using three different criteria. "
                    f"It's generally recommended to use the AICc criterion for most datasets, "
                    f"which selected <b>{best_model}</b>.</p>"
                    f"<p>The AICc model (<b>{best_model}</b>) has been automatically selected in the dropdown list.</p>"
                    f"<p>All models found have been added to the dropdown:</p>"
                    f"<ul>"
                    f"<li>AICc: <b>{model_results.get('AICc', 'Not found')}</b></li>"
                    f"<li>AIC: {model_results.get('AIC', 'Not found')}</li>"
                    f"<li>BIC: {model_results.get('BIC', 'Not found')}</li>"
                    f"</ul>"
                    f"<p><b>Command used:</b><br>{actual_cmd}</p>"
                )
                explanation.setWordWrap(True)
                result_layout.addWidget(explanation)
                
                # Add close button
                close_button = QPushButton("Close")
                close_button.clicked.connect(result_dialog.accept)
                result_layout.addWidget(close_button)
                
                # Show the dialog
                result_dialog.exec()
            else:
                # If model finding failed, try to parse the file directly as a last resort
                output_file = f"{input_file}.models.out"
                if os.path.exists(output_file):
                    try:
                        # Try to read and parse the output file directly
                        logger.info(f"Attempting direct parsing of output file: {output_file}")
                        models = {"BIC": None, "AIC": None, "AICc": None}
                        
                        # Read the content
                        with open(output_file, 'r') as f:
                            content = f.read()
                        
                        # Look for each criterion using grep-like approach
                        for criterion in ["BIC", "AIC", "AICc"]:
                            cmd = ["grep", "-A", "3", f"Best model according to {criterion}", output_file]
                            result = subprocess.run(cmd, capture_output=True, text=True)
                            if result.returncode == 0 and result.stdout:
                                lines = result.stdout.strip().split('\n')
                                for i, line in enumerate(lines):
                                    if line.startswith("Model:"):
                                        model_value = line.split(":", 1)[1].strip()
                                        models[criterion] = model_value
                                        logger.info(f"Direct grep found {criterion} model: {model_value}")
                                        break
                        
                        # If we found any models, show them to the user
                        if any(models.values()):
                            best_model = models["AICc"] or models["AIC"] or models["BIC"]
                            
                            # Show success dialog with model info
                            QMessageBox.information(self, "Model Found", 
                                               f"ModelTest-NG found the best model: {best_model} (AICc criterion)\n\n"
                                               f"BIC: {models['BIC'] or 'Not found'}\n"
                                               f"AIC: {models['AIC'] or 'Not found'}\n"
                                               f"AICc: {models['AICc'] or 'Not found'}\n\n"
                                               f"The AICc model ({best_model}) has been automatically selected in the dropdown list.\n"
                                               "The AICc criterion is recommended for most phylogenetic analyses.")
                            
                            # Add the model to the dropdown and select it
                            self.model_combobox.addItem(best_model)
                            self.model_combobox.setCurrentText(best_model)
                            
                            # Return early - we've successfully recovered
                            return
                    except Exception as e:
                        logger.error(f"Error in direct output parsing: {e}")
                
                # If all else fails, show an error with the exact command that was run
                QMessageBox.warning(self, "Model Finding Failed", 
                                  "ModelTest-NG was unable to find the best model. "
                                  "Check that your alignment is valid and try again.\n\n"
                                  f"Command that was run:\n{actual_cmd}\n\n"
                                  "The ModelTest-NG output file may exist at:\n"
                                  f"{input_file}.models.out")
        except Exception as e:
            import traceback
            logger.error(f"Error running ModelTest-NG: {e}")
            logger.error(traceback.format_exc())
            
            # Print more details for debugging
            logger.error(f"Variables at error time: input_file={input_file if 'input_file' in locals() else 'NOT DEFINED'}")
            logger.error(f"Error location: {traceback.format_exc()}")
            
            QMessageBox.critical(self, "Error", f"An error occurred while running ModelTest-NG: {str(e)}")
        finally:
            # Clean up temporary file
            try:
                if os.path.exists(input_file):
                    os.unlink(input_file)
            except Exception as cleanup_error:
                logger.error(f"Error cleaning up temp file: {cleanup_error}")
                pass
    
    def on_iqtree_clicked(self):
        """Handle IQ-TREE button click to find best model immediately"""
        # Store in class variables
        self.find_model = True
        self.model_tool = "iqtree"
        
        # Update button appearance
        self.findmodel_button.setStyleSheet("font-weight: normal;")
        self.iqtree_button.setStyleSheet("font-weight: bold;")
        
        # Get logger for detailed information
        logger = logging.getLogger("treecraft.raxml_dialog")
        
        # Find IQ-TREE path
        iqtree_path = find_external_tool("iqtree")
        if not iqtree_path:
            QMessageBox.warning(self, "IQ-TREE Not Found", 
                              "IQ-TREE was not found in your PATH. Please install IQ-TREE and make sure it's in your PATH.")
            return
            
        # Create a custom progress dialog with a label for the command
        from PyQt6.QtWidgets import QVBoxLayout, QDialog, QProgressBar, QPushButton, QLabel
        
        progress_dialog = QDialog(self)
        progress_dialog.setWindowTitle("Finding Best Model")
        progress_dialog.setModal(True)
        progress_dialog.resize(700, 200)
        
        # Create layout
        progress_layout = QVBoxLayout(progress_dialog)
        
        # Add status label
        status_label = QLabel("Running IQ-TREE to find best model...")
        progress_layout.addWidget(status_label)
        
        # Add command label
        cmd_label = QLabel(f"Command: {iqtree_path} -s [alignment] -m MF")
        cmd_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        cmd_label.setWordWrap(True)
        progress_layout.addWidget(cmd_label)
        
        # Add progress bar
        progress_bar = QProgressBar()
        progress_bar.setRange(0, 100)
        progress_bar.setValue(10)
        progress_layout.addWidget(progress_bar)
        
        # Add cancel button
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(progress_dialog.reject)
        progress_layout.addWidget(cancel_button)
        
        # Show the dialog
        progress_dialog.show()
        QApplication.processEvents()
        
        # Create a proxy progress object that updates our custom dialog
        class ProgressProxy:
            def __init__(self, dialog, status_label, cmd_label, progress_bar):
                self.dialog = dialog
                self.status_label = status_label
                self.cmd_label = cmd_label
                self.progress_bar = progress_bar
                self.wasCanceled = lambda: dialog.result() == QDialog.DialogCode.Rejected
                
            def setValue(self, value):
                self.progress_bar.setValue(value)
                QApplication.processEvents()
                
            def setLabelText(self, text):
                # Split text to separate status from command if there's a command
                parts = text.split("\n\nCommand:", 1)
                if len(parts) > 1:
                    self.status_label.setText(parts[0])
                    self.cmd_label.setText("Command: " + parts[1])
                else:
                    self.status_label.setText(text)
                QApplication.processEvents()
                
            def close(self):
                self.dialog.accept()
        
        # Create our proxy progress object
        progress = ProgressProxy(progress_dialog, status_label, cmd_label, progress_bar)
        
        # Create temporary files for the alignment
        import tempfile
        temp_dir = tempfile.mkdtemp(prefix="iqtree_")
        input_file = os.path.join(temp_dir, "alignment.fasta")
        
        # We need to get the current alignment from the main window
        main_window = self.parent()
        alignment = None
        if hasattr(main_window, 'current_alignment') and main_window.current_alignment:
            alignment = main_window.current_alignment
            
            # Write the alignment to the temporary file
            with open(input_file, 'w') as f:
                for record in alignment:
                    f.write(f">{record.id}\n{record.seq}\n")
                    
            # Run IQ-TREE model finder
            from treecraft.utils.external_tools import run_iqtree_modeltest
            
            # Update progress with actual command
            actual_cmd = f"{iqtree_path} -s {input_file} -m MF -nt auto -pre {input_file}"
            progress.setValue(20)
            progress.status_label.setText("Running IQ-TREE Model Finder (this may take 1-2 minutes)...")
            progress.cmd_label.setText(f"Command: {actual_cmd}")
            progress.cmd_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            QApplication.processEvents()
            
            # Create a timer to update progress periodically to show activity
            from PyQt6.QtCore import QTimer
            timer = QTimer()
            current_progress = 20
            
            def update_progress_animation():
                nonlocal current_progress
                current_progress = (current_progress + 1) % 80
                if current_progress < 20:
                    current_progress = 20
                progress.setValue(current_progress)
                progress.status_label.setText(f"Running IQ-TREE Model Finder (this may take 1-2 minutes)...{current_progress}%")
                QApplication.processEvents()
            
            # Start the timer to update progress every 0.5 seconds
            timer.timeout.connect(update_progress_animation)
            timer.start(500)  # Update every 500 ms
            
            # Run IQ-TREE
            best_model = run_iqtree_modeltest(input_file, None, progress)
            
            # Stop the timer once we have the result
            timer.stop()
            
            # Now that we have the best model, update the combobox
            if best_model:
                logger.info(f"Found best model with IQ-TREE: {best_model}")
                
                # Check if the model is in our list, if not, add it
                model_found = False
                for i in range(self.model_combobox.count()):
                    if self.model_combobox.itemText(i) == best_model:
                        self.model_combobox.setCurrentIndex(i)
                        model_found = True
                        break
                
                if not model_found:
                    self.model_combobox.addItem(best_model)
                    self.model_combobox.setCurrentText(best_model)
                
                # Look for specific log files with model information
                model_info = ""
                iqtree_log_file = f"{input_file}.iqtree"
                
                if os.path.exists(iqtree_log_file):
                    try:
                        with open(iqtree_log_file, 'r') as f:
                            log_content = f.read()
                            
                            # Extract the most relevant information
                            model_sections = []
                            
                            # Find ModelFinder section
                            modelfinder_pos = log_content.find("ModelFinder")
                            if modelfinder_pos > -1:
                                model_sections.append(log_content[modelfinder_pos:modelfinder_pos+500])
                                
                            # Find best model section
                            best_model_pos = log_content.find("Best-fit model")
                            if best_model_pos > -1:
                                model_sections.append(log_content[best_model_pos:best_model_pos+500])
                                
                            # Format the info for display
                            model_info = "\n".join([section.split("\n")[0] for section in model_sections if section])
                    except Exception as e:
                        logger.error(f"Error reading IQ-TREE log: {e}")
                
                # Create a dialog with more detailed information
                detail_dialog = QDialog(self)
                detail_dialog.setWindowTitle("IQ-TREE Results")
                detail_dialog.resize(600, 400)
                
                detail_layout = QVBoxLayout(detail_dialog)
                
                # Add a label with the result
                result_label = QLabel(f"<b>Best evolutionary model found: {best_model}</b>")
                result_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                detail_layout.addWidget(result_label)
                
                # Add a button to show the full log
                log_button = QPushButton("Show IQ-TREE Log")
                detail_layout.addWidget(log_button)
                
                # Add the command line used
                cmd_label = QLabel("<b>Command used:</b>")
                detail_layout.addWidget(cmd_label)
                
                cmd_text = QTextEdit()
                cmd_text.setReadOnly(True)
                cmd_text.setPlainText(actual_cmd)
                cmd_text.setMaximumHeight(80)
                detail_layout.addWidget(cmd_text)
                
                # Add the output
                output_label = QLabel("<b>IQ-TREE Output:</b>")
                detail_layout.addWidget(output_label)
                
                # Get the log file content from multiple possible locations
                log_files = [
                    f"{input_file}.iqtree",  # Primary output file
                    f"{input_file}.log",     # Log file
                    f"{input_file}.iqtree.log" # Alternative log location
                ]
                
                log_content = ""
                for log_file_path in log_files:
                    if os.path.exists(log_file_path):
                        try:
                            with open(log_file_path, 'r') as f:
                                log_content = f.read()
                                logger.info(f"Found IQ-TREE log file: {log_file_path}")
                                break
                        except Exception as e:
                            logger.error(f"Error reading log file {log_file_path}: {e}")
                
                # Add the summary output to the dialog
                output_text = QTextEdit()
                output_text.setReadOnly(True)
                
                # Extract the most important parts of the log
                important_lines = []
                if log_content:
                    lines = log_content.splitlines()
                    for i, line in enumerate(lines):
                        # Look for key information in the log
                        if any(keyword in line for keyword in ["BIC:", "AIC:", "Best-fit model", "ModelFinder", "SUBSTITUTION", "TPM", "GTR", "+F+", "+G", "+R"]):
                            # Include some context (5 lines before and after)
                            start = max(0, i - 5)
                            end = min(len(lines), i + 5)
                            important_lines.extend(lines[start:end])
                            important_lines.append("")  # Add a blank line between sections
                
                output_text.setPlainText("\n".join(important_lines))
                detail_layout.addWidget(output_text)
                
                # Function to show the full log
                def show_full_log():
                    log_dialog = QDialog(detail_dialog)
                    log_dialog.setWindowTitle("IQ-TREE Full Log")
                    log_dialog.resize(800, 600)
                    
                    log_layout = QVBoxLayout(log_dialog)
                    
                    # Add a note about the model detection
                    note_label = QLabel(f"<b>Selected Model:</b> {best_model}")
                    log_layout.addWidget(note_label)
                    
                    # Add information about which log file we're displaying
                    source_label = QLabel("<b>Log file source:</b>")
                    for log_file_path in log_files:
                        if os.path.exists(log_file_path):
                            source_label.setText(source_label.text() + f" {os.path.basename(log_file_path)}")
                            break
                    log_layout.addWidget(source_label)
                    
                    # Add the full log content with monospace font
                    full_log_text = QTextEdit()
                    full_log_text.setReadOnly(True)
                    full_log_text.setPlainText(log_content)
                    full_log_text.setStyleSheet("font-family: monospace;")
                    
                    # Try to find and highlight the best model line
                    if "Best-fit model" in log_content:
                        cursor = full_log_text.textCursor()
                        # Save current format
                        current_format = QTextCharFormat()
                        current_format.setFontFamily("monospace")
                        
                        # Create a highlighted format
                        highlight_format = QTextCharFormat()
                        highlight_format.setBackground(QColor(255, 255, 0, 100))  # Light yellow background
                        highlight_format.setFontWeight(QFont.Weight.Bold)
                        highlight_format.setFontFamily("monospace")
                        
                        # Find the line with the best model
                        full_log_text.setPlainText("")  # Clear text
                        for line in log_content.splitlines():
                            if "Best-fit model" in line:
                                # Add this line with highlighting
                                cursor.insertText(line + "\n", highlight_format)
                            else:
                                # Add normal line
                                cursor.insertText(line + "\n", current_format)
                    else:
                        # Just use the plain log if we can't find the model line
                        full_log_text.setPlainText(log_content)
                                
                    log_layout.addWidget(full_log_text)
                    
                    close_button = QPushButton("Close")
                    close_button.clicked.connect(log_dialog.accept)
                    log_layout.addWidget(close_button)
                    
                    log_dialog.exec()
                
                # Connect the log button
                log_button.clicked.connect(show_full_log)
                
                # Add a close button
                close_button = QPushButton("Close")
                close_button.clicked.connect(detail_dialog.accept)
                detail_layout.addWidget(close_button)
                
                # Show the dialog
                detail_dialog.exec()
            else:
                # Check if the iqtree logfile exists regardless of whether a model was found
                iqtree_log_file = f"{input_file}.iqtree"
                best_model_file = f"{input_file}.bestmodel.txt"
                logfile_exists = os.path.exists(iqtree_log_file)
                model_file_exists = os.path.exists(best_model_file)
                
                if logfile_exists:
                    # The process may have run correctly but we failed to parse the model
                    # Try to extract the model directly from the log file for display
                    model_from_log = None
                    try:
                        with open(iqtree_log_file, 'r') as f:
                            for line in f:
                                if "Best-fit model according to BIC:" in line:
                                    model_from_log = line.split("Best-fit model according to BIC:", 1)[1].strip()
                                    break
                    except Exception as e:
                        logger.error(f"Error extracting model from log file: {e}")
                    
                    if model_from_log:
                        # Found the model - show it and add it to the dropdown
                        QMessageBox.information(self, "Model Found Manually", 
                                           f"IQ-TREE found the best model: {model_from_log}\n\n"
                                           "This model has been added to the dropdown list. You may now\n"
                                           "continue building your tree with this model.")
                        
                        # Add the model to the dropdown and select it
                        self.model_combobox.addItem(model_from_log)
                        self.model_combobox.setCurrentText(model_from_log)
                    elif model_file_exists:
                        # We have a best model file with additional info
                        model_info = ""
                        try:
                            with open(best_model_file, 'r') as f:
                                model_info = f.read()
                        except:
                            model_info = "Check the model file for details."
                        
                        QMessageBox.warning(self, "Model Parsing Failed", 
                                          "IQ-TREE ran successfully but we couldn't automatically identify the best model.\n\n"
                                          f"Log file was created at: {iqtree_log_file}\n\n"
                                          "To find the model manually, open the log file and search for this text:\n"
                                          "  'Best-fit model according to BIC:'\n\n"
                                          f"Additional info:\n{model_info}")
                    else:
                        # Basic error with just the log file path
                        QMessageBox.warning(self, "Model Parsing Failed", 
                                          "IQ-TREE ran successfully but we couldn't automatically identify the best model.\n\n"
                                          f"Log file was created at: {iqtree_log_file}\n\n"
                                          "Please open this file and look for a line that says:\n"
                                          "  'Best-fit model according to BIC:'\n\n"
                                          "Once you find the model, you can manually enter it in the dropdown.")
                else:
                    # IQ-TREE likely encountered an error
                    QMessageBox.warning(self, "Model Finding Failed", 
                                      "IQ-TREE did not complete the model finding process correctly.\n\n"
                                      "The process may have been interrupted or timed out. The model finding can take 1-2 minutes to complete.\n\n"
                                      "You can try running the iqtree command manually in a terminal window.")
        else:
            progress.close()
            QMessageBox.warning(self, "No Alignment", 
                               "No alignment available. Please align sequences first.")
                               
    def show_gblocks_context_menu(self, position):
        """Show context menu for Gblocks button with advanced options"""
        from PyQt6.QtWidgets import QMenu, QFileDialog
        import os
        
        # Create logger for this method
        logger = logging.getLogger("treecraft.raxml_dialog")
        logger.info("Opening Gblocks context menu")
        
        # Create context menu
        context_menu = QMenu()
        
        # Add option to set custom Gblocks path
        set_path_action = context_menu.addAction("Set Custom Gblocks Path...")
        
        # Add option to check current Gblocks path
        check_path_action = context_menu.addAction("Check Current Gblocks Path")
        
        # Add refresh option
        refresh_action = context_menu.addAction("Refresh Gblocks Detection")
        
        # Show the menu and get the selected action
        action = context_menu.exec(self.gblocks_button.mapToGlobal(position))
        
        # Handle the selected action
        if action == set_path_action:
            # Open file dialog to select Gblocks executable
            file_path, _ = QFileDialog.getOpenFileName(
                self, 
                "Select Gblocks Executable",
                "",
                "All Files (*)"
            )
            
            if file_path:
                logger.info(f"User selected custom Gblocks path: {file_path}")
                
                # Verify if the selected file is executable
                if os.path.isfile(file_path) and os.access(file_path, os.X_OK):
                    # Success! Use this path
                    self.gblocks_path = file_path
                    self.has_gblocks = True
                    
                    # Update tooltip with the custom path
                    self.gblocks_button.setToolTip(f"Gblocks: Available (Custom Path: {file_path})")
                    
                    # Enable the button if trim is checked
                    if self.trim_check.isChecked():
                        self.gblocks_button.setEnabled(True)
                    
                    # Show confirmation to the user
                    QMessageBox.information(self, "Custom Gblocks Path", 
                                          f"Custom Gblocks path set to:\n{file_path}")
                    
                    logger.info(f"Custom Gblocks path set successfully to {file_path}")
                else:
                    # Not a valid executable
                    QMessageBox.warning(self, "Invalid Executable", 
                                      f"The selected file '{file_path}' is not an executable.\nPlease select a valid Gblocks executable.")
                    logger.warning(f"Selected file is not a valid executable: {file_path}")
        
        elif action == check_path_action:
            # Show the current Gblocks path
            if hasattr(self, 'gblocks_path') and self.gblocks_path:
                # Check if the path is valid
                if os.path.isfile(self.gblocks_path) and os.access(self.gblocks_path, os.X_OK):
                    status = "Valid executable"
                else:
                    status = "INVALID - Not an executable file"
                
                QMessageBox.information(self, "Current Gblocks Path", 
                                      f"Current Gblocks path: {self.gblocks_path}\nStatus: {status}")
                logger.info(f"Current Gblocks path: {self.gblocks_path}, Status: {status}")
            else:
                QMessageBox.warning(self, "No Gblocks Path", 
                                  "No Gblocks path is currently set.")
                logger.warning("No Gblocks path is currently set")
        
        elif action == refresh_action:
            # Re-detect Gblocks
            from treecraft.utils.external_tools import find_external_tool
            logger.info("Refreshing Gblocks detection")
            
            fresh_gblocks_path = find_external_tool("gblocks")
            logger.info(f"Refreshed Gblocks path: {fresh_gblocks_path}")
            
            # Update our fields
            old_path = getattr(self, 'gblocks_path', None)
            old_has_gblocks = getattr(self, 'has_gblocks', False)
            
            self.gblocks_path = fresh_gblocks_path
            self.has_gblocks = fresh_gblocks_path is not None
            
            # Update the button state if the trim checkbox is checked
            if self.trim_check.isChecked():
                self.gblocks_button.setEnabled(self.has_gblocks)
            
            # Update tooltip
            if self.has_gblocks:
                self.gblocks_button.setToolTip(f"Gblocks: Available ({self.gblocks_path})")
            else:
                self.gblocks_button.setToolTip("Gblocks: Not found in PATH")
            
            # Show results to the user
            if self.has_gblocks:
                if old_has_gblocks and old_path != fresh_gblocks_path:
                    QMessageBox.information(self, "Gblocks Detection", 
                                          f"Gblocks detected at a different path:\n{fresh_gblocks_path}")
                elif not old_has_gblocks:
                    QMessageBox.information(self, "Gblocks Detection", 
                                          f"Gblocks found at:\n{fresh_gblocks_path}")
                else:
                    QMessageBox.information(self, "Gblocks Detection", 
                                          f"Gblocks confirmed at:\n{fresh_gblocks_path}")
            else:
                QMessageBox.warning(self, "Gblocks Not Found", 
                                  "Gblocks was not found in your PATH.")
                                  
            logger.info(f"Gblocks detection refresh complete. Has Gblocks: {self.has_gblocks}, Path: {self.gblocks_path}")
        
    def get_parameters(self):
        """Get the selected RAxML parameters"""
        # Determine which RAxML binary to use
        raxml_binary = "auto"
        if self.raxml_ng_radio.isChecked():
            raxml_binary = "raxml-ng"
        elif self.raxml_pthreads_radio.isChecked():
            raxml_binary = "raxmlHPC-PTHREADS-AVX"
        elif self.raxml_std_radio.isChecked():
            raxml_binary = "raxml"
        
        # Get other parameters
        params = {
            "model": self.model_combobox.currentText(),
            "bootstrap": self.bootstrap_check.isChecked(),
            "bootstrap_replicates": self.bootstrap_replicates.value() if self.bootstrap_check.isChecked() else 0,
            "trim_alignment": self.trim_check.isChecked(),
            "gap_threshold": self.gap_threshold.value(),
            "min_conservation": self.min_conservation.value(),
            "partition_codons": self.partition_check.isChecked(),
            "ascertainment_bias": self.ascertainment_check.isChecked(),
            "outgroup": self.outgroup_edit.text(),
            "raxml_binary": raxml_binary,
            "threads": self.threads_spinbox.value(),
            "seed": self.seed_spinbox.value(),
            "additional_args": self.additional_args.text()
        }
        
        # Add trim tool information if set
        if hasattr(self, 'trim_tool'):
            params["trim_tool"] = self.trim_tool
        
        # Add model finding flags if set
        if hasattr(self, 'find_model') and self.find_model:
            params["find_model"] = True
            params["model_tool"] = getattr(self, 'model_tool', "modeltest-ng")
        
        return params