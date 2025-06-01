"""
This file has been programmatically modified to improve handling of sequence names with spaces
in RaxML trees and tree node rename/delete operations.
"""

import os
import logging
import tempfile
import re
from datetime import datetime
from io import StringIO

from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QSplitter, QToolBar, 
    QLabel, QListWidget, QStatusBar, QMenu, QMenuBar, QFileDialog, QApplication,
    QMessageBox, QInputDialog, QListWidgetItem, QDialog, QComboBox, QPushButton,
    QFormLayout, QDialogButtonBox, QToolButton, QCheckBox, QSpinBox,
    QTabWidget, QScrollArea, QGridLayout, QTextEdit, QFrame, QColorDialog,
    QSystemTrayIcon
)
from PyQt6.QtCore import Qt, QSize, QThread, pyqtSignal, QObject, QTimer # QApplication, Qt are used
from PyQt6.QtGui import QAction, QActionGroup, QIcon, QFont, QColor, QPixmap

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Align import MultipleSeqAlignment

from treecraft.gui.sequence_list import SequenceListWidget
from treecraft.gui.phylo_canvas import PhyloCanvas
from treecraft.gui.sequence_dialog import AddSequenceDialog as SequenceDialog
from treecraft.gui.edit_sequence_dialog import EditSequenceDialog
from treecraft.gui.alignment_dialog import SequenceAlignmentDialog
from treecraft.gui.tree_dialog import TreeBuildingDialog as TreeBuildDialog
from treecraft.core.tree_builder import build_tree_with_params
from treecraft.gui.export_dialog import ExportDialog, ImageExportOptionsDialog
from treecraft.gui.raxml_dialog import RaxMLDialog
from treecraft.gui.debug_console import DebugConsole

# Add imports for external tools
from treecraft.utils.external_tools import (
    run_trimal, run_gblocks,
    run_modeltest, run_iqtree_modeltest,
    run_raxml, find_tree_file
)

class SequenceNameViewer(QTextEdit):
    """Widget for displaying sequence names"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setFontFamily("Courier New")
        self.setFontPointSize(9)
        self.dark_mode = False
        self.horizontal_scrollbar = self.horizontalScrollBar()
        
    def set_dark_mode(self, dark_mode):
        """Update dark mode setting"""
        self.dark_mode = dark_mode
        self.update_style()
        
    def update_style(self):
        """Update the styling based on dark mode"""
        if self.dark_mode:
            self.setStyleSheet("QTextEdit { background-color: #2d2d2d; color: #e0e0e0; }")
        else:
            self.setStyleSheet("QTextEdit { background-color: #ffffff; color: #000000; }")
            
    def set_names(self, names):
        """Set the sequence names to display"""
        if not names:
            self.setPlainText("No sequences")
            return
            
        # Create HTML representation of the sequence names
        html = "<pre style='margin: 0; padding: 0; font-family: Courier New'>"
        
        for name in names:
            # Display up to 30 characters of the name with ellipsis if needed
            if len(name) > 30:
                formatted_name = name[:27] + "..."
            else:
                formatted_name = name
            
            # Ensure minimum width of 14 characters
            formatted_name = formatted_name.ljust(14)
            html += f"{formatted_name}<br>"
        
        html += "</pre>"
        self.setHtml(html)


class SequenceContentViewer(QTextEdit):
    """Widget for displaying sequence content with colored nucleotides"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setFontFamily("Courier New")
        self.setFontPointSize(9)
        self.dark_mode = False
        self.horizontal_scrollbar = self.horizontalScrollBar()
        
    def set_dark_mode(self, dark_mode):
        """Update dark mode setting"""
        self.dark_mode = dark_mode
        self.update_style()
        
    def update_style(self):
        """Update the styling based on dark mode"""
        if self.dark_mode:
            self.setStyleSheet("QTextEdit { background-color: #2d2d2d; color: #e0e0e0; }")
        else:
            self.setStyleSheet("QTextEdit { background-color: #ffffff; color: #000000; }")
            
    def set_sequences(self, sequences, colors):
        """Set the sequences to display"""
        if not sequences:
            self.setPlainText("No sequences")
            return
            
        # Create HTML representation of the sequence content
        html = "<pre style='margin: 0; padding: 0; font-family: Courier New'>"
        
        for sequence in sequences:
            # Add colored nucleotides
            for nucleotide in sequence:
                # Get the proper color for this nucleotide
                nucleotide = nucleotide.upper()
                color = colors.get(nucleotide, QColor(200, 200, 200))  # Default color for unknown nucleotides
                
                # Format as colored text
                html += f"<span style='color: {color.name()}'>{nucleotide}</span>"
            
            html += "<br>"
        
        html += "</pre>"
        self.setHtml(html)


class AlignmentViewer(QTextEdit):
    """Widget for viewing sequence alignments with colored nucleotides"""
    
    def __init__(self, alignment=None, parent=None):
        super().__init__(parent)
        self.alignment = alignment
        self.dark_mode = False
        self.setReadOnly(True)
        self.setFontFamily("Courier New")
        self.setFontPointSize(9)
        
        # Hide vertical scrollbar as it's redundant with sequence list scrollbar
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        
        self.refresh_alignment()
    
    def set_dark_mode(self, dark_mode):
        """Update dark mode setting"""
        self.dark_mode = dark_mode
        self.refresh_alignment()
    
    def set_alignment(self, alignment):
        """Set the alignment to display"""
        self.alignment = alignment
        self.refresh_alignment()
        
    def refresh_alignment(self):
        """Update the alignment display"""
        if not self.alignment: # Catches None or empty list/MSA
            self.setPlainText("No alignment to display")
            return

        # Define colors for nucleotides based on dark/light mode
        if self.dark_mode:
            colors = {
                'A': QColor(255, 100, 100),  # Red
                'C': QColor(100, 100, 255),  # Blue
                'G': QColor(255, 255, 255),  # White instead of black for dark mode
                'T': QColor(100, 255, 100),  # Green
                '-': QColor(150, 150, 150)   # Gray
            }
            self.setStyleSheet("QTextEdit { background-color: #2d2d2d; color: #e0e0e0; }")
        else:
            colors = {
                'A': QColor(255, 0, 0),       # Red
                'C': QColor(0, 0, 255),       # Blue
                'G': QColor(0, 0, 0),         # Black for light mode
                'T': QColor(0, 180, 0),       # Green
                '-': QColor(128, 128, 128)    # Gray
            }
            self.setStyleSheet("QTextEdit { background-color: #ffffff; color: #000000; }")
        
        html = "<pre style='margin: 0; padding: 0; font-family: Courier New'>"
        
        # self.alignment can be a MultipleSeqAlignment or a list of SeqRecords
        records_to_display = []
        if isinstance(self.alignment, MultipleSeqAlignment):
            records_to_display = self.alignment
        elif isinstance(self.alignment, list) and all(isinstance(rec, SeqRecord) for rec in self.alignment):
            records_to_display = self.alignment
        else:
            self.setPlainText("Invalid alignment data type")
            logger.error(f"AlignmentViewer received invalid alignment type: {type(self.alignment)}")
            return

        for record in records_to_display:
            for nucleotide in str(record.seq):
                nucleotide = nucleotide.upper()
                color = colors.get(nucleotide, QColor(200, 200, 200))
                html += f"<span style='color: {color.name()}'>{nucleotide}</span>"
            html += "<br>"
        
        html += "</pre>"
        self.setHtml(html)


class MainWindow(QMainWindow):
    """Main application window"""
    
    # Static variable to keep a strong reference to the debug console
    _debug_console_instance = None
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TreeCraft")
        self.resize(1200, 800)
        
        # Set application icon
        # Use a relative path to the application resources
        self.app_icon = QIcon("resources/treecraft_icon_32px.png")
        self.setWindowIcon(self.app_icon)
        # Set up system tray icon
        self.tray_icon = QSystemTrayIcon(self)
        self.tray_icon.setIcon(self.app_icon)
        self.tray_icon.setToolTip("TreeCraft")
        self.tray_icon.activated.connect(self.tray_icon_activated)
        self.tray_icon.show()
        
        # Set up logging
        self.setup_logging()
        logger = logging.getLogger("treecraft")
        logger.info("Starting TreeCraft application")
        
        # Initialize member variables
        self.sequence_list = SequenceListWidget()
        self.current_alignment = None # Can be MultipleSeqAlignment or list of SeqRecords
        self.current_tree = None
        self.temp_dir = tempfile.mkdtemp(prefix="treecraft_")
        logger.debug(f"Created temporary directory: {self.temp_dir}")
        self.raxml_name_map = {} # For RAxML name handling
        
        # Debug console window (created on demand)
        self.debug_console = None
        
        # Default to dark mode to ensure tree is visible
        self._dark_mode = True  # Use private variable to avoid triggering setter during init
        
        # Create UI components
        self.create_actions()
        self.create_menubar()
        self.create_toolbar()
        self.create_central_widget()
        self.create_statusbar()
        
        # Connect signals
        self.connect_signals()
        
        # Now set the dark mode properly to apply it to all components
        self.dark_mode_action.setChecked(True)
        self.on_toggle_dark_mode(True)
        
        # Show status message
        self.statusbar.showMessage("Ready")
        
        # Check for saved settings
        self.load_settings()
            
    def setup_logging(self):
        """Set up logging"""
        # Logging is now configured at the application level in treecraft.py
        # No need to add handlers here, which could lead to duplicate logs
        pass
        
    def load_settings(self):
        """Load application settings"""
        # TODO: Implement settings loading (theme, etc.)
        pass
        
    def create_actions(self):
        """Create application actions"""
        # File actions
        self.open_action = QAction("&Open Sequences...", self)
        self.open_action.setShortcut("Ctrl+O")
        self.open_action.setStatusTip("Open sequence file")
        
        self.export_sequences_action = QAction("Export Sequences...", self)
        self.export_sequences_action.setStatusTip("Export sequences to file")
        
        self.export_alignment_action = QAction("Export Alignment...", self)
        self.export_alignment_action.setStatusTip("Export aligned sequences to file")
        
        self.export_tree_action = QAction("Export Tree...", self)
        self.export_tree_action.setStatusTip("Export tree to file")
        
        self.exit_action = QAction("E&xit", self)
        self.exit_action.setShortcut("Ctrl+Q")
        self.exit_action.setStatusTip("Exit application")
        self.exit_action.triggered.connect(self.close)
        
        # Edit actions
        self.add_sequence_action = QAction("Add Sequence...", self)
        self.add_sequence_action.setStatusTip("Add a new sequence")
        
        self.edit_sequence_action = QAction("Edit Selected Sequence...", self)
        self.edit_sequence_action.setStatusTip("Edit the selected sequence")
        
        self.delete_sequence_action = QAction("Delete Selected Sequence", self)
        self.delete_sequence_action.setStatusTip("Delete the selected sequence")
        
        # Tree actions
        self.align_action = QAction("Align Sequences...", self)
        self.align_action.setStatusTip("Align sequences")
        
        self.trim_alignment_action = QAction("Trim Alignment...", self)
        self.trim_alignment_action.setStatusTip("Remove poorly aligned regions")
        
        self.build_tree_action = QAction("Build Tree...", self)
        self.build_tree_action.setStatusTip("Build phylogenetic tree")
        
        self.raxml_tree_action = QAction("Build RAxML Tree...", self)
        self.raxml_tree_action.setStatusTip("Build maximum likelihood tree with RAxML")
        
        # View actions
        self.dark_mode_action = QAction("Dark Mode", self)
        self.dark_mode_action.setCheckable(True)
        self.dark_mode_action.setStatusTip("Toggle dark mode")
        
        self.font_action = QAction("Change Font...", self)
        self.font_action.setStatusTip("Change the font used for tree labels")
        self.font_action.triggered.connect(self.on_change_font)
        
        self.sequence_list_action = QAction("Show Sequence List", self)
        self.sequence_list_action.setCheckable(True)
        self.sequence_list_action.setChecked(True)  # Default to visible
        self.sequence_list_action.setStatusTip("Show/hide the sequence list panel")
        
        self.alignment_viewer_action = QAction("Show Alignment Viewer", self)
        self.alignment_viewer_action.setCheckable(True)
        self.alignment_viewer_action.setChecked(True)  # Default to visible
        self.alignment_viewer_action.setStatusTip("Show/hide the alignment viewer panel")
        
        # Help actions
        self.about_action = QAction("&About", self)
        self.about_action.setStatusTip("Show about information")
        
        # Debug console action
        self.debug_console_action = QAction("Debug &Console", self)
        self.debug_console_action.setStatusTip("Show debug console with log messages")

        self.open_tree_action = QAction("&Open Tree File...", self)
        self.open_tree_action.setStatusTip("Open a tree file (Newick, Nexus, PhyloXML)")
        # self.open_tree_action.setShortcut("Ctrl+Shift+O") # Optional shortcut
    
    def create_menubar(self):
        """Create the menu bar"""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(self.open_action)
        file_menu.addAction(self.open_tree_action) # Added here
        file_menu.addSeparator()
        file_menu.addAction(self.export_sequences_action)
        file_menu.addAction(self.export_alignment_action)
        file_menu.addAction(self.export_tree_action)
        file_menu.addSeparator()
        file_menu.addAction(self.exit_action)
        
        # Edit menu
        edit_menu = menubar.addMenu("&Edit")
        edit_menu.addAction(self.add_sequence_action)
        edit_menu.addAction(self.edit_sequence_action)
        edit_menu.addAction(self.delete_sequence_action)
        
        # Tree menu
        tree_menu = menubar.addMenu("&Tree")
        tree_menu.addAction(self.align_action)
        tree_menu.addAction(self.trim_alignment_action)
        tree_menu.addSeparator()
        tree_menu.addAction(self.build_tree_action)
        tree_menu.addAction(self.raxml_tree_action)
        
        # View menu
        view_menu = menubar.addMenu("&View")
        view_menu.addAction(self.dark_mode_action)
        view_menu.addSeparator()
        
        # Line width submenu
        line_width_menu = view_menu.addMenu("Line Width")
        self.line_width_group = QActionGroup(self)
        self.line_width_group.setExclusive(True)
        
        for width in range(1, 6):
            action = QAction(f"{width} px", self)
            action.setCheckable(True)
            if width == 1:  # Default width
                action.setChecked(True)
            action.setData(width)
            action.triggered.connect(self.on_line_width_changed)
            self.line_width_group.addAction(action)
            line_width_menu.addAction(action)
        
        # Font selection
        view_menu.addAction(self.font_action)
        view_menu.addSeparator()
        
        view_menu.addAction(self.sequence_list_action)
        view_menu.addAction(self.alignment_viewer_action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        help_menu.addAction(self.about_action)
        help_menu.addSeparator()
        help_menu.addAction(self.debug_console_action)
    
    def create_toolbar(self):
        """Create the toolbar"""
        toolbar = QToolBar("Main Toolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        
        # Add actions to toolbar
        toolbar.addAction(self.open_action)
        toolbar.addSeparator()
        toolbar.addAction(self.align_action)
        toolbar.addAction(self.build_tree_action)
        toolbar.addSeparator()
        toolbar.addAction(self.export_tree_action)
    
    def create_central_widget(self):
        """Create the central widget"""
        # Main widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Create splitter for resizable panes
        splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left pane - sequence list with hide button
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        
        # Header with title and hide button
        seq_header = QWidget()
        seq_header_layout = QHBoxLayout(seq_header)
        seq_header_layout.setContentsMargins(0, 0, 0, 0)
        
        # Title label with sequence count
        self.seq_label = QLabel("Sequences: 0")
        seq_header_layout.addWidget(self.seq_label)
        
        # Add spacer to push the hide button to the right
        from PyQt6.QtWidgets import QSpacerItem, QSizePolicy
        spacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        seq_header_layout.addItem(spacer)
        
        # Hide button for sequence list
        hide_seq_btn = QPushButton("×")
        hide_seq_btn.setToolTip("Hide Sequence List")
        hide_seq_btn.setFixedSize(20, 20)
        hide_seq_btn.clicked.connect(lambda: self.sequence_list_action.setChecked(False))
        seq_header_layout.addWidget(hide_seq_btn)
        
        left_layout.addWidget(seq_header)
        left_layout.addWidget(self.sequence_list)
        splitter.addWidget(left_widget)
        
        # Middle pane - alignment view with hide button
        middle_widget = QWidget()
        middle_layout = QVBoxLayout(middle_widget)
        
        # Header with title and hide button
        align_header = QWidget()
        align_header_layout = QHBoxLayout(align_header)
        align_header_layout.setContentsMargins(0, 0, 0, 0)
        align_label = QLabel("Alignment:")
        align_header_layout.addWidget(align_label)
        
        # Add spacer to push the hide button to the right
        align_spacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        align_header_layout.addItem(align_spacer)
        
        # Hide button for alignment viewer
        hide_align_btn = QPushButton("×")
        hide_align_btn.setToolTip("Hide Alignment Viewer")
        hide_align_btn.setFixedSize(20, 20)
        hide_align_btn.clicked.connect(lambda: self.alignment_viewer_action.setChecked(False))
        align_header_layout.addWidget(hide_align_btn)
        
        middle_layout.addWidget(align_header)
        self.alignment_viewer = AlignmentViewer()
        # NOTE: We're now using the scrollbar for bidirectional sync, so keep it visible
        # This allows users to scroll either the sequence list or alignment pane
        # and have them stay in sync
        middle_layout.addWidget(self.alignment_viewer)
        splitter.addWidget(middle_widget)
        
        # Right pane - tree view
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.addWidget(QLabel("Phylogenetic Tree:"))
        self.tree_canvas = PhyloCanvas()
        right_layout.addWidget(self.tree_canvas)
        splitter.addWidget(right_widget)
        
        # Set initial sizes for the panes
        splitter.setSizes([300, 400, 500])
    
    def create_statusbar(self):
        """Create the status bar"""
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
    
    def connect_signals(self):
        """Connect signals to slots"""
        # File menu
        self.open_action.triggered.connect(self.on_open)
        self.open_tree_action.triggered.connect(self.on_open_tree_file) # Added here
        self.export_sequences_action.triggered.connect(self.on_export_sequences)
        self.export_alignment_action.triggered.connect(self.on_export_alignment)
        self.export_tree_action.triggered.connect(self.on_export_tree)
        
        # Edit menu
        self.add_sequence_action.triggered.connect(self.on_add_sequence)
        self.edit_sequence_action.triggered.connect(self.on_edit_sequence)
        self.delete_sequence_action.triggered.connect(self.on_delete_sequence)
        
        # View menu
        self.dark_mode_action.triggered.connect(self.on_toggle_dark_mode)
        self.sequence_list_action.toggled.connect(self.toggle_sequence_list)
        self.alignment_viewer_action.toggled.connect(self.toggle_alignment_viewer)
        
        # Tree menu
        self.align_action.triggered.connect(self.on_align_sequences)
        self.trim_alignment_action.triggered.connect(self.on_trim_alignment)
        self.build_tree_action.triggered.connect(self.on_build_tree)
        self.raxml_tree_action.triggered.connect(self.on_build_raxml_tree)
        
        # Help menu
        self.about_action.triggered.connect(self.show_about)
        self.debug_console_action.triggered.connect(self.show_debug_console)
        
        # Sequence list signals
        self.sequence_list.itemDoubleClicked.connect(self.on_sequence_double_clicked)
        self.sequence_list.rename_signal.connect(self.on_edit_sequence)
        self.sequence_list.delete_signal.connect(self.on_delete_sequence)
        
        # Synchronize the scrolling between sequence list and alignment viewer
        self.sync_sequence_alignment_scrollbars()
        
        # Tree canvas signals
        self.tree_canvas.node_double_clicked.connect(self.rename_tree_node)
        self.tree_canvas.node_right_clicked.connect(self.show_node_context_menu)
        self.tree_canvas.rename_signal.connect(self.rename_tree_node)
        self.tree_canvas.delete_signal.connect(self.delete_tree_node)
        self.tree_canvas.update_status.connect(self.statusbar.showMessage)
    
    def on_toggle_dark_mode(self, checked):
        """Toggle dark/light mode"""
        self.dark_mode = checked
        logger = logging.getLogger("treecraft")
        logger.info(f"Dark mode: {self.dark_mode}")
        
        # Ensure dark mode action is checked/unchecked to match state
        self.dark_mode_action.setChecked(checked)
        
        # Apply dark mode to components
        if self.dark_mode:
            self.setStyleSheet("""
                QWidget { background-color: #2d2d2d; color: #e0e0e0; }
                QMenuBar { background-color: #353535; color: #e0e0e0; }
                QMenuBar::item:selected { background-color: #5a5a5a; }
                QMenu { background-color: #353535; color: #e0e0e0; }
                QMenu::item:selected { background-color: #5a5a5a; }
                QToolBar { background-color: #353535; color: #e0e0e0; }
                QStatusBar { background-color: #353535; color: #e0e0e0; }
                QListWidget { background-color: #2d2d2d; color: #e0e0e0; border: 1px solid #5a5a5a; }
                QListWidget::item:selected { background-color: #4a4a4a; }
                QSplitter::handle { background-color: #5a5a5a; }
                QLabel { color: #e0e0e0; }
                QMessageBox { background-color: #2d2d2d; color: #e0e0e0; }
                QPushButton { background-color: #353535; color: #e0e0e0; border: 1px solid #5a5a5a; border-radius: 2px; padding: 5px; }
                QPushButton:hover { background-color: #454545; }
                QPushButton:pressed { background-color: #5a5a5a; }
                QComboBox { background-color: #353535; color: #e0e0e0; border: 1px solid #5a5a5a; border-radius: 2px; padding: 1px 18px 1px 3px; min-width: 6em; }
                QComboBox:editable { background-color: #353535; }
                QComboBox::drop-down { border-color: #5a5a5a; }
                QComboBox QAbstractItemView { background-color: #353535; color: #e0e0e0; selection-background-color: #5a5a5a; }
                QLineEdit { background-color: #353535; color: #e0e0e0; border: 1px solid #5a5a5a; border-radius: 2px; padding: 2px; }
                QTextEdit { background-color: #2d2d2d; color: #e0e0e0; border: 1px solid #5a5a5a; }
                QCheckBox { color: #e0e0e0; }
                QRadioButton { color: #e0e0e0; }
                QSpinBox { background-color: #353535; color: #e0e0e0; border: 1px solid #5a5a5a; border-radius: 2px; }
                QTabWidget::pane { border: 1px solid #5a5a5a; }
                QTabBar::tab { background-color: #353535; color: #e0e0e0; border: 1px solid #5a5a5a; border-bottom-color: none; border-top-left-radius: 4px; border-top-right-radius: 4px; padding: 5px; }
                QTabBar::tab:selected { background-color: #2d2d2d; border-bottom-color: #2d2d2d; }
                QScrollBar:vertical { background-color: #2d2d2d; width: 12px; margin: 16px 0 16px 0; }
                QScrollBar::handle:vertical { background-color: #5a5a5a; min-height: 20px; border-radius: 3px; }
                QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { background: none; }
                QScrollBar:horizontal { background-color: #2d2d2d; height: 12px; margin: 0 16px 0 16px; }
                QScrollBar::handle:horizontal { background-color: #5a5a5a; min-width: 20px; border-radius: 3px; }
                QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { background: none; }
            """)
        else:
            self.setStyleSheet("")  # Reset to default style
        
        # Update other components
        self.alignment_viewer.set_dark_mode(self.dark_mode)
        self.tree_canvas.set_dark_mode(self.dark_mode)
    
    def on_line_width_changed(self):
        """Change the line width for tree branches"""
        # Get the width from the sender action
        action = self.sender()
        if action and action.data() is not None:
            width = action.data()
            logger = logging.getLogger("treecraft")
            logger.info(f"Changing tree line width to {width}px")
            
            # Update the tree canvas with new line width
            self.tree_canvas.set_branch_width(width)
            
            # Update status
            self.statusbar.showMessage(f"Tree branch width set to {width}px", 3000)
    
    def on_change_font(self):
        """Change the font used for tree labels"""
        logger = logging.getLogger("treecraft")
        logger.info("Opening font selection dialog")
        
        # Get current font from tree canvas
        current_font = self.tree_canvas.sequence_font if hasattr(self.tree_canvas, 'sequence_font') else QFont("Arial", 9)
        
        # Show font dialog
        from PyQt6.QtWidgets import QFontDialog
        font, ok = QFontDialog.getFont(current_font, self, "Select Tree Label Font")
        
        if ok:
            # Update font in tree canvas
            logger.info(f"Setting tree font to {font.family()}, {font.pointSize()}pt")
            self.tree_canvas.sequence_font = font
            
            # Update status
            self.statusbar.showMessage(f"Tree label font set to {font.family()}, {font.pointSize()}pt", 3000)
        
    def toggle_sequence_list(self, checked):
        """Toggle visibility of the sequence list panel"""
        logger = logging.getLogger("treecraft")
        logger.info(f"Sequence list visible: {checked}")
        
        # Find the left widget (sequence list) in the splitter
        central_widget = self.centralWidget()
        splitter = central_widget.layout().itemAt(0).widget()
        
        # The left widget (index 0) is the sequence list widget
        left_widget = splitter.widget(0)
        
        # Toggle visibility
        if left_widget:
            left_widget.setVisible(checked)
            
        # Update splitter sizes to distribute space properly
        if checked:
            # Get current visibility of alignment viewer
            alignment_visible = self.alignment_viewer_action.isChecked()
            width = splitter.width()
            
            if alignment_visible:
                # Both sequence list and alignment viewer are visible
                splitter.setSizes([width//3, width//3, width//3])
            else:
                # Only sequence list and tree view are visible
                splitter.setSizes([width//2, 0, width//2])
        else:
            # Sequence list is hidden
            width = splitter.width()
            alignment_visible = self.alignment_viewer_action.isChecked()
            
            if alignment_visible:
                # Only alignment viewer and tree view are visible
                splitter.setSizes([0, width//2, width//2])
            else:
                # Only tree view is visible
                splitter.setSizes([0, 0, width])
    
    def toggle_alignment_viewer(self, checked):
        """Toggle visibility of the alignment viewer panel"""
        logger = logging.getLogger("treecraft")
        logger.info(f"Alignment viewer visible: {checked}")
        
        # Find the middle widget (alignment viewer) in the splitter
        central_widget = self.centralWidget()
        splitter = central_widget.layout().itemAt(0).widget()
        
        # The middle widget (index 1) is the alignment viewer widget
        middle_widget = splitter.widget(1)
        
        # Toggle visibility
        if middle_widget:
            middle_widget.setVisible(checked)
            
        # Update splitter sizes to distribute space properly
        if checked:
            # Get current visibility of sequence list
            sequence_list_visible = self.sequence_list_action.isChecked()
            width = splitter.width()
            
            if sequence_list_visible:
                # Both sequence list and alignment viewer are visible
                splitter.setSizes([width//3, width//3, width//3])
            else:
                # Only alignment viewer and tree view are visible
                splitter.setSizes([0, width//2, width//2])
        else:
            # Alignment viewer is hidden
            width = splitter.width()
            sequence_list_visible = self.sequence_list_action.isChecked()
            
            if sequence_list_visible:
                # Only sequence list and tree view are visible
                splitter.setSizes([width//2, 0, width//2])
            else:
                # Only tree view is visible
                splitter.setSizes([0, 0, width])
    
    def parse_malformed_fasta(self, file_path):
        """
        Custom parser for FASTA files that can handle malformed files with sequences 
        without proper line breaks (e.g. multiple sequences smooshed together)
        """
        logger = logging.getLogger("treecraft")
        records = []
        
        try:
            # Read the file
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Split the content at all '>' characters
            parts = content.split('>')
            
            # Process each part (the first one might be empty if the file starts with '>')
            for part in parts:
                if not part.strip():
                    continue
                
                # Split the part into lines
                lines = part.splitlines()
                
                # First line is the header
                header = lines[0].strip()
                
                # Parse the ID and description from the header
                header_parts = header.split(None, 1)
                if len(header_parts) >= 1:
                    seq_id = header_parts[0]
                    # Use FULL header as description for proper display in list
                    desc = header  # This ensures the full header will be the description
                else:
                    # Skip if we don't have a proper header
                    logger.warning(f"Skipping segment with invalid header: {header[:20]}...")
                    continue
                
                # Join the rest as the sequence, checking for embedded sequences
                seq_part = ""
                for i, line in enumerate(lines[1:]):
                    if '>' in line:
                        # Found an embedded sequence marker
                        pos = line.find('>')
                        
                        # If there's sequence data before the '>', add it to the current sequence
                        if pos > 0:
                            seq_part += line[:pos].strip()
                        
                        # Create a record for the current sequence
                        if seq_part:
                            try:
                                record = SeqRecord(Seq(seq_part), id=seq_id, description=desc)
                                records.append(record)
                                logger.info(f"Added sequence: {seq_id} (length: {len(seq_part)})")
                            except Exception as e:
                                logger.error(f"Error creating record for {seq_id}: {e}")
                        
                        # The rest of this line plus all remaining lines form a new embedded sequence
                        # Handle it by re-splitting and processing
                        new_content = line[pos:] + '\n' + '\n'.join(lines[i+2:])
                        logger.warning(f"Found embedded sequence at line {i+2}. Processing as new segment.")
                        
                        # Process the embedded segments by recursively calling the function
                        # on the modified content
                        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
                            tmp.write(new_content)
                            tmp_path = tmp.name
                        
                        try:
                            # Try to parse the embedded content
                            embedded_records = list(SeqIO.parse(tmp_path, "fasta"))
                            if embedded_records:
                                records.extend(embedded_records)
                                logger.info(f"Added {len(embedded_records)} embedded sequences")
                            
                            # Clean up temporary file
                            os.unlink(tmp_path)
                        except Exception as e:
                            logger.error(f"Error parsing embedded sequence: {e}")
                            os.unlink(tmp_path)
                        
                        # We've handled the embedded sequence, so break out of this loop
                        break
                    else:
                        # Regular sequence line
                        seq_part += line.strip()
                else:
                    # No embedded '>' found, add the complete sequence
                    if seq_part:
                        try:
                            record = SeqRecord(Seq(seq_part), id=seq_id, description=desc)
                            records.append(record)
                            logger.info(f"Added normal sequence: {seq_id} (length: {len(seq_part)})")
                        except Exception as e:
                            logger.error(f"Error creating record for {seq_id}: {e}")
            
            return records
                
        except Exception as e:
            logger.error(f"Error in custom FASTA parser: {e}")
            # If our custom parser fails, fall back to the standard parser
            logger.info("Falling back to standard Bio.SeqIO parser")
            
            try:
                return list(SeqIO.parse(file_path, "fasta"))
            except Exception as inner_e:
                logger.error(f"Standard parser also failed: {inner_e}")
                # Last resort - our own simple fix by replacing embedded '>' with newline + '>'
                try:
                    with open(file_path, 'r') as f:
                        content = f.read()
                    
                    # Replace embedded '>' with a newline plus '>'
                    fixed_content = re.sub(r'([ACGTN])>', r'\1\n>', content, flags=re.IGNORECASE)
                    
                    # Create a temporary file with the fixed content
                    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
                        tmp.write(fixed_content)
                        tmp_path = tmp.name
                    
                    # Try parsing the fixed content
                    fixed_records = list(SeqIO.parse(tmp_path, "fasta"))
                    
                    # Clean up
                    os.unlink(tmp_path)
                    
                    if fixed_records:
                        logger.info(f"Last resort parser succeeded with {len(fixed_records)} sequences")
                        return fixed_records
                    else:
                        return []
                except Exception as last_e:
                    logger.error(f"Last resort parser also failed: {last_e}")
                    return []
    
    def on_open(self):
        """Open a sequence file"""
        # Show file dialog
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Sequence Files (*.fasta *.fa *.fastq *.fq *.gb *.gbk *.genbank *.phy *.phylip);;All Files (*)")
        file_dialog.setViewMode(QFileDialog.ViewMode.List)
        
        if not file_dialog.exec():
            return
            
        selected_files = file_dialog.selectedFiles()
        if not selected_files:
            return
            
        file_path = selected_files[0]
        
        # Try to detect format from file extension
        format_name = "fasta"  # Default format
        ext = os.path.splitext(file_path)[1].lower()
        
        if ext in (".fq", ".fastq"):
            format_name = "fastq"
        elif ext in (".gb", ".gbk", ".genbank"):
            format_name = "genbank"
        elif ext in (".phy", ".phylip"):
            format_name = "phylip"
        
        # Load sequences
        try:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor) # Set wait cursor
            logger = logging.getLogger("treecraft")
            logger.info(f"Loading sequences from file: {file_path}")
            
            # For FASTA files, use our custom parser that can handle malformed files
            if format_name == "fasta":
                records = self.parse_malformed_fasta(file_path)
                logger.info(f"Used custom parser for FASTA file")
            else:
                records = list(SeqIO.parse(file_path, format_name))
            
            if not records:
                QMessageBox.warning(self, "Warning", "No sequences found in the file")
                # QApplication.restoreOverrideCursor() # Ensure cursor is restored - already in finally
                return # Return is important here, before finally if we want to avoid other operations
                
            logger.info(f"Loaded {len(records)} sequences")
            
            # Get existing sequence count for reporting
            existing_count = len(self.sequence_list.sequences)
            
            # Do NOT clear existing sequences - append new ones instead
            duplicate_count = 0
            new_count = 0
            
            for record in records:
                # Extract the sequence string and description
                sequence = str(record.seq)
                
                # Get a proper description and avoid duplicates
                sequence_id = record.id.strip()  # Use the ID as the sequence_id
                
                # For the display name, ALWAYS use the full description
                if hasattr(record, 'description') and record.description:
                    # Check for duplicate ID in the description and fix it
                    if record.description.startswith(record.id) and record.description.count(record.id) > 1:
                        # Remove duplicated IDs from description
                        parts = record.description.split()
                        seen = set()
                        unique_parts = []
                        for part in parts:
                            if part not in seen:
                                seen.add(part)
                                unique_parts.append(part)
                        desc = " ".join(unique_parts)
                        logger.debug(f"Fixed duplicate ID in description: {record.description} -> {desc}")
                    else:
                        # Use the full description as-is
                        desc = record.description
                else:
                    # If no description, use the ID as the display name too
                    desc = record.id
                    
                # Make sure description is never shorter than the ID
                if len(desc) < len(sequence_id):
                    desc = sequence_id
                
                # Check if this sequence already exists
                if sequence_id in self.sequence_list.sequences:
                    duplicate_count += 1
                    # Skip duplicate sequences or optionally replace them
                    # For now, we'll just skip them to preserve existing sequences
                    logger.debug(f"Skipping duplicate sequence: {sequence_id}")
                    continue
                    
                # Add to our sequence list - the ID and display name should be the same initially
                # unless the record had a different description
                self.sequence_list.add_sequence(sequence_id, sequence, desc)
                new_count += 1
            
            # Update status with counts of new and duplicate sequences
            if duplicate_count > 0:
                self.statusbar.showMessage(
                    f"Added {new_count} new sequences from {os.path.basename(file_path)} "
                    f"({duplicate_count} duplicates skipped). Total: {existing_count + new_count} sequences."
                )
            else:
                self.statusbar.showMessage(
                    f"Added {new_count} sequences from {os.path.basename(file_path)}. "
                    f"Total: {existing_count + new_count} sequences."
                )
            
            # Update sequence count in header
            self.update_sequence_count()
            
            # Create unaligned records to display in the alignment viewer
            self.update_alignment_display_from_sequences()
            
            # Clear tree if we added any new sequences
            if new_count > 0:
                self.current_tree = None
                self.tree_canvas.set_tree(None) # This will internally handle layout_is_dirty
            
        except Exception as e:
            logger.error(f"Error loading sequences: {e}")
            QMessageBox.critical(self, "Error", f"Failed to load sequences: {str(e)}")
        finally:
            QApplication.restoreOverrideCursor() # Restore cursor in finally block

    def on_open_tree_file(self):
        logger = logging.getLogger("treecraft")
        logger.info("Attempting to open a tree file...")

        # Ensure necessary imports are within the method or at class/module level
        # from PyQt6.QtWidgets import QFileDialog, QApplication, QMessageBox # Already imported
        # from PyQt6.QtCore import Qt # Already imported
        # from Bio import Phylo # Imported locally below
        # import os # Already imported

        dialog = QFileDialog(self)
        dialog.setWindowTitle("Open Tree File")
        dialog.setNameFilter("Tree Files (*.nwk *.newick *.nex *.nexus *.xml *.phyloxml);;Newick (*.nwk *.newick);;Nexus (*.nex *.nexus);;PhyloXML (*.xml *.phyloxml);;All Files (*)")
        dialog.setViewMode(QFileDialog.ViewMode.List)
        dialog.setFileMode(QFileDialog.FileMode.ExistingFile)

        tree_object = None
        tree_file_path_for_status = "Unknown" # For status bar message

        if dialog.exec():
            selected_files = dialog.selectedFiles()
            if selected_files:
                tree_file_path = selected_files[0]
                tree_file_path_for_status = os.path.basename(tree_file_path)
                file_ext = os.path.splitext(tree_file_path)[1].lower()

                format_hint = "newick" # Default
                if file_ext in (".nex", ".nexus"):
                    format_hint = "nexus"
                elif file_ext in (".xml", ".phyloxml"):
                    format_hint = "phyloxml"

                logger.info(f"User selected tree file: {tree_file_path} (format hint: {format_hint})")

                try:
                    QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                    from Bio import Phylo
                    tree_object = Phylo.read(tree_file_path, format_hint)

                    if tree_object:
                        # Process quoted names
                        for clade in tree_object.get_terminals() + list(tree_object.get_nonterminals()):
                            if clade.name:
                                if clade.name.startswith("'") and clade.name.endswith("'"):
                                    clade.name = clade.name[1:-1]
                                elif clade.name.startswith('"') and clade.name.endswith('"'):
                                    clade.name = clade.name[1:-1]
                        logger.info(f"Successfully parsed tree from {tree_file_path_for_status}")
                except Exception as e:
                    logger.error(f"Failed to load or parse tree file '{tree_file_path}': {e}")
                    QMessageBox.critical(self, "Error Opening Tree", f"Failed to load tree file: {str(e)}")
                    tree_object = None
                finally:
                    QApplication.restoreOverrideCursor()

        if tree_object:
            self.current_tree = tree_object

            logger.info("Clearing existing sequence list and alignment for new tree display.")
            if hasattr(self.sequence_list, 'clear_all_sequences') and callable(getattr(self.sequence_list, 'clear_all_sequences')):
                self.sequence_list.clear_all_sequences()
            else:
                logger.info("Manually clearing sequence_list (no clear_all_sequences method found).")
                self.sequence_list.sequences.clear()
                self.sequence_list.descriptions.clear()
                self.sequence_list.clear()

            self.current_alignment = None
            if self.alignment_viewer:
                self.alignment_viewer.set_alignment(None)

            self.update_sequence_count() # Update "Sequences: 0" label

            self.tree_canvas.set_tree(self.current_tree)
            self.statusbar.showMessage(f"Opened tree file: {tree_file_path_for_status}", 5000)
            logger.info(f"Tree '{tree_file_path_for_status}' set to display.")
        else:
            # Only show "failed" if a file was actually selected.
            if dialog.result() == QDialog.DialogCode.Accepted and tree_file_path_for_status != "Unknown":
                self.statusbar.showMessage("Failed to open or parse tree file.", 3000)
            else: # Dialog was cancelled or no file selected
                self.statusbar.showMessage("Open tree file cancelled.", 3000)
            logger.info("Tree loading failed, was cancelled, or no file selected.")
            
    def update_sequence_count(self):
        """Update the sequence count in the header"""
        count = len(self.sequence_list.sequences)
        self.seq_label.setText(f"Sequences: {count}")
    
    def update_alignment_display_from_sequences(self):
        """Create and display alignment from current sequences (even if unaligned)"""
        if not self.sequence_list.sequences:
            self.current_alignment = None
            self.alignment_viewer.set_alignment(None)
            return
        
        # Get logger instance
        import logging # Ensure logger is available
        logger = logging.getLogger("treecraft") # Ensure logger is available
            
        # Extract the sequence order from the list widget to match the sequence list display
        sequence_ids_in_order = []
        for i in range(self.sequence_list.count()):
            item = self.sequence_list.item(i)
            if item:
                # Get the sequence ID from the item's user data
                seq_id = item.data(Qt.ItemDataRole.UserRole)
                if seq_id:
                    sequence_ids_in_order.append(seq_id)
        
        # Log the sequence order for debugging
        logger.debug(f"Sequence order for alignment display: {sequence_ids_in_order}")
        
        # Create SeqRecord objects from the current sequences in the same order as the list widget
        records = []
        
        # First add sequences in the order they appear in the list widget
        for seq_id in sequence_ids_in_order:
            if seq_id in self.sequence_list.sequences:
                sequence = self.sequence_list.sequences[seq_id]
                description = self.sequence_list.descriptions.get(seq_id, seq_id)
                record = SeqRecord(Seq(sequence), id=seq_id, description=description)
                records.append(record)
        
        # Add any remaining sequences not found in the list widget
        for seq_id, sequence in self.sequence_list.sequences.items():
            if seq_id not in sequence_ids_in_order:
                description = self.sequence_list.descriptions.get(seq_id, seq_id)
                record = SeqRecord(Seq(sequence), id=seq_id, description=description)
                records.append(record)

        # <<< LOGGING START >>>
        logger.debug(f"Entering update_alignment_display_from_sequences. Number of records: {len(records)}")
        if not records: # This check is technically redundant due to the one at the very top of the function
            logger.debug("Record list is empty (second check).")
            self.current_alignment = None
            self.alignment_viewer.set_alignment(None)
            QTimer.singleShot(100, lambda: self.sync_alignment_scroll_position())
            return
        else:
            for i, r in enumerate(records):
                logger.debug(f"Record {i}: ID='{r.id}', Length={len(r.seq)}, Description='{r.description}'")
                # logger.debug(f"Record {i} Seq Sample: {r.seq[:30]}...") # Optional: log sequence sample
        # <<< LOGGING END >>>

        lengths = [len(record.seq) for record in records]
        all_same_length = all(l == lengths[0] for l in lengths) if lengths else True # Keep `if lengths else True` for safety, though `if not records` above should catch empty.

        # <<< LOGGING START >>>
        logger.debug(f"Calculated lengths: {lengths}")
        logger.debug(f"Result of all_same_length check: {all_same_length}")
        # <<< LOGGING END >>>

        if all_same_length:
            # <<< LOGGING START >>>
            logger.debug("All sequences reported as same length. Attempting to create MultipleSeqAlignment object.")
            # <<< LOGGING END >>>
            alignment = MultipleSeqAlignment(records)
            self.current_alignment = alignment
            self.alignment_viewer.set_alignment(alignment)
            # logger.debug("Created MultipleSeqAlignment as all sequences are of the same length.") # Redundant with above
        else:
            logger.info("Displaying sequences in alignment view; sequences are not of uniform length and won't be treated as a formal alignment.")
            self.current_alignment = records
            self.alignment_viewer.set_alignment(records)

        # <<< LOGGING START >>>
        if self.current_alignment is None:
            logger.debug("self.current_alignment is None after update_alignment_display.")
        else:
            logger.debug(f"self.current_alignment is now type: {type(self.current_alignment)}, Number of sequences in current_alignment: {len(self.current_alignment) if self.current_alignment is not None else 0}")
        # <<< LOGGING END >>>

        # Reset the vertical scroll position to match the sequence list
        QTimer.singleShot(100, lambda: self.sync_alignment_scroll_position())

    def sync_alignment_scroll_position(self):
        """Force sync the alignment position to match the sequence list"""
        # If we're already in a synchronizing operation, avoid recursion
        if hasattr(self, 'is_synchronizing') and self.is_synchronizing:
            return
            
        try:
            # Set synchronizing flag to avoid recursion
            self.is_synchronizing = True
            
            # Get the sequence list scrollbar
            seq_scrollbar = self.sequence_list.verticalScrollBar()
            if seq_scrollbar:
                # Find first visible row in sequence list
                first_visible_row = None
                for i in range(self.sequence_list.count()):
                    item_rect = self.sequence_list.visualItemRect(self.sequence_list.item(i))
                    if item_rect.bottom() >= 0 and item_rect.top() < self.sequence_list.height():
                        first_visible_row = i
                        break
                
                if first_visible_row is not None:
                    # Calculate line height in alignment view
                    alignment_line_height = self.alignment_viewer.fontMetrics().height()
                    
                    # Calculate new position for alignment based on first visible sequence
                    new_pos = first_visible_row * alignment_line_height
                    
                    # Update the alignment viewer scroll position
                    self.alignment_viewer.verticalScrollBar().setValue(new_pos)
                    
                    # Log the sync operation
                    logger = logging.getLogger("treecraft")
                    logger.debug(f"Forced sync to row {first_visible_row}, alignment pos {new_pos}")
        finally:
            # Clear synchronizing flag
            self.is_synchronizing = False
    
    def on_add_sequence(self):
        """Add a new sequence"""
        dialog = SequenceDialog(self)
        if dialog.exec():
            seq_id, sequence, description = dialog.get_sequence_info()
            
            if seq_id and sequence:
                # Add to sequence list
                self.sequence_list.add_sequence(seq_id, sequence, description)
                
                # Update status
                self.statusbar.showMessage(f"Added sequence: {seq_id}")
                
                # Update sequence count
                self.update_sequence_count()
                
                # Update alignment display with the new sequence
                self.update_alignment_display_from_sequences()
                
                # Clear tree since it's now outdated
                self.current_tree = None
                self.tree_canvas.set_tree(None)
    
    def on_edit_sequence(self):
        """Edit the selected sequence"""
        selected_items = self.sequence_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Warning", "Please select a sequence to edit")
            return
            
        selected_item = selected_items[0]
        seq_id = selected_item.data(Qt.ItemDataRole.UserRole)
        
        if seq_id and seq_id in self.sequence_list.sequences:
            sequence = self.sequence_list.sequences[seq_id]
            description = self.sequence_list.descriptions.get(seq_id, seq_id)
            
            dialog = EditSequenceDialog(parent=self, sequence_id=seq_id, sequence=sequence, description=description)
            if dialog.exec():
                _, new_seq_id, new_sequence, new_description = dialog.get_sequence_data()
                
                if new_seq_id and new_sequence:
                    # Update sequence in the list
                    self.sequence_list.update_sequence(seq_id, new_seq_id, new_sequence, new_description)
                    
                    # Update status
                    self.statusbar.showMessage(f"Updated sequence: {new_seq_id}")
                    
                    # Update sequence count
                    self.update_sequence_count()
                    
                    # Update alignment display with the modified sequence
                    self.update_alignment_display_from_sequences()
                    
                    # Clear tree since it's now outdated
                    self.current_tree = None
                    self.tree_canvas.set_tree(None)
    
    def on_delete_sequence(self):
        """Delete the selected sequences"""
        selected_items = self.sequence_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Warning", "Please select at least one sequence to delete")
            return
            
        # Check if we need to delete multiple sequences
        if len(selected_items) > 1:
            # Confirm deletion of multiple sequences
            reply = QMessageBox.question(
                self,
                "Confirm Delete",
                f"Are you sure you want to delete {len(selected_items)} selected sequences?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # Collect sequence IDs to delete
                sequence_ids_to_delete = []
                for item in selected_items:
                    seq_id = item.data(Qt.ItemDataRole.UserRole)
                    if seq_id and seq_id in self.sequence_list.sequences:
                        sequence_ids_to_delete.append(seq_id)
                
                # Delete all selected sequences
                for seq_id in sequence_ids_to_delete:
                    self.sequence_list.remove_sequence(seq_id)
                
                # Update status
                self.statusbar.showMessage(f"Deleted {len(sequence_ids_to_delete)} sequences")
                
                # Update sequence count
                self.update_sequence_count()
                
                # Update alignment view with the remaining sequences
                self.update_alignment_display_from_sequences()
                    
                # If there's a tree and we have enough sequences, rebuild it
                if self.current_tree and len(self.sequence_list.sequences) >= 2:
                    # Simple rebuild using UPGMA for now
                    params = {"method": "upgma", "distance_model": "identity"}
                    self.current_tree = build_tree_with_params(self.current_alignment, params)
                    self.tree_canvas.set_tree(self.current_tree)
                else:
                    # Not enough sequences left for a tree
                    self.current_tree = None
                    self.tree_canvas.set_tree(None)
        else:
            # Single sequence deletion (original behavior)
            selected_item = selected_items[0]
            seq_id = selected_item.data(Qt.ItemDataRole.UserRole)
            
            if seq_id and seq_id in self.sequence_list.sequences:
                # Confirm deletion
                reply = QMessageBox.question(
                    self,
                    "Confirm Delete",
                    f"Are you sure you want to delete sequence '{seq_id}'?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No
                )
                
                if reply == QMessageBox.StandardButton.Yes:
                    # Remove from sequence list
                    description = self.sequence_list.descriptions.get(seq_id, seq_id)
                    self.sequence_list.remove_sequence(seq_id)
                    
                    # Update status
                    self.statusbar.showMessage(f"Deleted sequence: {seq_id}")
                    
                    # Update sequence count
                    self.update_sequence_count()
                    
                    # Update alignment view with the remaining sequences
                    self.update_alignment_display_from_sequences()
                        
                    # If there's a tree and we have enough sequences, rebuild it
                    if self.current_tree and len(self.sequence_list.sequences) >= 2:
                        # Simple rebuild using UPGMA for now
                        params = {"method": "upgma", "distance_model": "identity"}
                        self.current_tree = build_tree_with_params(self.current_alignment, params)
                        self.tree_canvas.set_tree(self.current_tree)
                    else:
                        # Not enough sequences left for a tree
                        self.current_tree = None
                        self.tree_canvas.set_tree(None)
    
    def on_sequence_double_clicked(self, item):
        """Handle double click on a sequence"""
        # Get the sequence ID and edit it
        seq_id = item.data(Qt.ItemDataRole.UserRole)
        if seq_id and seq_id in self.sequence_list.sequences:
            sequence = self.sequence_list.sequences[seq_id]
            description = self.sequence_list.descriptions.get(seq_id, seq_id)
            
            dialog = EditSequenceDialog(parent=self, sequence_id=seq_id, sequence=sequence, description=description)
            if dialog.exec():
                _, new_seq_id, new_sequence, new_description = dialog.get_sequence_data()
                
                if new_seq_id and new_sequence:
                    # Update sequence in the list
                    self.sequence_list.update_sequence(seq_id, new_seq_id, new_sequence, new_description)
                    
                    # Update status
                    self.statusbar.showMessage(f"Updated sequence: {new_seq_id}")
                    
                    # Update sequence count
                    self.update_sequence_count()
                    
                    # Update alignment display with the modified sequence
                    self.update_alignment_display_from_sequences()
                    
                    # Clear tree since it's now outdated
                    self.current_tree = None
                    self.tree_canvas.set_tree(None)
    
    def on_align_sequences(self):
        """Align sequences"""
        # Check if we have sequences
        if not self.sequence_list.sequences or len(self.sequence_list.sequences) < 2:
            QMessageBox.warning(self, "Warning", "Need at least two sequences to align")
            return
            
        # Create Bio.SeqRecord objects from the sequences
        records = []
        for seq_id, sequence in self.sequence_list.sequences.items():
            description = self.sequence_list.descriptions.get(seq_id, seq_id)
            record = SeqRecord(Seq(sequence), id=seq_id, description=description)
            records.append(record)
            
        # Show alignment dialog
        dialog = SequenceAlignmentDialog(self, records)
        if dialog.exec():
            # Get the aligned sequences
            aligned_records = dialog.get_aligned_sequences()
            
            if aligned_records:
                self.current_alignment = aligned_records
                self.alignment_viewer.set_alignment(aligned_records)
                
                # Scroll to the middle of the alignment horizontally
                scrollbar = self.alignment_viewer.horizontalScrollBar()
                if scrollbar:
                    # Scroll to half of the maximum value
                    QTimer.singleShot(100, lambda: scrollbar.setValue(scrollbar.maximum() // 2))
                
                # Update status
                self.statusbar.showMessage(f"Aligned {len(aligned_records)} sequences")
                
                # Clear any existing tree
                self.current_tree = None
                self.tree_canvas.set_tree(None)
    
    def on_trim_alignment(self):
        """Trim poorly aligned regions"""
        if not self.current_alignment or len(self.current_alignment) < 2:
            QMessageBox.warning(self, "Warning", "Need aligned sequences to trim")
            return
            
        # Ask user which trimming tool to use
        tools = ["trimAl", "Gblocks"]
        tool, ok = QInputDialog.getItem(
            self, "Trim Alignment", "Choose trimming method:", tools, 0, False
        )
        
        if not ok:
            return
            
        # Create temporary files
        import tempfile
        input_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        input_file.close()
        output_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        output_file.close()
        
        # Write alignment to input file
        with open(input_file.name, "w") as f:
            for record in self.current_alignment:
                f.write(f">{record.id}\n{record.seq}\n")
        
        # Set up progress dialog
        progress_dialog = QMessageBox(self)
        progress_dialog.setWindowTitle("Trimming Alignment")
        progress_dialog.setText("Trimming alignment... please wait.")
        progress_dialog.setStandardButtons(QMessageBox.StandardButton.NoButton)
        progress_dialog.setModal(True)
        progress_dialog.show()
        QApplication.processEvents()
        
        try:
            # Run the appropriate trimming tool
            if tool == "trimAl":
                success, info = run_trimal(input_file.name, output_file.name, 0.5, 0.6)
            else:  # Gblocks
                success, info = run_gblocks(input_file.name, output_file.name)
                
            if success:
                # Load the trimmed alignment
                trimmed_records = list(SeqIO.parse(output_file.name, "fasta"))
                
                if trimmed_records:
                    # Preserve the descriptions from the original alignment
                    desc_map = {record.id: record.description for record in self.current_alignment}
                    
                    for record in trimmed_records:
                        if record.id in desc_map:
                            record.description = desc_map[record.id]
                    
                    self.current_alignment = trimmed_records
                    self.alignment_viewer.set_alignment(trimmed_records)
                    
                    # Show results message
                    if "percent_kept" in info:
                        QMessageBox.information(
                            self,
                            "Trimming Results",
                            f"Alignment trimmed: {info['percent_kept']}% of positions kept\n"
                            f"Original: {info['original_count']} positions\n"
                            f"Trimmed: {info['trimmed_count']} positions"
                        )
                    else:
                        QMessageBox.information(self, "Trimming Results", info.get("message", "Alignment trimmed successfully"))
                    
                    # Update status
                    self.statusbar.showMessage(f"Trimmed alignment to {len(trimmed_records[0].seq)} positions")
                    
                    # Since the alignment changed, clear any existing tree
                    self.current_tree = None
                    self.tree_canvas.set_tree(None)
                    
                else:
                    QMessageBox.warning(self, "Warning", "No sequences in trimmed alignment")
            else:
                QMessageBox.critical(self, "Error", f"Trimming failed: {info.get('message', 'Unknown error')}")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error during trimming: {str(e)}")
            
        finally:
            # Close progress dialog
            progress_dialog.close()
            
            # Clean up temporary files
            try:
                os.unlink(input_file.name)
                os.unlink(output_file.name)
            except:
                pass
    
    def on_build_tree(self):
        """Build a phylogenetic tree"""
        if not self.current_alignment or len(self.current_alignment) < 2:
            # Check if we have sequences that need aligning
            if self.sequence_list.sequences and len(self.sequence_list.sequences) >= 2:
                reply = QMessageBox.question(
                    self,
                    "Align First?",
                    "Sequences need to be aligned before building a tree. Align now?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.Yes
                )
                
                if reply == QMessageBox.StandardButton.Yes:
                    self.on_align_sequences()
                    # If alignment was canceled or failed
                    if not self.current_alignment or len(self.current_alignment) < 2:
                        return
                else:
                    return
            else:
                QMessageBox.warning(self, "Warning", "Need at least two aligned sequences to build a tree")
                return
        
        # Show tree building dialog
        dialog = TreeBuildDialog(self)
        if dialog.exec():
            params = dialog.get_parameters()
            
            # Set up progress dialog with a proper QProgressDialog
            from PyQt6.QtWidgets import QProgressDialog
            progress_dialog = QProgressDialog("Building phylogenetic tree...", "Cancel", 0, 100, self)
            progress_dialog.setWindowTitle("Building Tree")
            progress_dialog.setWindowModality(Qt.WindowModality.WindowModal)
            progress_dialog.setValue(10)
            progress_dialog.setMinimumDuration(0)
            progress_dialog.show()
            QApplication.processEvents()
            
            try:
                # Build the tree with selected parameters
                tree = build_tree_with_params(self.current_alignment, params)
                
                if tree:
                    self.current_tree = tree
                    self.tree_canvas.set_tree(tree)
                    
                    # Update status
                    self.statusbar.showMessage(f"Built phylogenetic tree with {len(self.current_alignment)} sequences")
                else:
                    QMessageBox.warning(self, "Warning", "Failed to build tree")
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error building tree: {str(e)}")
                
            finally:
                # Close progress dialog
                progress_dialog.close()
    
    def on_build_raxml_tree(self):
        """Build a tree using RAxML"""
        # Initialize logger at the method level
        logger = logging.getLogger("treecraft")
        if not self.current_alignment or len(self.current_alignment) < 3:
            # RAxML requires at least 3 sequences
            if self.sequence_list.sequences and len(self.sequence_list.sequences) >= 3:
                reply = QMessageBox.question(
                    self,
                    "Align First?",
                    "Sequences need to be aligned before building a tree. Align now?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.Yes
                )
                
                if reply == QMessageBox.StandardButton.Yes:
                    self.on_align_sequences()
                    # If alignment was canceled or failed
                    if not self.current_alignment or len(self.current_alignment) < 3:
                        return
                else:
                    return
            else:
                QMessageBox.warning(self, "Warning", "Need at least three aligned sequences to build a tree with RAxML")
                return
        
        # Show RAxML dialog
        dialog = RaxMLDialog(self)
        if dialog.exec():
            params = dialog.get_parameters()
            
            # Set up progress dialog
            from PyQt6.QtWidgets import QProgressDialog
            progress = QProgressDialog("Building Maximum Likelihood Tree...", "Cancel", 0, 100, self)
            progress.setWindowTitle("RAxML Tree")
            progress.setMinimumDuration(0)
            progress.setWindowModality(Qt.WindowModality.WindowModal)
            progress.setValue(10)
            
            try:
                # Create temporary directory for RAxML
                import tempfile
                temp_dir = tempfile.mkdtemp(prefix="raxml_")
                
                # --- RAxML Name Handling: Pre-processing ---
                self.raxml_name_map = {} # Ensure it's fresh for this run (using self.raxml_name_map consistently)
                sanitized_records = []
                alignment_records = list(self.current_alignment) if isinstance(self.current_alignment, MultipleSeqAlignment) else self.current_alignment

                logger.info(f"Preparing alignment for RAxML. Original alignment has {len(alignment_records)} records.")

                for i, record in enumerate(alignment_records):
                    original_full_name = record.description if record.description and record.description.strip() else record.id
                    sanitized_id = f"s{i}" # Use simple 's0', 's1', ... style IDs
                    self.raxml_name_map[sanitized_id] = original_full_name
                    # Write FASTA with only the sanitized ID and empty description for RAxML
                    sanitized_records.append(SeqRecord(record.seq, id=sanitized_id, description=""))
                    # logger.debug(f"Mapping: '{sanitized_id}' -> '{original_full_name}'") # Logged below

                # Log first 5 entries of the map for debugging
                logged_map_entries = 0
                for k, v in self.raxml_name_map.items():
                    logger.debug(f"RAxML Name Map: '{k}' -> '{v}'")
                    logged_map_entries += 1
                    if logged_map_entries >= 5:
                        break
                if len(self.raxml_name_map) > 5:
                    logger.debug(f"... and {len(self.raxml_name_map) - 5} more entries.")

                aligned_file_path = os.path.join(temp_dir, "alignment_sanitized_for_raxml.fasta")
                SeqIO.write(sanitized_records, aligned_file_path, "fasta")
                logger.info(f"RAxML pre-processing: Wrote {len(sanitized_records)} records with sanitized names (e.g., s0, s1...) to {aligned_file_path}")
                # --- End RAxML Name Handling: Pre-processing ---
                
                # Run RAxML with progress dialog
                result_info = run_raxml(temp_dir, aligned_file_path, params, progress) # Pass the sanitized file
                
                if result_info:
                    tree_dir = result_info.get("raxml_working_dir", temp_dir)
                    tree_file = find_tree_file(tree_dir)
                    
                    if tree_file and os.path.exists(tree_file):
                        # Log raw tree sample
                        try:
                            with open(tree_file, 'r') as f:
                                raw_tree_sample = f.read(200)
                            logger.debug(f"RAxML raw tree sample: {raw_tree_sample.strip()}...")
                        except Exception as e:
                            logger.error(f"Could not read sample from tree file '{tree_file}': {e}")

                        from Bio import Phylo
                        tree = Phylo.read(tree_file, "newick")
                        
                        # --- RAxML Name Handling: Post-processing ---
                        # The map from result_info might be based on RAxML's internal sanitization if it differs.
                        # Forcing use of self.raxml_name_map as the primary source of truth.
                        # raxml_output_name_map = result_info.get("sequence_name_map", {})
                        # if raxml_output_name_map:
                        #     logger.info(f"RAxML's run_raxml returned a name map with {len(raxml_output_name_map)} entries. Prioritizing self.raxml_name_map.")
                            # Example: if RAxML changed 's1' to 'T1', raxml_output_name_map might be {'T1': 's1'}
                            # We need to map from what's in the tree file (e.g. 'T1') back to our original name.
                            # This requires careful chaining if RAxML internally renames our 'sX' IDs.
                            # For now, assume RAxML uses 'sX' as-is due to their simplicity.

                        logger.info(f"RAxML post-processing: Restoring original names using self.raxml_name_map ({len(self.raxml_name_map)} entries).")
                        restored_count = 0
                        not_found_count = 0
                        problematic_names_from_tree = []

                        for terminal_node in tree.get_terminals(): # Use a different variable name here
                            name_from_tree = terminal_node.name
                            logger.debug(f"Processing terminal from tree: '{name_from_tree}'")

                            original_name = self.raxml_name_map.get(name_from_tree)

                            if original_name:
                                terminal_node.name = original_name
                                restored_count += 1
                                logger.debug(f"Restored RAxML name: '{name_from_tree}' -> '{original_name}'")
                            else:
                                not_found_count +=1
                                problematic_names_from_tree.append(name_from_tree)
                                logger.warning(f"RAxML post-processing: Could not find original name for ID '{name_from_tree}' in self.raxml_name_map.")

                        logger.info(f"RAxML name restoration: {restored_count} names restored, {not_found_count} not found in map.")
                        if not_found_count > 0:
                             logger.warning(f"Problematic names from RAxML tree file (not found in map): {problematic_names_from_tree}")
                             logger.warning(f"RAxML Name Map dump for problematic tree (self.raxml_name_map has {len(self.raxml_name_map)} entries): First few: {dict(list(self.raxml_name_map.items())[:5])}")
                        # --- End RAxML Name Handling: Post-processing ---
                        
                        self.current_tree = tree
                        # Store that this is a RaxML tree to help with sequence matching
                        self.tree_source = 'raxml'
                        self.tree_canvas.set_tree(tree)
                        
                        # Collect detailed information for the success dialog
                        raxml_binary = "Unknown"
                        if "raxml_direct_path" in params:
                            raxml_binary = params["raxml_direct_path"]
                        elif "raxml_binary" in params:
                            raxml_binary = params["raxml_binary"]
                        
                        # Get execution time if available
                        runtime = "Unknown"
                        if "runtime" in result_info:
                            runtime = f"{result_info['runtime']:.2f} seconds"
                        
                        # Count nodes in tree
                        terminal_count = len(list(tree.get_terminals()))
                        internal_count = len(list(tree.get_nonterminals()))
                        
                        # Create a rich detailed message
                        detail_message = f"""
RAxML TREE BUILDING COMPLETED SUCCESSFULLY

COMMAND:
{result_info.get('command', 'Command information not available')}

TREE INFORMATION:
• Number of terminal nodes (leaves): {terminal_count}
• Number of internal nodes: {internal_count}
• Total nodes: {terminal_count + internal_count}
"""
                        
                        # Add log-likelihood if available
                        if "log_likelihood" in result_info:
                            detail_message += f"• Log-likelihood: {result_info['log_likelihood']}\n"

                        detail_message += f"""
PERFORMANCE:
• Execution time: {runtime}

PARAMETERS USED:
• Evolutionary model: {params.get('model', 'GTR+GAMMA')}
• RAxML binary: {raxml_binary}
• Number of CPU threads: {params.get('threads', '8')}
• Random seed: {params.get('seed', '12345')}
"""

                        # Add bootstrap info if applicable
                        if params.get("bootstrap", False):
                            detail_message += f"• Bootstrap replicates: {params.get('bootstrap_replicates', 100)}\n"
                            
                        # Add additional info if available
                        if "tree_file" in result_info:
                            detail_message += f"• Tree file: {result_info['tree_file']}\n"
                            
                        # Add sequence information
                        detail_message += f"\nSEQUENCE INFORMATION:\n• Number of sequences: {len(self.current_alignment)}\n"
                        
                        # Add alignment length if we can determine it
                        if self.current_alignment and len(self.current_alignment) > 0:
                            align_length = len(self.current_alignment[0].seq)
                            detail_message += f"• Alignment length: {align_length} sites\n"
                        
                        # Use different titles based on bootstrap
                        if params.get("bootstrap", False):
                            title = "RAxML Maximum Likelihood Tree with Bootstrap Support"
                        else:
                            title = "RAxML Maximum Likelihood Tree"
                            
                        # Show the enhanced success message
                        QMessageBox.information(
                            self,
                            title,
                            detail_message
                        )
                        
                        # Update status
                        self.statusbar.showMessage(f"Built RAxML tree with {len(self.current_alignment)} sequences")
                    else:
                        QMessageBox.warning(self, "Warning", "Failed to find output tree file")
                else:
                    error_msg = "Failed to build tree with RAxML"
                    if isinstance(result_info, dict) and "error" in result_info:
                        error_msg += f": {result_info['error']}"
                    QMessageBox.critical(self, "Error", error_msg)
                    
            except Exception as e:
                import traceback
                error_details = traceback.format_exc()
                QMessageBox.critical(
                    self, 
                    "Error", 
                    f"Error building RAxML tree: {str(e)}\n\nDetails: {error_details}"
                )
                
            finally:
                # Close progress dialog
                progress.close()
    
    def on_export_tree(self):
        """Export the tree to a file"""
        if not self.current_tree:
            QMessageBox.warning(self, "Warning", "No tree to export")
            return
            
        # Show export dialog
        dialog = ExportDialog(self, self.current_tree)
        if dialog.exec():
            # Get selected format
            format_ext = dialog.get_selected_format()
            
            # Determine file type based on format
            if format_ext in ("png", "jpg", "svg"):
                # For image formats, show image options dialog
                img_options_dialog = ImageExportOptionsDialog(self)
                if img_options_dialog.exec():
                    options = img_options_dialog.get_options()
                    # Get save path
                    file_path, _ = QFileDialog.getSaveFileName(
                        self, "Export Tree Image", "", f"Image Files (*.{format_ext})"
                    )
                    if file_path:
                        # Ensure file has the correct extension
                        if not file_path.lower().endswith(f".{format_ext}"):
                            file_path = f"{file_path}.{format_ext}"
                            
                        # Export as image using PhyloCanvas
                        success = self.tree_canvas.export_image(file_path, options)
                        
                        if success:
                            self.statusbar.showMessage(f"Tree exported to {file_path}", 3000)
                        else:
                            QMessageBox.critical(self, "Export Error", f"Failed to export tree to {file_path}")
            else:
                # For text-based formats (nwk, nex, xml)
                file_path, _ = QFileDialog.getSaveFileName(
                    self, "Export Tree", "", 
                    f"Tree Files (*.{format_ext})"
                )
                if file_path:
                    # Ensure file has the correct extension
                    if not file_path.lower().endswith(f".{format_ext}"):
                        file_path = f"{file_path}.{format_ext}"
                    
                    # Export tree in selected format
                    from Bio import Phylo
                    try:
                        # Create a deep copy of the tree to avoid modifying the original
                        import copy
                        tree_copy = copy.deepcopy(self.current_tree)
                        
                        # Write the tree to the file
                        Phylo.write(tree_copy, file_path, format_ext)
                        
                        # Verify the file was created and has content
                        import os
                        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                            self.statusbar.showMessage(f"Tree exported to {file_path}", 5000)
                        else:
                            QMessageBox.critical(self, "Export Error", f"Failed to export tree: File is empty")
                    except Exception as e:
                        QMessageBox.critical(self, "Export Error", f"Failed to export tree: {e}")
    
    def show_node_context_menu(self, node_name, point):
        """Show context menu for tree node"""
        import logging
        from PyQt6.QtCore import QTimer
        logger = logging.getLogger("treecraft.main_window")
        logger.debug(f"Showing context menu for node: {node_name}")
        
        menu = QMenu(self)
        rename_action = menu.addAction("Rename")
        delete_action = menu.addAction("Delete")
        
        # Add additional options if this is an internal node
        internal_node = False
        for clade in self.current_tree.get_nonterminals():
            if hasattr(clade, 'name') and clade.name == node_name:
                internal_node = True
                break
                
        if internal_node:
            menu.addSeparator()
            highlight_action = menu.addAction("Highlight Clade")
            collapse_action = menu.addAction("Collapse Clade")
            
            # Handle highlight/collapse
            highlight_action.triggered.connect(lambda: self.highlight_tree_clade(node_name))
            collapse_action.triggered.connect(lambda: self.collapse_tree_clade(node_name))
        
        # Connect actions directly with captured node name
        rename_action.triggered.connect(lambda checked=False, name=node_name: self.rename_tree_node(name))
        delete_action.triggered.connect(lambda checked=False, name=node_name: self.delete_tree_node(name))
        
        # Show the menu at the point and block until it returns
        action = menu.exec(self.tree_canvas.mapToGlobal(point))
        logger.debug(f"Context menu action: {action.text() if action else 'None'}")
        
        # Prevent tree dragging operations after menu
        if hasattr(self.tree_canvas, 'tree_canvas'):
            # Access the actual TreeCanvas if we're using a PhyloCanvas wrapper
            canvas = self.tree_canvas.tree_canvas
        else:
            canvas = self.tree_canvas
            
        # Reset dragging state on the actual canvas
        canvas.dragging = False
        canvas.panning = False
        canvas.drag_mode = None
        canvas.drag_clade = None
    
    def highlight_tree_clade(self, node_name):
        """Highlight a clade in the tree"""
        # Find the clade by name
        clade = None
        for node in self.current_tree.find_clades():
            if hasattr(node, 'name') and node.name == node_name:
                clade = node
                break
                
        if not clade:
            return
            
        # Show color dialog
        color = QColorDialog.getColor(QColor(255, 0, 0), self, "Select Highlight Color")
        if not color.isValid():
            return
            
        # Apply the color to the clade
        self.tree_canvas.highlight_clade(clade, color.name())
    
    def collapse_tree_clade(self, node_name):
        """Collapse a clade in the tree"""
        # Find the clade by name
        clade = None
        for node in self.current_tree.find_clades():
            if hasattr(node, 'name') and node.name == node_name:
                clade = node
                break
                
        if not clade:
            return
            
        # Toggle collapsed state
        self.tree_canvas.toggle_collapse_clade(clade)
        
        # Update display
        self.tree_canvas.set_tree(self.current_tree)
    
    def rename_tree_node(self, node_name):
        """Rename a sequence in the tree via double-click"""
        logger = logging.getLogger("treecraft")
        logger.debug(f"Rename request for node: {node_name}")
        
        # Check if this node exists in our sequences
        sequence_id = None
        
        # Special handling for RAxML tree nodes that may have been truncated
        # Example: "dirty" should be expanded to "MK281405 Inocybe 'lanatopurpurea PNW30' JLF5893 iNat66020010 OR (dirty)"
        original_node_name = node_name
        
        # For nodes named like "type NY", look for any full name that contains this fragment
        full_match_found = False
        if ' ' in node_name:  # Only try this for multi-word fragments
            for seq_id, desc in self.sequence_list.descriptions.items():
                if node_name in desc:
                    logger.info(f"Found full name '{desc}' containing '{node_name}'")
                    node_name = desc  # Use the full name 
                    full_match_found = True
                    break
        
        # Only proceed with pattern cleanup if we didn't find a full match
        if not full_match_found:
            # Fix duplicate names issue (like "KX897432.1_KX897432.1")
            # Detect and clean up self-repeating patterns in the node name
            if node_name:
                # Check for duplicated node names like KX897432.1_KX897432.1
                parts = node_name.split('_')
                if len(parts) > 1 and parts[0] == parts[1]:
                    logger.warning(f"Detected exact duplicated name parts in node name: {node_name}")
                    node_name = parts[0]  # Use just the first part
                    logger.info(f"Simplified duplicated name to: {node_name}")
                
                # Check for repeated space-separated segments 
                elif ' ' in node_name:
                    words = node_name.split()
                    # Check if same word appears multiple times (e.g., "KX897432.1 KX897432.1")
                    if len(set(words)) < len(words):
                        logger.warning(f"Detected repeated words in: {node_name}")
                        # Keep only unique parts to build the clean name
                        clean_name = ' '.join(dict.fromkeys(words))
                        logger.info(f"Cleaned name from '{node_name}' to '{clean_name}'")
                        node_name = clean_name
        
        # First, check if node_name matches a sequence ID directly
        if node_name in self.sequence_list.sequences:
            sequence_id = node_name
            logger.debug(f"Found direct ID match: {node_name}")
        else:
            # Try more matching approaches
            found = False
            
            # Try matching by description
            for seq_id, desc in self.sequence_list.descriptions.items():
                if desc == node_name or seq_id == node_name:
                    sequence_id = seq_id
                    found = True
                    logger.debug(f"Found exact match by description: {desc} -> {seq_id}")
                    break
                    
            # If not found, try matching just the first part (before space)
            if not found and ' ' in node_name:
                first_part = node_name.split()[0]
                logger.debug(f"Trying to match by first part: {first_part}")
                
                # Keep track of all potential matches
                potential_matches = []
                
                for seq_id, desc in self.sequence_list.descriptions.items():
                    # Check if this description starts with the first part or the ID does
                    if desc.startswith(first_part) or seq_id.startswith(first_part):
                        potential_matches.append((seq_id, desc))
                
                # If we found exactly one match, use it
                if len(potential_matches) == 1:
                    sequence_id = potential_matches[0][0]
                    found = True
                    logger.info(f"Found unique match by first part: {first_part} -> {sequence_id}")
                # If we found multiple matches, log them and use the first one
                elif len(potential_matches) > 1:
                    logger.warning(f"Found multiple matches by first part '{first_part}': {potential_matches}")
                    # Use the first match
                    sequence_id = potential_matches[0][0]
                    found = True
                    logger.info(f"Using first match: {sequence_id}")
            
            # If still not found, try a more fuzzy match based on contained strings
            if not found:
                logger.debug(f"Trying more fuzzy matches for: {node_name}")
                
                for seq_id, desc in self.sequence_list.descriptions.items():
                    # Check if the node name is contained within the description or vice versa
                    if (node_name in desc) or (desc in node_name):
                        sequence_id = seq_id
                        found = True
                        logger.info(f"Found fuzzy match: '{node_name}' and '{desc}' -> {seq_id}")
                        break
        
        # Check if we're dealing with a RaxML tree
        is_raxml_tree = False
        if hasattr(self, 'tree_source') and self.tree_source == 'raxml':
            is_raxml_tree = True
            logger.info("RaxML tree detected from tree_source attribute")
            
        # If we couldn't find the sequence ID but we're dealing with a RaxML tree
        if not sequence_id and is_raxml_tree:
            # Extract any full name that might be embedded in the node_name
            # Example: from "type NY" find the full name "KX897432.1 Inocybe griseoscabrosa PQ (type NY)"
            full_name = ""
            for desc in self.sequence_list.descriptions.values():
                if node_name in desc:
                    full_name = desc
                    logger.info(f"Found full name '{full_name}' containing '{node_name}'")
                    break
            
            # For RaxML trees, create a temporary ID based on the node name
            # This allows renaming even when we can't find the exact sequence
            logger.info(f"RaxML tree special handling: Using node name as temporary ID: {node_name}")
            sequence_id = f"raxml_node_{node_name}"
            
            # Add this as a new sequence entry in the list (or use an existing item with similar name)
            if not any(node_name in desc for desc in self.sequence_list.descriptions.values()):
                # Try to find a sequence with a prefix or suffix match
                found_id = None
                for seq_id, desc in self.sequence_list.descriptions.items():
                    if node_name in desc or desc in node_name:
                        found_id = seq_id
                        full_name = desc  # Use the matched description as full name
                        logger.info(f"Found similar sequence for RaxML node: {desc}")
                        break
                
                if found_id:
                    # Use the existing entry
                    sequence_id = found_id
                else:
                    # If not found at all, create a placeholder entry with the best name we have
                    display_name = full_name if full_name else node_name
                    logger.info(f"Creating placeholder entry for RaxML node: {display_name}")
                    sequence_text = ""  # Empty sequence since we don't have the actual sequence
                    self.sequence_list.add_sequence(sequence_id, sequence_text, display_name)
                    self.sequence_list.rebuild_sequence_list()
        elif not sequence_id:
            # Only show warning if it's not a RaxML tree
            QMessageBox.warning(self, "Rename Failed", f"Could not find sequence '{node_name}' in the current dataset.")
            return
            
        # Get the current description - ensure we get the full description
        current_desc = self.sequence_list.descriptions.get(sequence_id, sequence_id)
        logger.debug(f"Showing rename dialog for sequence '{sequence_id}' with description '{current_desc}'")
        
        # Use a more robust dialog approach - create the dialog but don't show it yet
        try:
            from PyQt6.QtWidgets import QInputDialog
            dialog = QInputDialog(self)
            dialog.setWindowTitle("Rename Sequence")
            dialog.setLabelText("Enter new name for the sequence:")
            dialog.setTextValue(current_desc)
            # Make the dialog wider to accommodate longer sequence names
            dialog.resize(500, dialog.height())
            
            # Connect to finished signal to avoid issues with double-processing
            dialog.finished.connect(lambda result: self.handle_rename_dialog_result(result, dialog, sequence_id))
            
            # Show the dialog non-modally
            dialog.open()
            
            # Return now to avoid the rest of the code executing immediately
            return
        except Exception as e:
            logger.error(f"Error showing rename dialog: {e}")
            QMessageBox.critical(self, "Error", f"Could not show rename dialog: {str(e)}")
            return
    
    def handle_rename_dialog_result(self, result, dialog, sequence_id):
        """Handle the result of a sequence rename dialog"""
        logger = logging.getLogger("treecraft")
        
        # Check if the dialog was accepted
        if result == QDialog.DialogCode.Accepted:
            # Get the new name from the dialog
            new_name = dialog.textValue()
            
            if new_name and sequence_id in self.sequence_list.sequences:
                # Update the description only (keeping the same sequence_id)
                logger.debug(f"Updating sequence description for {sequence_id} to: {new_name}")
                self.sequence_list.descriptions[sequence_id] = new_name
                
                # Update the list widget
                self.sequence_list.rebuild_sequence_list()
                
                # Update the current tree
                if self.current_tree:
                    self.update_tree_node_names(sequence_id, new_name)
                
                # Update status
                self.statusbar.showMessage(f"Renamed sequence '{sequence_id}' to '{new_name}'", 3000)
                
                # Update the alignment viewer
                self.update_alignment_display_from_sequences()
        else:
            logger.debug("Rename dialog cancelled")
    
    def update_tree_node_names(self, sequence_id, new_name):
        """Update tree node names when a sequence is renamed"""
        logger = logging.getLogger("treecraft")
        
        if not self.current_tree:
            return
            
        # Find all terminal nodes that might match this sequence
        renamed_count = 0
        for terminal in self.current_tree.get_terminals():
            # Check for various forms of the name
            terminal_name = terminal.name if terminal.name else ""
            
            # Check if this is the node we want to rename
            is_match = False
            
            # Direct match with sequence_id
            if terminal_name == sequence_id:
                is_match = True
            # Match with old description in duplicated format
            elif '_' in terminal_name and terminal_name.startswith(sequence_id + '_'):
                is_match = True
            # Check for prefix match
            elif terminal_name.startswith(sequence_id + ' '):
                is_match = True
                
            if is_match:
                # Update the node name with the new display name
                old_name = terminal.name
                terminal.name = new_name
                renamed_count += 1
                logger.info(f"Renamed tree node from '{old_name}' to '{new_name}'")
                
        if renamed_count > 0:
            # Names within self.current_tree (which is self.tree_canvas.tree_canvas.tree) are updated.
            # We need to trigger a repaint of the TreeCanvas, not a full set_tree.
            # The TreeCanvas's draw_clade will pick up the new names.
            # No need to set layout_is_dirty as positions haven't changed.
            if hasattr(self.tree_canvas, 'tree_canvas'):
                self.tree_canvas.tree_canvas.update() # Call update on the actual TreeCanvas instance
            else:
                self.tree_canvas.update() # Fallback if direct tree_canvas not found (should not happen)
            logger.info(f"Updated {renamed_count} node names in tree and triggered TreeCanvas repaint.")
    
    def delete_tree_node(self, node_name):
        """Delete a sequence from the tree and dataset via right-click"""
        logger = logging.getLogger("treecraft")
        logger.debug(f"Delete request for node: {node_name}")

        sequence_id_to_delete = None

        # Primary Search: Match node_name (expected to be the full display name/description)
        # with the values in self.sequence_list.descriptions.
        if self.sequence_list.descriptions:
            for seq_id, desc in self.sequence_list.descriptions.items():
                if desc == node_name:
                    sequence_id_to_delete = seq_id
                    logger.debug(f"Found sequence by description match: '{desc}' -> ID: '{seq_id}'")
                    break
        
        # Secondary Search: If not found by description, try matching node_name
        # with the sequence IDs (keys in self.sequence_list.sequences).
        if sequence_id_to_delete is None:
            if node_name in self.sequence_list.sequences:
                sequence_id_to_delete = node_name
                logger.debug(f"Found sequence by ID match: '{node_name}'")

        if sequence_id_to_delete is None:
            # If still not found after both primary and secondary search
            # This block replaces the complex fallback logic that was previously here.
            # The assumption is that node_name *should* be the full description.
            logger.warning(f"Delete node: Could not find sequence corresponding to displayed name '{node_name}' in sequence list using primary (description) or secondary (ID) search.")
            QMessageBox.warning(self, "Delete Failed", f"Could not find sequence data for node '{node_name}'. The sequence may have been modified or the name is inconsistent.")
            return

        # Get the description for confirmation using the found sequence_id_to_delete
        current_desc = self.sequence_list.descriptions.get(sequence_id_to_delete, sequence_id_to_delete)
        
        # Check if we're in delete mode from phylo_canvas
        in_delete_mode = hasattr(self.tree_canvas, 'delete_mode') and self.tree_canvas.delete_mode
        
        # Class attribute to store the "Don't ask again" preference
        if not hasattr(MainWindow, 'skip_delete_confirmation'):
            MainWindow.skip_delete_confirmation = False
        
        # In delete mode, we ALWAYS skip the confirmation dialog
        # If not in delete mode, show confirmation ONLY if "Don't ask again" is not checked
        if in_delete_mode:
            logger.debug(f"In delete mode - automatically confirming deletion without dialog")
            reply = QMessageBox.StandardButton.Yes
        elif MainWindow.skip_delete_confirmation:
            logger.debug(f"'Don't ask again' is checked - skipping confirmation dialog")
            reply = QMessageBox.StandardButton.Yes
        else:
            # Create a custom dialog with "Don't ask again" checkbox
            from PyQt6.QtWidgets import QCheckBox, QVBoxLayout, QHBoxLayout, QDialog, QPushButton, QLabel
            
            dialog = QDialog(self)
            dialog.setWindowTitle("Confirm Delete")
            
            layout = QVBoxLayout(dialog)
            
            # Message
            message = QLabel(f"Are you sure you want to delete sequence '{current_desc}'?")
            layout.addWidget(message)
            
            # "Don't ask again" checkbox
            dont_ask_checkbox = QCheckBox("Don't ask again")
            layout.addWidget(dont_ask_checkbox)
            
            # Buttons
            button_layout = QHBoxLayout()
            yes_button = QPushButton("Yes")
            no_button = QPushButton("No")
            button_layout.addWidget(yes_button)
            button_layout.addWidget(no_button)
            layout.addLayout(button_layout)
            
            # Connect buttons
            yes_button.clicked.connect(lambda: dialog.done(QMessageBox.StandardButton.Yes))
            no_button.clicked.connect(lambda: dialog.done(QMessageBox.StandardButton.No))
            
            # Execute dialog
            reply = dialog.exec()
            
            # Store "Don't ask again" preference if user clicked Yes
            if reply == QMessageBox.StandardButton.Yes and dont_ask_checkbox.isChecked():
                logger.debug("User checked 'Don't ask again', storing preference")
                MainWindow.skip_delete_confirmation = True
        
        if reply == QMessageBox.StandardButton.Yes:
            logger.info(f"Deleting sequence: ID={sequence_id_to_delete}, Description={current_desc}")
            
            # Remove from sequence list
            if sequence_id_to_delete in self.sequence_list.sequences:
                del self.sequence_list.sequences[sequence_id_to_delete]
            if sequence_id_to_delete in self.sequence_list.descriptions:
                del self.sequence_list.descriptions[sequence_id_to_delete]
                
            # Update the list widget
            self.sequence_list.rebuild_sequence_list()
            
            # Update sequence count
            self.update_sequence_count()
            
            # Remove from the alignment if it exists
            if self.current_alignment:
                # Find ONLY the exact record that matches this sequence ID
                # This prevents deleting sequences with similar names
                new_alignment = []
                deleted_record = None
                
                for record in self.current_alignment:
                    # Use EXACT match for sequence_id_to_delete to avoid deleting similar sequences
                    if record.id == sequence_id_to_delete:
                        logger.info(f"Removing exact match record from alignment: {record.id}")
                        deleted_record = record
                    # Special case for when record.id is formatted differently from sequence_id_to_delete
                    elif hasattr(record, 'description') and record.description == current_desc: # current_desc is from the node to be deleted
                        logger.info(f"Removing description match record from alignment: {record.id} (desc: {record.description})")
                        deleted_record = record
                    else:
                        new_alignment.append(record)
                
                # Ensure we only removed ONE record
                if deleted_record is not None:
                    logger.info(f"Successfully removed record with ID: {deleted_record.id}")
                    # If no records were removed, that's an issue worth logging
                    if len(new_alignment) == len(self.current_alignment):
                        logger.warning(f"No records were removed from alignment for sequence ID: {sequence_id_to_delete}")
                else:
                    logger.warning(f"No matching record found in alignment for sequence ID: {sequence_id_to_delete}")
                
                self.current_alignment = new_alignment
                
                # Update the alignment display
                self.alignment_viewer.set_alignment(new_alignment)
                
            # Force the tree canvas to exit delete mode if it's active
            if hasattr(self.tree_canvas, 'delete_mode') and self.tree_canvas.delete_mode:
                logger.debug("Disabling tree canvas delete mode after successful deletion")
                self.tree_canvas.delete_mode = False
                # Reset the cursor to default
                self.tree_canvas.setCursor(Qt.CursorShape.ArrowCursor)
            
            # Rebuild the tree with better error handling
            try:
                if self.current_alignment and len(self.current_alignment) >= 2:
                    # Simple rebuild using current parameters
                    params = {"method": "upgma", "distance_model": "identity"}
                    logger.info(f"Rebuilding tree with {len(self.current_alignment)} sequences")
                    self.current_tree = build_tree_with_params(self.current_alignment, params)
                    
                    # Check if tree build was successful
                    if self.current_tree:
                        logger.info("Tree rebuilt successfully, updating display")
                        self.tree_canvas.set_tree(self.current_tree)
                    else:
                        logger.error("Failed to rebuild tree")
                        QMessageBox.warning(self, "Warning", "Failed to rebuild tree after sequence deletion.")
                        self.current_tree = None
                        self.tree_canvas.set_tree(None)
                else:
                    # Not enough sequences left for a tree
                    logger.info("Not enough sequences left for a tree")
                    self.current_tree = None
                    self.tree_canvas.set_tree(None)
            except Exception as e:
                logger.error(f"Error rebuilding tree after deletion: {str(e)}")
                import traceback
                logger.error(traceback.format_exc())
                QMessageBox.warning(self, "Error", f"Error rebuilding tree: {str(e)}")
                self.current_tree = None
                self.tree_canvas.set_tree(None)
                
            # Update alignment display after deletion
            self.update_alignment_display_from_sequences()
            
            # Update status
            self.statusbar.showMessage(f"Deleted sequence '{current_desc}'")
    
    def show_about(self):
        """Show about information"""
        QMessageBox.about(
            self,
            "About TreeCraft",
            "TreeCraft - Phylogenetic Tree Visualization Tool\n\n"
            "Version 1.0 by Alan Rockefeller\n\n"
            "April 21, 2025\n\n"
            "https://github.com/AlanRockefeller/Treecraft\n\n"
            "TreeCraft allows you to create, visualize, and manipulate "
            "phylogenetic trees with an easy-to-use interface.\n\n"
            "Features include:\n"
            "• Loading FASTA and FASTQ sequence files\n"
            "• Building trees with multiple methods\n"
            "• Rearranging tree visualization\n"
            "• Zoom and pan controls with mouse wheel support\n"
            "• Highlighting branches and adding labels\n"
            "• Exporting trees in various formats"
        )
        
    def show_debug_console(self):
        """Show the debug console window"""
        # Create debug console if it doesn't exist
        if not self.debug_console:
            self.debug_console = DebugConsole(self)
            # Keep a strong reference to prevent garbage collection
            MainWindow._debug_console_instance = self.debug_console
            
        # Show the console
        self.debug_console.show()
        self.debug_console.raise_()
        self.debug_console.activateWindow()
        
    def sync_sequence_alignment_scrollbars(self):
        """Synchronize the scrolling between sequence list and alignment viewer"""
        logger = logging.getLogger("treecraft")
        logger.info("Setting up one-way scrollbar synchronization from sequence list to alignment viewer")
        
        # Get scrollbar from sequence list widget (alignment viewer's scrollbar is hidden)
        sequence_vscroll = self.sequence_list.verticalScrollBar()
        
        if not sequence_vscroll:
            logger.warning("Could not find sequence list vertical scrollbar")
            return
            
        # Flag to prevent infinite recursion between scroll events
        self.is_synchronizing = False
        
        # Function to find the first visible sequence list item
        def find_first_visible_row():
            first_visible_row = None
            for i in range(self.sequence_list.count()):
                item_rect = self.sequence_list.visualItemRect(self.sequence_list.item(i))
                # Item is visible if its top edge is within or above viewport and bottom edge is within or below
                if item_rect.bottom() >= 0 and item_rect.top() < self.sequence_list.height():
                    first_visible_row = i
                    break
            return first_visible_row
        
        # Function to sync alignment content position with sequence list scroll position
        def sync_alignment_with_sequence(value):
            """Sync alignment position when sequence list scrollbar changes"""
            if self.is_synchronizing:
                return
                
            try:
                self.is_synchronizing = True
                
                # Find first visible row in sequence list
                first_visible_row = find_first_visible_row()
                if first_visible_row is None:
                    self.is_synchronizing = False
                    return
                
                # Calculate line height in alignment view
                alignment_line_height = self.alignment_viewer.fontMetrics().height()
                
                # Calculate new position for alignment based on first visible sequence
                new_pos = first_visible_row * alignment_line_height
                
                # Set the alignment position directly using scrollContents since we're hiding the scrollbar
                self.alignment_viewer.verticalScrollBar().setValue(new_pos)
                
                logger.debug(f"Synced alignment to row {first_visible_row}, pos {new_pos}")
            except Exception as e:
                logger.error(f"Error in sync_alignment_with_sequence: {e}")
            finally:
                self.is_synchronizing = False
        
        # Connect the sequence list scrollbar value change signal
        sequence_vscroll.valueChanged.connect(sync_alignment_with_sequence)
        
        # Initial synchronization
        sync_alignment_with_sequence(sequence_vscroll.value())
        
        logger.info("One-way scrollbar synchronization setup complete")
        
    def tray_icon_activated(self, reason):
        """Handle system tray icon activation"""
        # QSystemTrayIcon.ActivationReason.Trigger = 3 (single click)
        # QSystemTrayIcon.ActivationReason.DoubleClick = 2
        if reason == 3 or reason == 2:  # Single or double click
            self.show()
            self.raise_()
            self.activateWindow()
    
    def closeEvent(self, event):
        """Clean up resources before closing the main window"""
        logger = logging.getLogger("treecraft")
        logger.info("Application closing, cleaning up resources")
        
        # Clean up debug console if it exists
        if self.debug_console:
            # Close the debug console properly to ensure its closeEvent is called
            self.debug_console.close()
            self.debug_console = None
            MainWindow._debug_console_instance = None
        
        # Accept the close event
        event.accept()
        
    def on_export_sequences(self):
        """Export the current sequences to a file"""
        if not self.sequence_list.sequences or len(self.sequence_list.sequences) == 0:
            QMessageBox.warning(self, "Warning", "No sequences to export")
            return
            
        # Create SeqRecord objects for the current sequences with original sequences
        # Make sure to use the actual sequences without padding/alignment dashes
        records = []
        for seq_id, sequence in self.sequence_list.sequences.items():
            description = self.sequence_list.descriptions.get(seq_id, seq_id)
            
            # Remove any trailing padding dashes that might have been added
            clean_sequence = sequence.rstrip('-')
            
            record = SeqRecord(Seq(clean_sequence), id=seq_id, description=description)
            records.append(record)
            
        # Define supported formats, with FASTA as default
        formats = "FASTA (*.fasta *.fa);;FASTQ (*.fastq *.fq);;GENBANK (*.gb *.genbank);;All Files (*)"
        
        # Open save dialog with FASTA as default
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Export Sequences", "", formats
        )
        
        if not file_path:
            return
            
        # Determine format based on selected filter or file extension
        format_name = "fasta"  # Default format
        
        if "FASTQ" in selected_filter:
            format_name = "fastq"
        elif "GENBANK" in selected_filter:
            format_name = "genbank"
        else:
            # Try to determine format from file extension
            ext = os.path.splitext(file_path)[1].lower()
            if ext in (".fq", ".fastq"):
                format_name = "fastq"
            elif ext in (".gb", ".genbank"):
                format_name = "genbank"
                
        # Add appropriate extension if missing
        extensions = {
            "fasta": [".fasta", ".fa"],
            "fastq": [".fastq", ".fq"],
            "genbank": [".gb", ".genbank"]
        }
        
        has_correct_ext = False
        for ext in extensions.get(format_name, []):
            if file_path.lower().endswith(ext):
                has_correct_ext = True
                break
                
        if not has_correct_ext:
            file_path += extensions[format_name][0]
            
        try:
            # Write the original sequences without padding directly to the file
            with open(file_path, "w") as f:
                if format_name == "fasta":
                    # Write FASTA format manually to ensure no alignment/padding
                    for record in records:
                        f.write(f">{record.id} {record.description}\n")
                        f.write(f"{record.seq}\n")
                    self.statusbar.showMessage(f"Exported {len(records)} sequences to {os.path.basename(file_path)} in FASTA format")
                elif format_name in ["fastq", "genbank"]:
                    # For other formats, use SeqIO but with clean sequences
                    SeqIO.write(records, file_path, format_name)
                    self.statusbar.showMessage(f"Exported {len(records)} sequences to {os.path.basename(file_path)} in {format_name.upper()} format")
                else:
                    QMessageBox.warning(self, "Format Error", f"The selected format ({format_name}) is not supported for unaligned sequences.")
                    return
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export sequences: {str(e)}")
    
    def on_export_alignment(self):
        """Export the current alignment to a file"""
        if not self.current_alignment or len(self.current_alignment) == 0:
            QMessageBox.warning(self, "Warning", "No alignment to export")
            return
        
        # Define supported formats with FASTA as default
        formats = "FASTA (*.fasta *.fa);;CLUSTAL (*.aln);;PHYLIP (*.phy);;NEXUS (*.nex);;STOCKHOLM (*.sth *.sto);;All Files (*)"
        
        # Open save dialog with FASTA as default
        file_path, selected_filter = QFileDialog.getSaveFileName(
            self, "Export Alignment", "", formats
        )
        
        if not file_path:
            return
            
        # Determine format based on selected filter or file extension
        format_name = "fasta"  # Default format
        
[end of treecraft/gui/main_window.py]
