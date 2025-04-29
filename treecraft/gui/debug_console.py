import sys
import logging
import re
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QTextEdit, QHBoxLayout, 
                           QCheckBox, QPushButton, QLabel, QDialog)
from PyQt6.QtCore import Qt, QObject, pyqtSignal, pyqtSlot, QTimer
from PyQt6.QtGui import QColor

class QTextEditLogger(logging.Handler, QObject):
    """Custom logging handler that emits logs to a QTextEdit widget in a thread-safe way"""
    
    # Signal to safely emit log messages to the GUI thread
    log_signal = pyqtSignal(str, str)
    
    def __init__(self):
        # Initialize both parent classes
        super().__init__()
        QObject.__init__(self)
        
        # This handler no longer has a direct reference to any widget
        # All access is done via signals to ensure thread safety
        
        # Format string for log messages
        self.format = logging.Formatter('%(asctime)s - [%(levelname)s] - %(name)s - %(message)s')
        
        # Set to prevent crash during shutdown
        self.flushOnClose = False
    
    def emit(self, record):
        try:
            # Convert logging level to string for consistent matching
            level = record.levelname  
            formatted_message = self.format.format(record)
            
            # Use log_signal to safely emit from any thread to GUI thread
            # This will be caught by the DebugConsole which owns this handler
            self.log_signal.emit(formatted_message, level)
        except Exception:
            self.handleError(record)


class DebugConsole(QDialog):
    """Debug console window for displaying log messages
    
    IMPORTANT: This class is designed to be thread-safe. It receives log messages from
    any thread via signals, and only interacts with GUI components in the main thread.
    This prevents crashes when logging from background threads like the alignment thread.
    
    The window properly cleans up its logging handler when closed to prevent "RuntimeError: 
    wrapped C/C++ object of type QTextEditLogger has been deleted" exceptions at exit.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Debug Console")
        self.resize(800, 600)
        
        # Track if handler is installed to safely remove on close
        self.handler_installed = False
        
        # Create main layout
        layout = QVBoxLayout(self)
        
        # Create text area for logs
        self.log_display = QTextEdit()
        self.log_display.setMinimumWidth(700)
        self.log_display.setReadOnly(True)
        self.log_display.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        layout.addWidget(self.log_display)
        
        # Set a nice monospaced font for log display
        font = self.log_display.font()
        font.setFamily("Courier New")
        self.log_display.setFont(font)
        
        # Store all log messages for filtering
        self.all_messages = []
        
        # Create log handler (no direct widget reference)
        self.log_handler = QTextEditLogger()
        self.log_handler.setLevel(logging.DEBUG)
        
        # Connect the log handler's signal to our slot to display messages
        self.log_handler.log_signal.connect(self.process_log_message)
        
        # Add handler to root logger
        root_logger = logging.getLogger()
        root_logger.addHandler(self.log_handler)
        self.handler_installed = True
        
        # Create checkbox layout
        checkbox_layout = QHBoxLayout()
        
        # Create level checkboxes
        self.debug_checkbox = QCheckBox("DEBUG")
        self.debug_checkbox.setChecked(True)
        self.info_checkbox = QCheckBox("INFO")
        self.info_checkbox.setChecked(True)
        self.warning_checkbox = QCheckBox("WARNING")
        self.warning_checkbox.setChecked(True)
        self.error_checkbox = QCheckBox("ERROR")
        self.error_checkbox.setChecked(True)
        
        # Add checkboxes to layout
        checkbox_layout.addWidget(QLabel("Show log levels:"))
        checkbox_layout.addWidget(self.debug_checkbox)
        checkbox_layout.addWidget(self.info_checkbox)
        checkbox_layout.addWidget(self.warning_checkbox)
        checkbox_layout.addWidget(self.error_checkbox)
        checkbox_layout.addStretch()
        
        # Add checkbox layout to main layout
        layout.addLayout(checkbox_layout)
        
        # Add clear button
        button_layout = QHBoxLayout()
        self.clear_button = QPushButton("Clear Console")
        self.clear_button.clicked.connect(self.clear_console)
        button_layout.addStretch()
        button_layout.addWidget(self.clear_button)
        layout.addLayout(button_layout)
        
        # Connect checkbox signals
        self.debug_checkbox.toggled.connect(self.update_visible_levels)
        self.info_checkbox.toggled.connect(self.update_visible_levels)
        self.warning_checkbox.toggled.connect(self.update_visible_levels)
        self.error_checkbox.toggled.connect(self.update_visible_levels)
        
        # Indicate that we're now receiving log messages
        logging.info("Debug console initialized - now capturing log messages")
    
    def process_log_message(self, message, level):
        """Process a log message from the signal (thread-safe)"""
        # Add to our stored messages
        self.all_messages.append((message, level))
        
        # Only display if the level is checked
        if ((level == "DEBUG" and self.debug_checkbox.isChecked()) or
            (level == "INFO" and self.info_checkbox.isChecked()) or 
            (level == "WARNING" and self.warning_checkbox.isChecked()) or
            ((level == "ERROR" or level == "CRITICAL") and self.error_checkbox.isChecked())):
            
            # Set color based on level
            if level == "DEBUG":
                self.log_display.setTextColor(QColor(128, 128, 128))
            elif level == "INFO":
                self.log_display.setTextColor(QColor(0, 0, 255))
            elif level == "WARNING":
                self.log_display.setTextColor(QColor(255, 165, 0))
            elif level == "ERROR" or level == "CRITICAL":
                self.log_display.setTextColor(QColor(255, 0, 0))
            else:
                self.log_display.setTextColor(QColor(0, 0, 0))
            
            # Append message to the display
            self.log_display.append(message)
            
            # Scroll to bottom
            scrollbar = self.log_display.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())
    
    def clear_console(self):
        """Clear the console display"""
        self.log_display.clear()
        self.all_messages = []
    
    def update_visible_levels(self, *args):
        """Filter log messages based on which checkboxes are checked"""
        # Clear display
        self.log_display.clear()
        
        # Get checked levels
        visible_levels = []
        if self.debug_checkbox.isChecked():
            visible_levels.append("DEBUG")
        if self.info_checkbox.isChecked():
            visible_levels.append("INFO")
        if self.warning_checkbox.isChecked():
            visible_levels.append("WARNING")
        if self.error_checkbox.isChecked():
            visible_levels.append("ERROR")
            visible_levels.append("CRITICAL")
        
        # Re-add only the lines with the checked levels from our stored messages
        for message, level in self.all_messages:
            if level in visible_levels:
                # Set color based on level - same as in process_log_message
                if level == "DEBUG":
                    self.log_display.setTextColor(QColor(128, 128, 128))
                elif level == "INFO":
                    self.log_display.setTextColor(QColor(0, 0, 255))
                elif level == "WARNING":
                    self.log_display.setTextColor(QColor(255, 165, 0))
                elif level == "ERROR" or level == "CRITICAL":
                    self.log_display.setTextColor(QColor(255, 0, 0))
                else:
                    self.log_display.setTextColor(QColor(0, 0, 0))
                
                # Add the message
                self.log_display.append(message)
        
        # Scroll to bottom
        scrollbar = self.log_display.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
        
    def closeEvent(self, event):
        """Clean up logging handler when window is closed"""
        if self.handler_installed:
            try:
                # Get the root logger and remove our handler
                root_logger = logging.getLogger()
                root_logger.removeHandler(self.log_handler)
                self.handler_installed = False
                
                # Explicitly disconnect signals
                self.log_handler.log_signal.disconnect()
                
                # Break circular references
                self.log_handler = None
                
                logging.debug("Debug console logging handler removed")
            except Exception as e:
                logging.error(f"Error removing debug console handler: {e}")
                
        # Accept the close event
        event.accept()