#!/usr/bin/env python3
import sys
import logging
from PyQt6.QtWidgets import QApplication
from treecraft.gui.main_window import MainWindow

def main():
    # Get the debug mode flag from the parent app
    logger = logging.getLogger("treecraft")
    
    # Register a cleanup function to ensure proper shutdown
    import atexit
    def cleanup_resources():
        """Clean up any resources that might cause issues during shutdown"""
        logger.debug("Running atexit cleanup handler")
        # Remove all handlers from the root logger to prevent Qt object deletion errors
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
            # Check if it's our custom handler by checking for flushOnClose attribute
            if hasattr(handler, 'flushOnClose'):
                # Skip our custom handler that has flushOnClose=False
                continue
            try:
                root_logger.removeHandler(handler)
            except:
                pass
    
    # Register the cleanup function
    atexit.register(cleanup_resources)
    
    # Add QApplication arguments from command line
    app = QApplication(sys.argv)
    
    # Disable debug messages from Qt unless in debug mode
    if not any(arg == '--debug' for arg in sys.argv):
        # Suppress Qt debug messages
        from PyQt6.QtCore import QtMsgType, qInstallMessageHandler
        def qt_message_handler(msg_type, context, msg):
            # Compare the enum values correctly
            if msg_type in [QtMsgType.QtWarningMsg, QtMsgType.QtCriticalMsg, QtMsgType.QtFatalMsg]:
                print(f"Qt: {msg}", file=sys.stderr)
                
        qInstallMessageHandler(qt_message_handler)
    
    # Create and show the main window
    window = MainWindow()
    window.show()
    
    # Run the application
    return app.exec()

if __name__ == "__main__":
    main()
