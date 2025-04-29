#!/usr/bin/python3

# Treecraft version 1.0 by Alan Rockefeller - April 21, 2025
# https://github.com/AlanRockefeller/Treecraft

import sys
import os
import logging
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='TreeCraft - Phylogenetic Tree Visualization Tool')
parser.add_argument('--debug', action='store_true', help='Enable debug output')
parser.add_argument('--version', action='store_true', help='Show version information')
parser.add_argument('--show-commandlines', '--show-commandline', '--show', action='store_true', 
                    help='Show all command lines for external tools')
args, unknown_args = parser.parse_known_args()

# Configure logging based on debug flag
if args.debug:
    logging.basicConfig(
        level=logging.DEBUG,  # Capture all logs
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
else:
    # Set up more sophisticated logging in normal mode
    # Root logger captures everything for debug console
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    
    # But console handler only shows warnings and errors
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    
    # Clear existing handlers (basicConfig might have added some)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        
    # Add our custom handler
    root_logger.addHandler(console_handler)

logger = logging.getLogger("treecraft")

# Display version information if requested
if args.version:
    print("TreeCraft v1.0")
    sys.exit(0)

# Add the parent directory to sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Create a global variable to indicate if we should show command lines
show_commandlines = args.show_commandlines

if args.debug:
    # Print debug information only in debug mode
    logger.debug("Starting TreeCraft application in debug mode")
    logger.debug(f"Python executable: {sys.executable}")
    logger.debug(f"Python version: {sys.version}")
    logger.debug(f"Python path: {sys.path}")
    logger.debug(f"Show command lines: {show_commandlines}")
else:
    logger.info("Starting TreeCraft application")
    if show_commandlines:
        print("TreeCraft: Command lines for external tools will be displayed")

# Set the flag in the environment so modules can access it
os.environ["TREECRAFT_SHOW_COMMANDLINES"] = "1" if show_commandlines else "0"

try:
    # Try importing PyQt6 directly
    import PyQt6
    if args.debug:
        logger.debug(f"Successfully imported PyQt6 from {PyQt6.__file__}")
except ImportError as e:
    logger.error(f"Failed to import PyQt6: {e}")
    print(f"Error: TreeCraft requires PyQt6, but it could not be imported. Please install PyQt6 and try again.")
    sys.exit(1)

# Import and run the main function
try:
    from treecraft.__main__ import main
    if args.debug:
        logger.debug("Successfully imported main function")
    
    # Pass remaining arguments to the application
    sys.argv[1:] = unknown_args
    main()
except ImportError as e:
    logger.error(f"Failed to import main: {e}")
    print(f"Error: Failed to start TreeCraft: {e}")
    sys.exit(1)
