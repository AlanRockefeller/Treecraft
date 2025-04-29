from PyQt6.QtWidgets import QWidget, QScrollArea, QPushButton, QHBoxLayout, QVBoxLayout, QFrame, QMenu, QInputDialog, QFontDialog
from PyQt6.QtGui import QPainter, QPen, QColor, QFont, QIcon, QTransform, QFontMetrics, QCursor
from PyQt6.QtCore import Qt, QRect, QPoint, QSize, pyqtSignal, QTimer

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os
import subprocess
import tempfile
import shutil
import time
import logging

logger = logging.getLogger("treecraft.phylo_canvas")

class TreeCanvas(QWidget):
    """Widget for drawing the phylogenetic tree"""
    
    # Add signals for node interactions
    node_double_clicked = pyqtSignal(str)
    node_right_clicked = pyqtSignal(str, QPoint)
    rename_signal = pyqtSignal(str)  # Signal for renaming node
    delete_signal = pyqtSignal(str)  # Signal for deleting node
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.tree = None
        self.highlighted_branches = {}  # Dictionary of highlighted branches
        self.clade_labels = {}  # Dictionary of clade labels
        self.setMinimumSize(400, 300)
        self.node_positions = {}  # Store node positions for interaction
        self._dark_mode = True  # Initialize private attribute for dark mode
        self.delete_mode = False  # Track delete mode
        self.reroot_mode = False  # Track reroot mode
        self.hovered_sequence = None  # Track which sequence is being hovered in delete mode
        self.branch_width = 1  # Default branch width
        self.setAutoFillBackground(True)

        # Import QFont here to ensure it's defined before use
        from PyQt6.QtGui import QFont
        # Set default font for sequence names - create a proper font object
        default_font = QFont("Arial", 9)
        self.sequence_font = default_font
        # Log for debugging
        logger.debug(f"TreeCanvas initialized with font: {default_font.family()}, {default_font.pointSize()}pt")
        
        self.show_bootstrap = True
        
        # For clade dragging
        self.drag_clade = None
        self.drag_start = None
        self.dragging = False
        self.setMouseTracking(True)

        # For better interaction feedback
        self.hover_node = None
        # Use OpenHand cursor by default to indicate panning is available
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        self.setMouseTracking(True)  # Enable mouse tracking for hover effects
        self.setToolTip("Click and drag nodes to manipulate the tree.\n"
                        "Click and drag the background to pan the view.\n"
                        "Double-click a leaf node to rename.\n"
                        "Right-click a leaf node for more options.")
                        
        # For rerooting
        self.reroot_mode = False  # Flag to indicate whether we're in reroot mode
        
        # For zooming
        self.scale_factor = 1.0
        self.min_scale = 0.5  # Increase minimum scale to prevent tree reversal
        self.max_scale = 5.0
        
        # For branch spacing
        self.vertical_spacing_factor = 0.9  # Start slightly condensed for vertical spacing
        self.horizontal_spacing_factor = 1.0
        self.min_spacing = 0.15  # Allow much more compression of tree
        # Set different max values for horizontal and vertical to prevent runaway expansion
        self.max_horizontal_spacing = 3.0  # Allow more expansion as well
        self.max_vertical_spacing = 1.25  # Extremely conservative limit for vertical spacing
        
        # For panning
        self.panning = False
        self.pan_start = None
        self.drag_mode = None  # Can be "pan" or "manipulate"
        
        # For node context menu
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)
        
        # Enable mouse tracking for better hover effects and interaction
        self.setMouseTracking(True)
        
        # Signals already defined at class level

    def mousePressEvent(self, event):
        """Handle mouse press events for clade manipulation, panning, rerooting, or deletion"""
        if not self.tree:
            return
            
        pos = event.position().toPoint()
        
        # If we're in delete mode, don't do anything on mouse press - we'll handle deletion on mouse release
        if hasattr(self, 'delete_mode') and self.delete_mode:
            event.accept()
            return
        
        # Find clade at this position
        clade = self.find_clade_at_position(pos)
        
        # Special handling for reroot mode - works for both left and right clicks
        if self.reroot_mode and clade:
            # If in reroot mode, reroot the tree at this clade
            logger.debug(f"Rerooting tree at clade: {clade.name if hasattr(clade, 'name') and clade.name else 'unnamed'}")
            self.reroot_tree_at_clade(clade)
            # Prevent further handling of this click
            event.accept()
            return
            
        from PyQt6.QtCore import Qt
        if event.button() == Qt.MouseButton.LeftButton:
            # If we found a clade and not in reroot mode, start dragging it
            if clade and not self.reroot_mode:
                # Log the clade we're starting to drag
                logger.debug(f"Starting drag on clade: {clade.name if hasattr(clade, 'name') and clade.name else 'unnamed'}")
                
                # Store clade for dragging
                self.drag_clade = clade
                self.drag_start = pos
                self.dragging = True
                self.drag_mode = "manipulate"
                logger.debug(f"Drag mode set to: {self.drag_mode}")
                
                # Set cursor to indicate dragging is active
                self.setCursor(Qt.CursorShape.ClosedHandCursor)
                
                # Store clade name or generate a consistent ID for highlighting
                highlight_key = None
                if hasattr(clade, 'name') and clade.name:
                    highlight_key = clade.name
                else:
                    # For unnamed clades, use internal ID consistently
                    highlight_key = f"internal_{id(clade)}"
                
                # Highlight the clade being dragged
                self.highlighted_branches[highlight_key] = QColor(0, 255, 0)
                
                # Change cursor to indicate dragging
                if hasattr(clade, 'clades') and clade.clades:
                    # Internal node with children - can be manipulated
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
                else:
                    # Leaf node - different cursor
                    self.setCursor(Qt.CursorShape.DragMoveCursor)
                
                # Mark event as handled
                event.accept()
                self.update()
                return
            else:
                # Start panning (if not in reroot mode)
                if not self.reroot_mode:
                    self.pan_start = pos
                    self.panning = True
                    self.drag_mode = "pan"
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
            
        elif event.button() == Qt.MouseButton.RightButton:
            # Check if we're right-clicking on a text label
            for name, (x, y) in self.node_positions.items():
                # Consider all nodes for right-click menu
                # Check if click is in the label area
                if pos.x() > x + 5 and pos.x() <= x + 800 and abs(pos.y() - y) < 15:
                    # Find the actual node name
                    node_name = None
                    
                    # Try to find the terminal node that corresponds to this label
                    for terminal in self.tree.get_terminals():
                        terminal_name = getattr(terminal, 'name', None)
                        # Check if this is the terminal we're looking for
                        if terminal_name and (
                            name == f"leaf_{terminal_name}" or
                            name == f"unnamed_leaf_{terminal_name}" or
                            terminal_name == name
                        ):
                            node_name = terminal_name
                            break
                    
                    # If we found a node name, emit the right-click signal
                    if node_name:
                        logger.debug(f"Right-clicked on node: {node_name}")
                        self.node_right_clicked.emit(node_name, pos)
                        # Accept the event to prevent further handling
                        event.accept()
                        return
            
            # If event wasn't handled above and we have a clade and not in reroot mode
            if clade and not self.reroot_mode:
                # Log detailed info about the clade being dragged
                logger.debug(f"Starting drag on clade: {clade}")
                logger.debug(f"Clade name: {clade.name if clade.name else 'unnamed'}")
                logger.debug(f"Clade has children: {hasattr(clade, 'clades') and bool(clade.clades)}")
                if hasattr(clade, 'clades'):
                    logger.debug(f"Number of child clades: {len(clade.clades)}")
                
                # Start clade manipulation
                self.drag_clade = clade
                self.drag_start = pos
                self.dragging = True
                self.drag_mode = "manipulate"
                
                # Store clade name or generate a consistent ID for highlighting
                highlight_key = None
                if clade.name:
                    highlight_key = clade.name
                else:
                    # For unnamed clades, use internal ID consistently
                    highlight_key = f"internal_{id(clade)}"
                
                # Highlight the clade being dragged - use highlight_key for consistency
                self.highlighted_branches[highlight_key] = QColor(0, 255, 0)
                
                # Change cursor to indicate dragging
                if hasattr(clade, 'clades') and clade.clades:
                    # Internal node with children - can be manipulated
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
                else:
                    # Leaf node - different cursor
                    self.setCursor(Qt.CursorShape.DragMoveCursor)
            elif not self.reroot_mode:
                # Start panning only when we click on the background
                self.pan_start = pos
                self.panning = True
                self.drag_mode = "pan"
                self.setCursor(Qt.CursorShape.ClosedHandCursor)  # Use closed hand to indicate active panning
                
                # Mark the event as accepted to prevent propagation
                event.accept()
            
            self.update()
            
    def mouseMoveEvent(self, event):
        """Handle mouse move events for clade manipulation, panning or hover effects"""
        pos = event.position().toPoint()
        
        if self.dragging and self.drag_clade and self.drag_start and self.drag_mode == "manipulate":
            # Clade manipulation
            dx = pos.x() - self.drag_start.x()
            dy = pos.y() - self.drag_start.y()
            
            if abs(dx) > abs(dy):
                self.setCursor(Qt.CursorShape.SizeHorCursor)  # Horizontal drag
            else:
                self.setCursor(Qt.CursorShape.SizeVerCursor)  # Vertical drag
            
            self.update()
        elif self.panning and self.pan_start and self.drag_mode == "pan":
            # Calculate how much we've moved
            dx = pos.x() - self.pan_start.x()
            dy = pos.y() - self.pan_start.y()
            
            # Get the parent scroll area directly from the hierarchy
            parent = self
            while parent and not isinstance(parent, QScrollArea):
                parent = parent.parent()
                
            if parent and isinstance(parent, QScrollArea):
                # Now we have the scroll area, adjust scrollbars
                hbar = parent.horizontalScrollBar()
                vbar = parent.verticalScrollBar()
                
                # Move scrollbars in the direction of drag
                hbar.setValue(hbar.value() - dx)
                vbar.setValue(vbar.value() - dy)
            else:
                logger.debug("No QScrollArea parent found in hierarchy")
                
            # Update pan start point for continuous motion regardless
            self.pan_start = pos
            
            # Update cursor to show panning is active
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            
            # Prevent event propagation
            event.accept()
        else:
            # Keep track of old hover state for update optimization
            old_hover = self.hover_node
            old_hovered_sequence = self.hovered_sequence
            
            # Now self.over_text_label and self.text_label_node are populated by find_clade_at_position
            # So we can use these values directly
            if hasattr(self, 'over_text_label') and self.over_text_label:
                # We're over text - don't trigger node hover effect
                self.hover_node = None
                
                # Store the hovered sequence name for highlighting in delete mode
                if self.delete_mode and hasattr(self, 'text_label_node') and self.text_label_node:
                    self.hovered_sequence = self.text_label_node
                    # Don't change cursor in delete mode - keep the trash can cursor
                else:
                    self.hovered_sequence = None
                    self.setCursor(Qt.CursorShape.IBeamCursor)  # Text cursor for labels
                
                # Show tooltip for the label
                if hasattr(self, 'node_name_map') and self.text_label_node in self.node_name_map:
                    display_name = self.node_name_map[self.text_label_node]
                else:
                    # Extract name from the text_label_node format
                    if self.text_label_node.startswith("leaf_"):
                        parts = self.text_label_node.split("_", 2)
                        if len(parts) >= 2:
                            display_name = parts[1]
                        else:
                            display_name = self.text_label_node
                    else:
                        display_name = self.text_label_node
                
                if self.delete_mode:
                    self.setToolTip(f"Click to delete: {display_name}")
                else:
                    self.setToolTip(f"Double-click to rename: {display_name}")
            else:
                # Find if we're hovering over a node (but not a label)
                self.hover_node = self.find_clade_at_position(pos)
            
                # Always update cursor on mouse move for better feedback
                if self.reroot_mode:
                    # In reroot mode, we use a crosshair cursor for all nodes
                    if self.hover_node:
                        # Highlight potential reroot targets
                        self.setCursor(Qt.CursorShape.CrossCursor)
                        # Show info about the branch that would become the new root
                        node_name = self.hover_node.name if self.hover_node.name else "this branch"
                        self.setToolTip(f"Click to reroot at {node_name}")
                    else:
                        # Not over a node but still in reroot mode
                        self.setCursor(Qt.CursorShape.CrossCursor)
                        self.setToolTip("Click on a branch to make it the new root")
                elif self.hover_node:
                    # Normal mode interaction with nodes
                    if hasattr(self.hover_node, 'clades') and self.hover_node.clades:
                        # Internal node - can be dragged
                        self.setCursor(Qt.CursorShape.PointingHandCursor)
                        self.setToolTip("Drag to rearrange tree")
                    else:
                        # Leaf node - can be renamed or deleted
                        self.setCursor(Qt.CursorShape.PointingHandCursor)  # Hand cursor for leaf nodes
                        node_name = self.hover_node.name if self.hover_node.name else "unnamed"
                        self.setToolTip(f"Drag to rearrange: {node_name}")
                else:
                    # When hovering over background and not panning,
                    # use the open hand cursor to indicate panning is possible
                    self.setCursor(Qt.CursorShape.OpenHandCursor)
                    self.setToolTip("")  # No tooltip when over empty areas
                
            # Only trigger repaint if hover state changed to avoid unnecessary redraws
            if self.hover_node != old_hover:
                self.update()
    
    def mouseReleaseEvent(self, event):
        """Handle mouse release events for clade manipulation, panning, or node deletion"""
        import logging
        logger = logging.getLogger("treecraft")
        
        logger.debug(f"MouseReleaseEvent triggered with button: {event.button()}")
        
        from PyQt6.QtCore import Qt
        if event.button() == Qt.MouseButton.LeftButton:
            logger.debug(f"LEFT button released, drag_mode: {getattr(self, 'drag_mode', None)}, dragging: {getattr(self, 'dragging', False)}")
            # Check if we're in delete mode
            if hasattr(self, 'delete_mode') and self.delete_mode:
                pos = event.position().toPoint()
                logger.debug(f"MouseReleaseEvent in delete_mode at position {pos.x()}, {pos.y()}")
                
                # For delete mode, first try to identify a node directly under the cursor
                # This is a more direct approach than searching through all nodes
                
                # Try to find the text label first - it's easier to click on
                for name, (x, y) in self.node_positions.items():
                    # For leaf nodes, check if we clicked on the text area
                    if name.startswith("leaf_"):
                        # Check if cursor is near this label position for deletion
                        label_x = int(x + 5)  # Offset from node position
                        label_y = int(y) 
                        
                        # Use simple distance calculation instead of QRect for reliability
                        # Check if cursor is within 800 pixels horizontally and 10 pixels vertically
                        if (pos.x() >= label_x and pos.x() <= label_x + 800 and 
                            abs(pos.y() - label_y) <= 10):
                            logger.debug(f"Direct hit on text label: {name}")
                            # Extract the terminal name from the leaf_ prefix format - handle position suffix
                            # Format is typically "leaf_TERMINAL_NAME_POSITION" 
                            if "_" in name:
                                parts = name.split("_")
                                # Special case for RAxML trees with location data suffixed (e.g. "_Washington_US_35")
                                if parts[0] == "leaf" and len(parts) > 3:
                                    # For RAxML trees location data may be appended with underscore
                                    # First try using our node_name_map which should have the exact original name
                                    if hasattr(self, 'node_name_map') and name in self.node_name_map:
                                        terminal_name = self.node_name_map[name]
                                        logger.debug(f"Using mapped name from node_name_map: {terminal_name}")
                                    else:
                                        # Otherwise, search for common location identifiers
                                        location_parts = []
                                        for i, part in enumerate(parts):
                                            if i > 0:  # Skip the "leaf_" part
                                                # Check if this part looks like a location identifier
                                                if part in ["US", "EU", "CA", "UK", "China", "IN", "JP", "AU"] or part in ["Washington", "California", "Oregon", "BC", "Indiana", "Florida", "Texas", "NY", "Alaska", "Hawaii"] or part.isdigit():
                                                    location_parts.append(i)
                                        
                                        if location_parts:
                                            # We found location parts, use everything before the first location
                                            first_loc = min(location_parts)
                                            terminal_name = "_".join(parts[1:first_loc])
                                            logger.debug(f"Extracted terminal name before location: '{terminal_name}'")
                                        else:
                                            # Regular leaf node - take all parts between first and last
                                            terminal_name = "_".join(parts[1:-1]) if len(parts) > 2 else parts[1]
                                            logger.debug(f"Extracted terminal name using standard format: '{terminal_name}'")
                                else:
                                    if len(parts) >= 2:
                                        # Normal case with position suffix - extract just the name part
                                        # Join all parts between the first and last with underscores
                                        terminal_name = "_".join(parts[1:-1]) if len(parts) > 2 else parts[1]
                                        logger.debug(f"Extracted terminal name '{terminal_name}' from '{name}'")
                                    else:
                                        terminal_name = name
                            else:
                                terminal_name = name.replace("leaf_", "", 1)
                            
                            # Check the node name map for more accurate name resolution
                            if hasattr(self, 'node_name_map') and name in self.node_name_map:
                                mapped_name = self.node_name_map[name]
                                logger.debug(f"Found mapped name '{mapped_name}' for node '{name}'")
                                terminal_name = mapped_name
                            
                            # Check if we have a valid tree before trying to access terminals
                            exists = False
                            actual_terminal_name = None  # Will store the actual name found in the tree
                            
                            if self.tree and hasattr(self.tree, 'get_terminals'):
                                try:
                                    # Log more detailed information
                                    logger.debug(f"Searching tree terminals for name '{terminal_name}'")
                                    
                                    # First, try exact matching
                                    for terminal in self.tree.get_terminals():
                                        term_name = getattr(terminal, 'name', '')
                                        
                                        # Check for exact match first
                                        if term_name == terminal_name:
                                            exists = True
                                            actual_terminal_name = term_name
                                            logger.debug(f"Found exact terminal match: '{terminal_name}'")
                                            break
                                    
                                    # If not found, use more flexible matching
                                    if not exists:
                                        for terminal in self.tree.get_terminals():
                                            term_name = getattr(terminal, 'name', '')
                                            
                                            # Special case for accession numbers with version (e.g., KX897432.1)
                                            if '.' in terminal_name and '.' in term_name:
                                                term_base = term_name.split('.')[0]
                                                terminal_base = terminal_name.split('.')[0]
                                                
                                                if term_base == terminal_base:
                                                    exists = True
                                                    actual_terminal_name = term_name
                                                    logger.debug(f"Matched on accession base: '{terminal_base}' → '{term_name}'")
                                                    break
                                            
                                            # Try substring matching as a fallback
                                            elif (terminal_name in term_name) or (term_name in terminal_name):
                                                exists = True
                                                actual_terminal_name = term_name
                                                logger.debug(f"Found substring match: '{terminal_name}' ↔ '{term_name}'")
                                                break
                                        
                                except Exception as e:
                                    logger.error(f"Error accessing tree terminals: {e}")
                                    exists = False
                            else:
                                logger.debug("Tree is not available - cannot verify terminals")
                            
                            # Always proceed with deletion in delete mode, even if verification fails
                            if exists:
                                logger.debug(f"Verified terminal '{actual_terminal_name}' exists, deleting directly")
                                # Use the actual terminal name we found in the tree
                                final_name = actual_terminal_name
                            else:
                                # In delete mode, proceed with deletion even without exact match
                                logger.warning(f"No exact match found for node '{terminal_name}' but continuing with deletion anyway")
                                final_name = terminal_name
                            
                            # Emit the delete signal directly
                            self.delete_signal.emit(final_name)
                            # Also try direct deletion through the main window
                            QTimer.singleShot(10, lambda name=final_name: self.direct_delete_node(name))
                            # Force an immediate update of the view
                            self.update()
                            return
                
                # Only attempt node finding if we have a valid tree
                if not self.tree or not hasattr(self.tree, 'get_terminals'):
                    logger.debug("No tree available for delete operation, ignoring click")
                    return
                
                # Find the closest node to this position as fallback
                closest_node = None
                closest_distance = float('inf')
                
                logger.debug(f"Looking for nodes near position {pos.x()}, {pos.y()}")
                leaf_nodes_count = 0
                
                # Verify we have positions to check
                if not hasattr(self, 'node_positions') or not self.node_positions:
                    logger.error("No node_positions available for deletion")
                    return
                    
                for name, (x, y) in self.node_positions.items():
                    # Only consider leaf nodes for deletion
                    if not (name.startswith("leaf_") or name.startswith("unnamed_leaf_")):
                        continue
                        
                    leaf_nodes_count += 1
                    
                    # Calculate distance
                    distance = ((pos.x() - x) ** 2 + (pos.y() - y) ** 2) ** 0.5
                    logger.debug(f"Leaf node {name} at ({x}, {y}), distance: {distance:.1f}")
                    
                    # In delete mode use a much larger click area for better usability (100 pixels)
                    click_radius = 100 if self.delete_mode else 40
                    if distance < click_radius and distance < closest_distance:
                        closest_node = name
                        closest_distance = distance
                        logger.debug(f"Found closest node so far: {name} at distance {distance:.1f}")
                
                # Add debug info
                logger.debug(f"Found {leaf_nodes_count} leaf nodes total")
                
                # Also check if we're clicking directly on a text label
                if not closest_node and pos.x() > 0 and hasattr(self, 'over_text_label') and hasattr(self, 'text_label_node'):
                    logger.debug(f"Text label attributes: over_text_label={getattr(self, 'over_text_label', None)}, text_label_node={getattr(self, 'text_label_node', None)}")
                    if getattr(self, 'over_text_label', False) and getattr(self, 'text_label_node', None):
                        logger.debug(f"Clicking directly on text label: {self.text_label_node}")
                        if self.text_label_node.startswith("leaf_"):
                            closest_node = self.text_label_node
                            logger.debug(f"Using text label as target: {closest_node}")
                
                # If we found a node to delete
                if closest_node:
                    logger.debug(f"Selected node for deletion: {closest_node}")
                    # Find the corresponding clade in the tree - with safe checks
                    try:
                        terminal_found = False
                        terminal_name_to_delete = None
                        
                        for terminal in self.tree.get_terminals():
                            terminal_name = getattr(terminal, 'name', None)
                            # If we found the node to delete
                            if terminal_name and (terminal_name == closest_node or 
                                                f"leaf_{terminal_name}" == closest_node or
                                                (closest_node and closest_node.startswith("leaf_") and terminal_name in closest_node)):
                                # Found a match
                                terminal_found = True
                                terminal_name_to_delete = terminal_name
                                break
                        
                        if terminal_found and terminal_name_to_delete:
                            # Emit the delete signal with the node name
                            logger.debug(f"Deleting node: {terminal_name_to_delete}")
                            logger.debug(f"Emitting delete_signal from TreeCanvas for node: {terminal_name_to_delete}")
                            
                            try:
                                # Add a direct call to simulate right-click delete in the main window
                                # This is a fallback in case the signal chain isn't working properly
                                QTimer.singleShot(10, lambda name=terminal_name_to_delete: self.direct_delete_node(name))
                            except Exception as e:
                                logger.error(f"Error setting up QTimer: {e}")
                            
                            # Also emit the normal signal
                            self.delete_signal.emit(terminal_name_to_delete)
                            logger.debug(f"Signal emitted successfully")
                            
                            # Force update to show deletion immediately
                            self.update()
                    except AttributeError as e:
                        logger.error(f"Error accessing terminals: {e}")
                        return  # Skip further processing
                
                # Continue with normal cursor behavior
                return
            
            # Normal behavior for drag and manipulate
            if self.dragging and self.drag_clade and self.drag_mode == "manipulate":
                logger.debug(f"Processing drag and manipulate on release")
                logger.debug(f"Dragged clade: {getattr(self.drag_clade, 'name', 'unnamed')}")
                
                # Handle clade manipulation
                end_pos = event.position().toPoint()
                
                # Calculate drag distance and direction
                dx = end_pos.x() - self.drag_start.x()
                dy = end_pos.y() - self.drag_start.y()
                logger.debug(f"Drag distance: dx={dx}, dy={dy}")
                
                # Only manipulate if there was a significant drag
                if abs(dx) > 5 or abs(dy) > 5:
                    # More forgiving direction detection - just use the larger component
                    # Log the drag attempt
                    logger.debug(f"Attempting to manipulate clade via drag (dx={dx}, dy={dy})")
                    logger.debug(f"Drag clade: {self.drag_clade}")
                    logger.debug(f"Drag clade name: {self.drag_clade.name if self.drag_clade.name else 'unnamed'}")
                    
                    # Only perform the operation if this is an internal node with children
                    if hasattr(self.drag_clade, 'clades') and self.drag_clade.clades:
                        logger.debug(f"Clade has {len(self.drag_clade.clades)} children - can be manipulated")
                        
                        if abs(dx) > abs(dy):
                            # Horizontal drag - swap clade order in parent
                            logger.debug("Horizontal drag detected - attempting to swap with siblings")
                            parent_clade = self.find_parent_clade(self.tree.root, self.drag_clade)
                            if parent_clade:
                                logger.debug(f"Found parent clade: {parent_clade}")
                                
                                # Find the dragged clade's index in parent's children
                                try:
                                    index = parent_clade.clades.index(self.drag_clade)
                                    logger.debug(f"Dragged clade at index {index} of {len(parent_clade.clades)} children")
                                    
                                    # Simplify to a direct swap - more reliable
                                    logger.debug(f"Swapping all children in parent")
                                    success = self.swap_clades(parent_clade)
                                    
                                    # Force layout recalculation (redundant but being extra safe)
                                    self.calculate_tree_layout()
                                    # Force immediate redraw
                                    self.update()
                                except ValueError:
                                    # Fallback to traditional swap
                                    logger.debug("Cannot find index of dragged clade, using traditional swap")
                                    success = self.swap_clades(parent_clade)
                                
                                # Update tree display
                                self.update()
                                logger.debug(f"Swap clades result: {success}")
                            else:
                                logger.debug("No parent clade found - cannot swap")
                        else:
                            # Vertical drag - rotate child clades
                            logger.debug("Vertical drag detected - attempting to rotate children")
                            
                            # Simplify to a direct rotation call
                            if self.drag_clade and hasattr(self.drag_clade, 'clades') and len(self.drag_clade.clades) > 1:
                                # Always rotate the same way for consistency
                                logger.debug("Performing standard rotation on clade")
                                success = self.rotate_clade(self.drag_clade)
                                
                                # Force update (redundant but being extra safe)
                                self.update()
                            else:
                                logger.debug("Cannot rotate: clade has no children")
                                success = False
                                
                            logger.debug(f"Rotate clade result: {success}")
                    else:
                        logger.debug("Clade has no children or clades attribute - cannot manipulate")
                        
                # Clear the drag state and highlighting - handle both named and unnamed clades
                if self.drag_clade:
                    # Generate the same consistent key we used when highlighting
                    highlight_key = self.drag_clade.name if self.drag_clade.name else f"internal_{id(self.drag_clade)}"
                    
                    # Remove this specific highlight
                    if highlight_key in self.highlighted_branches:
                        del self.highlighted_branches[highlight_key]
                        logger.debug(f"Removed highlight for {highlight_key}")
                    else:
                        logger.debug(f"No highlight found for {highlight_key}")
                
                # Also check for any other highlights that might be stuck
                self.clear_highlights()  # This ensures all highlights are cleared
                
                # Reset cursor to default or hover state
                if self.hover_node and self.hover_node.clades:
                    self.setCursor(Qt.CursorShape.PointingHandCursor)
                else:
                    self.setCursor(Qt.CursorShape.OpenHandCursor)
            
            elif self.panning and self.drag_mode == "pan":
                # End panning - switch back to open hand to indicate pan is possible
                self.setCursor(Qt.CursorShape.OpenHandCursor)
            elif not self.dragging and not self.panning:
                # Released without dragging, ensure we have the right cursor
                if self.hover_node and self.hover_node.clades:
                    self.setCursor(Qt.CursorShape.PointingHandCursor)
                else:
                    self.setCursor(Qt.CursorShape.OpenHandCursor)
            
            # Reset drag/pan state
            self.dragging = False
            self.panning = False
            self.drag_clade = None
            self.drag_start = None
            self.pan_start = None
            self.drag_mode = None
            
            self.update()
                
    def mouseDoubleClickEvent(self, event):
        """Handle mouse double-click events for node renaming"""
        if not self.tree:
            return
            
        pos = event.position().toPoint()
        
        # Check if we're double-clicking on a text label
        is_over_text_label = False
        label_node_name = None
        renamed_clade = None
        
        # Search for a text label at this position
        for name, (x, y) in self.node_positions.items():
            # Only consider leaf nodes
            if not (name.startswith("leaf_") or name.startswith("unnamed_leaf_")):
                continue
                
            # Check if click is in the label area (text part)
            if pos.x() > x + 5 and pos.x() <= x + 800 and abs(pos.y() - y) < 15:
                is_over_text_label = True
                label_node_name = name
                
                # Extract the actual name to pass to the signal
                node_name = None
                
                # Find the terminal node that corresponds to this label
                for terminal in self.tree.get_terminals():
                    terminal_name = getattr(terminal, 'name', None)
                    # Check if this is the terminal we're looking for
                    if terminal_name and (
                        name == f"leaf_{terminal_name}" or
                        name == f"unnamed_leaf_{terminal_name}" or
                        terminal_name == label_node_name
                    ):
                        node_name = terminal_name
                        break
                
                # If we found a node name, emit the signal
                if node_name:
                    logger.debug(f"Double-clicked on node: {node_name}")
                    self.node_double_clicked.emit(node_name)
                    return
                
                # Find the actual clade for this label
                for terminal in self.tree.get_terminals():
                    # Get the original name from our mapping
                    original_name = None
                    if hasattr(self, 'node_name_map') and name in self.node_name_map:
                        original_name = self.node_name_map[name]
                    else:
                        # Extract from the label_name format
                        parts = name.split("_", 2)
                        if len(parts) >= 2:
                            original_name = parts[1]
                            
                    if terminal.name == original_name:
                        renamed_clade = terminal
                        break
                        
                break
        
        # If we didn't find a label, try the standard node finder
        if not is_over_text_label:
            renamed_clade = self.find_clade_at_position(pos)
        
        # Log position and found clade using debug level
        logger.debug(f"Double-click at position {pos.x()}, {pos.y()}")
        if renamed_clade:
            logger.debug(f"Found clade to rename: {renamed_clade.name}")
        else:
            logger.debug("No clade found at click position")
            
        # Only allow renaming leaf nodes (terminal nodes)
        if renamed_clade and hasattr(renamed_clade, 'clades') and not renamed_clade.clades:
            # This is a leaf node, trigger rename signal if connected
            if self.rename_signal:
                logger.debug(f"Emitting rename signal for: {renamed_clade.name}")
                self.rename_signal.emit(renamed_clade.name)
                # Accept the event to prevent it from propagating
                event.accept()
                return
        
        # Let the event propagate if not handled
        event.ignore()
    
    def show_context_menu(self, position):
        """Show context menu for right-clicked nodes"""
        if not self.tree:
            return
            
        # Convert to QPoint
        pos = position.toPoint() if hasattr(position, 'toPoint') else position
        
        # Log position for debugging
        logger.debug(f"Context menu requested at position {pos.x()}, {pos.y()}")
        
        # Special handling for right-clicking on leaf node text labels
        # since find_clade_at_position no longer returns leaf nodes for hover
        leaf_clade = None
        
        # Check if we're clicking on a leaf node label (similar to the text detection in mouseMoveEvent)
        for name, (x, y) in self.node_positions.items():
            # Only check leaf nodes for context menu
            is_leaf = name.startswith("leaf_") or name.startswith("unnamed_leaf_")
            
            if is_leaf:
                # Check if we're clicking on node or its label
                if pos.x() >= x and pos.x() <= x + 300 and abs(pos.y() - y) < 15:
                    # Get the original clade name from our mapping
                    orig_name = None
                    if hasattr(self, 'node_name_map') and name in self.node_name_map:
                        orig_name = self.node_name_map[name]
                    else:
                        # Extract name from the format
                        if name.startswith("leaf_"):
                            parts = name.split("_", 2)
                            if len(parts) >= 2:
                                orig_name = parts[1]
                    
                    if orig_name:
                        # Find the actual clade in the tree
                        for term in self.tree.get_terminals():
                            if term.name == orig_name:
                                leaf_clade = term
                                logger.debug(f"Found leaf clade for context menu: {orig_name}")
                                break
                    
                    if leaf_clade:
                        break
        
        # If we didn't find a leaf node, try the regular node finder 
        # (which might return internal nodes)
        if not leaf_clade:
            leaf_clade = self.find_clade_at_position(pos)
        
        if leaf_clade:
            logger.debug(f"Found clade for menu: {leaf_clade.name}, has children: {hasattr(leaf_clade, 'clades') and bool(leaf_clade.clades)}")
        else:
            logger.debug("No clade found at context menu position")
            
        # Only show context menu for leaf nodes
        if leaf_clade and (not hasattr(leaf_clade, 'clades') or not leaf_clade.clades):
            logger.debug(f"Showing context menu for {leaf_clade.name}")
            
            # Create context menu
            context_menu = QMenu(self)
            
            # Add rename action
            rename_action = context_menu.addAction("Rename Sequence")
            
            # Add delete action
            delete_action = context_menu.addAction("Delete Sequence")
            
            # Reset dragging state before showing menu
            self.dragging = False
            self.panning = False
            self.drag_mode = None
            self.drag_clade = None
            
            # Show the menu
            action = context_menu.exec(self.mapToGlobal(pos))
            
            # Handle the action, using QTimer to add a slight delay
            from PyQt6.QtCore import QTimer
            if action == rename_action:
                logger.debug(f"Rename action selected for {leaf_clade.name}")
                QTimer.singleShot(10, lambda: self.rename_signal.emit(leaf_clade.name))
            elif action == delete_action:
                logger.debug(f"Delete action selected for {leaf_clade.name}")
                QTimer.singleShot(10, lambda: self.delete_signal.emit(leaf_clade.name))
            return True
        else:
            # Only show generic context menu for non-leaf nodes
            logger.debug("No context menu shown for non-leaf node or empty space")
            return False
    
    def keyPressEvent(self, event):
        """Handle key press events"""
        # Check for Escape key to exit reroot mode
        if event.key() == Qt.Key.Key_Escape and self.reroot_mode:
            logger.debug("Escape key pressed, exiting reroot mode")
            self.exit_reroot_mode()
            # Mark event as handled
            event.accept()
            return
        
        # For debugging purposes - will help implement keyboard shortcuts
        if 'TREECRAFT_DEBUG' in os.environ:
            logger.debug(f"Key press: {event.key()}")
        super().keyPressEvent(event)
    
    def wheelEvent(self, event):
        """Handle mouse wheel events for zooming"""
        # Get the vertical scroll delta directly from both X and Y
        delta_y = event.angleDelta().y()
        delta_x = event.angleDelta().x()
        
        # Skip verbose wheel event logging in normal operation
        
        # Check for Alt modifier
        if event.modifiers() & Qt.KeyboardModifier.AltModifier:
            # Alternative vertical scrolling implementation
            # Use keyboard UP/DOWN simulation
            
            # Assume scrollbar can be found in the main window
            from PyQt6.QtWidgets import QApplication
            main_window = QApplication.activeWindow()
            
            pass  # Skip verbose debug logging
            
            # Try to find the PhyloCanvas instance - it's the actual scrollable area
            from treecraft.gui.phylo_canvas import PhyloCanvas
            
            # Try to get the scroll area directly from this widget's parent chain
            # This should reliably find our parent PhyloCanvas
            parent = self.parent()
            while parent:
                if isinstance(parent, PhyloCanvas):
                    # Found the PhyloCanvas
                    scroll_area = parent
                    pass  # Skip verbose debug logging
                    
                    # Direct scrollbar access
                    vbar = scroll_area.verticalScrollBar()
                    current = vbar.value()
                    
                    # Hard-coded step size for consistent behavior
                    step = 100
                    
                    # Check both X and Y deltas for scrolling direction
                    # Your mouse is sending horizontal scroll events (delta_x)
                    if delta_x > 0:
                        # Horizontal scroll right = vertical scroll up
                        vbar.setValue(current - step)
                    elif delta_x < 0:
                        # Horizontal scroll left = vertical scroll down
                        vbar.setValue(current + step)
                    elif delta_y > 0:
                        # Vertical scroll up
                        vbar.setValue(current - step)
                    else:
                        # Vertical scroll down or default case
                        vbar.setValue(current + step)
                    
                    # Prevent further handling
                    event.accept()
                    return
                parent = parent.parent()
            
            logger.debug("Could not find parent PhyloCanvas for scrolling")
            return
        
        # No Alt modifier - use wheel for zooming
        if delta_y > 0:
            # Zoom in on wheel up
            scale_factor = 1.1  # Zoom in factor
            self.scale_factor = min(self.scale_factor * scale_factor, self.max_scale)
        else:
            # Zoom out on wheel down
            scale_factor = 0.9  # Zoom out factor
            
            # Calculate the new scale factor but don't apply it yet
            new_scale = self.scale_factor * scale_factor
            
            # Check if this would cause the tree to reverse
            # We consider a scale below min_scale as a potential reversal
            if new_scale < self.min_scale:
                logger.debug(f"Preventing scale below minimum: {new_scale:.2f} < {self.min_scale:.2f}")
                new_scale = self.min_scale
            
            # Apply the safe scale factor
            self.scale_factor = new_scale
        
        # Find parent QScrollArea to notify of size changes
        parent = self
        while parent and not isinstance(parent, QScrollArea):
            parent = parent.parent()
            
        if parent and isinstance(parent, QScrollArea):
            # Update the parent's widget size based on new scale
            parent_canvas = parent
            # The parent expects the scale factor to be accessible
            if hasattr(parent_canvas, "updateTreeCanvasSize"):
                parent_canvas.updateTreeCanvasSize()
        
        # Force layout recalculation and redraw
        self.calculate_tree_layout()
        self.update()
        event.accept()

    def set_tree(self, tree):
        """Set the tree to display"""
        logger.debug(f"Setting tree in TreeCanvas: {type(tree)}")
        try:
            # Store old tree for comparison
            old_tree = self.tree
            self.tree = tree
            
            # Initialize maps if they don't exist
            if not hasattr(self, 'node_name_map'):
                self.node_name_map = {}
            if not hasattr(self, 'object_id_map'):
                self.object_id_map = {}
            
            # Handle case where tree is None
            if tree is None:
                logger.warning("Setting None tree in TreeCanvas")
                # Clear existing data
                self.highlighted_branches = {}
                self.clade_labels = {}
                self.node_positions = {}
                self.node_name_map = {}
                self.object_id_map = {}
                self.update()
                return
                
            # Verify that this is a Bio.Phylo.BaseTree.Tree instance
            from Bio.Phylo.BaseTree import Tree
            if not isinstance(tree, Tree):
                logger.error(f"Object is not a Bio.Phylo.BaseTree.Tree: {type(tree)}")
                # Still allow it to be set, but log the warning
            
            # Verify tree structure and log information for debugging
            if hasattr(tree, 'root'):
                logger.debug(f"Tree has root: {tree.root}")
                
                # Diagnose any tree issues
                if not tree.root:
                    logger.error("Tree root exists but is None or empty!")
                elif not hasattr(tree.root, 'clades'):
                    logger.error("Tree root has no clades attribute!")
                elif not tree.root.clades:
                    logger.error("Tree root has empty clades list!")
                    
                # Check terminal nodes and their names
                try:
                    terminals = list(tree.get_terminals())
                    logger.debug(f"Terminal nodes: {len(terminals)}")
                    
                    # Log the first few terminal names to debug
                    for i, term in enumerate(terminals[:5]):  # Show first 5 terminals
                        name = getattr(term, 'name', '[no name]')
                        branch_length = getattr(term, 'branch_length', '[no branch_length]')
                        logger.debug(f"Terminal {i} name: '{name}', branch length: {branch_length}")
                    
                    # Warning for very few terminals
                    if len(terminals) < 2:
                        logger.warning(f"Tree has only {len(terminals)} terminal nodes - may not display properly")
                except Exception as term_err:
                    logger.error(f"Error getting terminals: {term_err}")
            else:
                logger.warning("Tree has no root attribute - may not be a valid tree object")
                
            # Clear previous data to ensure clean rendering
            self.highlighted_branches = {}
            self.clade_labels = {}
            self.node_positions = {}
            
            # Force a complete layout recalculation
            self.calculate_tree_layout()
            
            # Force redraw
            self.update()
            
            # Log whether this is a new tree or the same one
            if old_tree is not tree:
                logger.debug("New tree object set (different from previous)")
            else:
                logger.debug("Same tree object was reset")
                
        except Exception as e:
            import traceback
            logger.error(f"Error setting tree in TreeCanvas: {e}")
            logger.error(traceback.format_exc())
            # Don't clear the tree on error - it might be partially working
            self.update()
    
    def highlight_branch(self, clade_name, color=QColor(255, 0, 0)):
        """Highlight a branch with the given color"""
        self.highlighted_branches[clade_name] = color
        self.update()
    
    def add_clade_label(self, clade_name, label):
        """Add a label to a clade"""
        self.clade_labels[clade_name] = label
        self.update()
    
    def clear_highlights(self):
        """Clear all branch highlights"""
        self.highlighted_branches = {}
        self.update()
        
    def increase_vertical_spacing(self):
        """Increase the vertical spacing between branches"""
        # Only allow scaling up to the maximum value
        if self.vertical_spacing_factor >= self.max_vertical_spacing:
            logger.debug(f"Maximum vertical spacing reached: {self.vertical_spacing_factor:.1f}x")
            self.update_status_message("Maximum vertical spacing reached")
            return False
            
        # Use extremely small increments for vertical spacing to prevent excessive expansion
        if self.vertical_spacing_factor > 1.1:
            increment = 0.01  # Micro increment for large scales
        elif self.vertical_spacing_factor > 1.0:
            increment = 0.02  # Tiny increment for medium scales
        else:
            increment = 0.05  # Small increment for low scales
            
        # Apply increment with bounds checking
        old_value = self.vertical_spacing_factor
        self.vertical_spacing_factor += increment
        self.vertical_spacing_factor = min(self.vertical_spacing_factor, self.max_vertical_spacing)
        
        # Check if we're getting too large and limit growth
        if self.required_space_too_large():
            # Revert to previous value if we're getting too big
            logger.debug("Tree would become too large, limiting vertical growth")
            self.vertical_spacing_factor = old_value
            self.update_status_message("Maximum reasonable tree size reached")
            return False
            
        self.calculate_tree_layout()
        self.update()
        
        # Log spacing change for debugging
        logger.debug(f"Vertical spacing increased to {self.vertical_spacing_factor:.1f}x")
        
        # Update parent status bar if available
        self.update_status_with_spacing_info()
        return True
        
    def required_space_too_large(self):
        """Check if the tree would require too much space with current settings"""
        if not self.tree:
            return False
            
        # Estimate required height based on number of leaf nodes and spacing
        leaf_count = len(list(self.tree.get_terminals()))
        min_spacing = 25  # Increased minimum pixels between leaf nodes
        
        # Adaptive base spacing - more space for larger trees
        if leaf_count > 100:
            base_spacing = 30  # More space for very large trees
        elif leaf_count > 50:
            base_spacing = 25  # More space for large trees
        else:
            base_spacing = max(20, 600 / max(leaf_count, 1))  # Standard spacing
            
        estimated_height = (base_spacing * self.vertical_spacing_factor * leaf_count) + 200  # Added more padding
        
        # Calculate an adaptive limit based on number of sequences
        # Higher limits for larger trees to ensure they can be fully displayed
        if leaf_count > 100:
            max_height = 6000  # Much higher limit for very large trees
        elif leaf_count > 50:
            max_height = 5000  # Higher limit for large trees
        elif leaf_count > 20:
            max_height = 4000  # Moderate limit for medium trees
        else:
            max_height = 3000  # Standard limit for small trees
            
        # Return true if estimated height would be excessive
        return estimated_height > max_height
        
    def decrease_vertical_spacing(self):
        """Decrease the vertical spacing between branches"""
        # Only allow scaling down to the minimum value
        if self.vertical_spacing_factor <= self.min_spacing:
            logger.debug(f"Minimum vertical spacing reached: {self.vertical_spacing_factor:.1f}x")
            self.update_status_message("Minimum vertical spacing reached")
            return False
            
        # Use matching micro decrements for vertical spacing
        if self.vertical_spacing_factor > 1.1:
            increment = 0.01  # Micro decrement for large scales
        elif self.vertical_spacing_factor > 1.0:
            increment = 0.02  # Tiny decrement for medium scales
        else:
            increment = 0.05  # Small decrement for low scales
            
        self.vertical_spacing_factor -= increment
        self.vertical_spacing_factor = max(self.vertical_spacing_factor, self.min_spacing)
        self.calculate_tree_layout()
        self.update()
        
        # Log spacing change for debugging
        logger.debug(f"Vertical spacing decreased to {self.vertical_spacing_factor:.1f}x")
        
        # Update parent status bar if available
        self.update_status_with_spacing_info()
        return True
        
    def increase_horizontal_spacing(self):
        """Increase the horizontal branch lengths"""
        # Only allow scaling up to the maximum value
        if self.horizontal_spacing_factor >= self.max_horizontal_spacing:
            logger.debug(f"Maximum horizontal spacing reached: {self.horizontal_spacing_factor:.1f}x")
            self.update_status_message("Maximum horizontal spacing reached")
            return False
            
        # Increase by smaller increments at larger scales for better control
        if self.horizontal_spacing_factor > 2.0:
            increment = 0.1  # Small increment for fine control at large scales
        elif self.horizontal_spacing_factor < 0.5:
            increment = 0.05  # Very small increment when coming out of high compression
        else:
            increment = 0.2  # Standard increment for normal ranges
            
        self.horizontal_spacing_factor += increment
        self.horizontal_spacing_factor = min(self.horizontal_spacing_factor, self.max_horizontal_spacing)
        self.calculate_tree_layout()
        self.update()
        
        # Log spacing change for debugging
        logger.debug(f"Horizontal spacing increased to {self.horizontal_spacing_factor:.1f}x")
        
        # Update parent status bar if available
        self.update_status_with_spacing_info()
        return True
        
    def decrease_horizontal_spacing(self):
        """Decrease the horizontal branch lengths"""
        # Only allow scaling down to the minimum value
        if self.horizontal_spacing_factor <= self.min_spacing:
            logger.debug(f"Minimum horizontal spacing reached: {self.horizontal_spacing_factor:.1f}x")
            self.update_status_message("Minimum horizontal spacing reached")
            return False
            
        # Decrease by smaller increments at smaller scales, but use more aggressive scaling
        # for better compression control
        if self.horizontal_spacing_factor < 0.5:
            # Very small increment for fine control at high compression
            increment = 0.05
        elif self.horizontal_spacing_factor < 1.0:
            # Small increment for fine control at moderate compression
            increment = 0.1
        else:
            # Larger increment for faster compression from expanded state
            increment = 0.2
            
        self.horizontal_spacing_factor -= increment
        self.horizontal_spacing_factor = max(self.horizontal_spacing_factor, self.min_spacing)
        self.calculate_tree_layout()
        self.update()
        
        # Log spacing change for debugging
        logger.debug(f"Horizontal spacing decreased to {self.horizontal_spacing_factor:.1f}x")
        
        # Update parent status bar if available
        self.update_status_with_spacing_info()
        return True
        
    def update_status_with_spacing_info(self):
        """Update the main window status bar with spacing information if available"""
        try:
            # Find the main window to update status
            from PyQt6.QtWidgets import QApplication
            main_window = QApplication.activeWindow()
            
            if main_window and hasattr(main_window, 'statusbar'):
                # Create status message with current spacing factors
                h_status = f"{self.horizontal_spacing_factor:.1f}x"
                v_status = f"{self.vertical_spacing_factor:.1f}x"
                
                # Add indicators if we're at limits
                if self.horizontal_spacing_factor >= self.max_horizontal_spacing:
                    h_status += " (max)"
                elif self.horizontal_spacing_factor <= self.min_spacing:
                    h_status += " (min)"
                    
                if self.vertical_spacing_factor >= self.max_vertical_spacing:
                    v_status += " (max)"
                elif self.vertical_spacing_factor <= self.min_spacing:
                    v_status += " (min)"
                
                message = f"Tree spacing: Horizontal {h_status}, Vertical {v_status}"
                main_window.statusbar.showMessage(message, 2000)  # Show for 2 seconds
        except Exception as e:
            # Ignore any errors if status bar update fails
            logger.debug(f"Failed to update status bar: {e}")
            pass
    
    def enter_reroot_mode(self):
        """Enter reroot mode, changing cursor and enabling branch selection for rerooting"""
        if not self.tree:
            return False
            
        self.reroot_mode = True
        # Use a cross hair cursor to indicate selection mode
        self.setCursor(Qt.CursorShape.CrossCursor)
        
        # Update tooltip
        self.setToolTip("Click on a branch to reroot the tree at that point")
        
        # Update status bar
        self.update_status_message("Click on a branch to reroot the tree")
        return True
        
    def exit_reroot_mode(self):
        """Exit reroot mode, restoring normal interaction"""
        self.reroot_mode = False
        # Restore default cursor
        self.setCursor(Qt.CursorShape.OpenHandCursor)
        
        # Restore default tooltip
        self.setToolTip("Click and drag nodes to manipulate the tree.\n"
                       "Click and drag the background to pan the view.\n"
                       "Double-click a leaf node to rename.\n"
                       "Right-click a leaf node for more options.")
        
        # Update status bar
        self.update_status_message("Reroot mode exited")
        return True
        
    def reroot_tree_at_clade(self, clade):
        """Reroot the tree at the specified clade"""
        if not self.tree or not clade:
            logger.error("Cannot reroot: tree or clade is None")
            return False
            
        try:
            # Log the clade we're trying to reroot at
            logger.debug(f"Attempting to reroot tree at clade: {clade}")
            if hasattr(clade, 'name'):
                logger.debug(f"Clade name: {clade.name if clade.name else 'None'}")
            
            # Make sure we're using the right method for the tree type
            if hasattr(self.tree, 'root_with_outgroup'):
                logger.debug("Using root_with_outgroup method")
                self.tree.root_with_outgroup(clade)
            elif hasattr(self.tree, 'reroot_at_node'):
                logger.debug("Using reroot_at_node method")
                self.tree.reroot_at_node(clade)
            else:
                # Fallback to manual rerooting if needed
                logger.debug("No standard reroot method found, attempting manual reroot")
                
                from Bio.Phylo.BaseTree import Clade
                # Create a new root node
                new_root = Clade(branch_length=0.0)
                # Make the current root and the selected clade children of the new root
                new_root.clades.append(self.tree.root)
                new_root.clades.append(clade)
                # Set the new root
                self.tree.root = new_root
            
            # Recalculate layout and update display
            self.calculate_tree_layout()
            self.update()
            
            # Update status
            clade_name = clade.name if hasattr(clade, 'name') and clade.name else "selected branch"
            self.update_status_message(f"Tree rerooted at {clade_name}")
            
            # Exit reroot mode
            self.exit_reroot_mode()
            return True
        except Exception as e:
            import traceback
            logger.error(f"Failed to reroot tree: {e}")
            logger.error(traceback.format_exc())
            self.update_status_message(f"Failed to reroot tree: {str(e)}")
            self.exit_reroot_mode()
            return False
            
    def update_status_message(self, message):
        """Update the main window status bar with a message"""
        try:
            from PyQt6.QtWidgets import QApplication
            main_window = QApplication.activeWindow()
            
            if main_window and hasattr(main_window, 'statusbar'):
                main_window.statusbar.showMessage(message, 3000)  # Show for 3 seconds
        except Exception as e:
            # Ignore any errors if status bar update fails
            logger.debug(f"Failed to update status bar: {e}")
            pass

    def find_clade_at_position(self, pos):
        """Find a clade at the given position"""
        if not self.tree:
            return None

        # Find the closest node with special consideration for leaf node labels
        closest_node = None
        closest_dist = float('inf')
        leaf_label_width = 800  # Approximate width of label text area
        
        # First, check all nodes - for dragging, we want both internal and leaf nodes to be draggable
        # Use our node_positions mapping which contains unique IDs
        for name, (x, y) in self.node_positions.items():
            # Check distance to this node
            dist_to_node = ((pos.x() - x)**2 + (pos.y() - y)**2)**0.5
            
            # Use a larger hit area (15 pixel radius) for better usability
            if dist_to_node < 15:  # Radius for node hover/click
                # If this is closer than any previous match, use it
                if dist_to_node < closest_dist:
                    closest_dist = dist_to_node
                    closest_node = name
                    logger.debug(f"Click on node: {name} at ({x},{y}), click at ({pos.x()},{pos.y()})")
                    # Don't break here - we want the closest node
                
        # Special check for leaf nodes - only for text labels, NOT for hover effects
        # This is separate from the hover checking
        self.over_text_label = False
        self.text_label_node = None
        
        # Check all leaf nodes for text labels - this is only for cursor changes, not hover effects
        for name, (x, y) in self.node_positions.items():
            # Only check leaf nodes for text hovering
            is_leaf = name.startswith("leaf_") or name.startswith("unnamed_leaf_")
            
            if is_leaf:
                # Check if cursor is over the text part of a label
                if pos.x() > x + 10 and pos.x() <= x + leaf_label_width and abs(pos.y() - y) < 15:
                    # Store that we're over text, but don't set a hover node
                    self.over_text_label = True
                    self.text_label_node = name
                    logger.debug(f"Cursor over leaf text, not showing hover: {name}")
                    break

        if closest_node:
            # Enhanced clade finder with better debugging and handling of internal nodes
            def find_clade_by_name(clade, name):
                # Check if we have a direct mapping from node name to clade object
                if hasattr(self, 'object_id_map'):
                    # Invert the map to find clades by their unique node name
                    for obj_id, node_name in self.object_id_map.items():
                        if node_name == name:
                            # Now search for the clade with this ID
                            stack = [self.tree.root]
                            while stack:
                                current = stack.pop()
                                if id(current) == obj_id:
                                    logger.debug(f"Found clade with ID {obj_id} matching node name {name}")
                                    return current
                                
                                if hasattr(current, 'clades'):
                                    for child in current.clades:
                                        stack.append(child)
                
                # Handle legacy names
                if name.startswith("leaf_") and "_" in name:
                    # Extract original name from leaf_Name_Position format
                    parts = name.split("_", 2)
                    if len(parts) >= 3:
                        original_name = parts[1]
                        # Find terminal with this name
                        for terminal in self.tree.get_terminals():
                            if terminal.name == original_name:
                                logger.debug(f"Found leaf node with name {original_name}")
                                return terminal
                
                # Handle internal_named nodes
                if name.startswith("internal_named_") and "_" in name:
                    # Format is internal_named_Name_ID - extract original name
                    parts = name.split("_", 3)
                    if len(parts) >= 4:
                        original_name = parts[2]
                        # Try to match by name
                        stack = [self.tree.root]
                        while stack:
                            current = stack.pop()
                            if current.name == original_name and hasattr(current, 'clades') and current.clades:
                                logger.debug(f"Found internal node with name {original_name}")
                                return current
                            
                            if hasattr(current, 'clades'):
                                for child in current.clades:
                                    stack.append(child)
                
                # Is this a legacy internal node ID?
                if name.startswith("internal_") and '_' in name:
                    # Extract the ID from the name
                    try:
                        clade_id_str = name.split("_", 1)[1]
                        if clade_id_str.isdigit():
                            # This is a old-style numeric ID like internal_123_456
                            logger.debug(f"Looking for internal clade with old-style ID: {name}")
                            # For these, we try to compare by string representation 
                            # which can sometimes help locate the node
                            if str(id(clade)) in name:
                                logger.debug(f"Found internal clade with ID: {name}")
                                return clade
                        else:
                            # This is the new-style direct object ID like internal_139983751998496
                            # Convert the string representation of the ID to an integer
                            clade_id = int(clade_id_str)
                            # Check if this clade's ID matches
                            if id(clade) == clade_id:
                                logger.debug(f"Found exact match for internal clade ID: {clade_id}")
                                return clade
                    except (ValueError, IndexError) as e:
                        # If we can't parse the ID, fall back to name comparison
                        logger.debug(f"Error parsing internal ID '{name}': {e}")
                        pass
                
                # Fallback to regular name comparison
                if clade.name == name:
                    logger.debug(f"Found clade with exact name match: {name}")
                    return clade
                    
                # Recursively check children
                if hasattr(clade, 'clades'):
                    for child in clade.clades:
                        result = find_clade_by_name(child, name)
                        if result:
                            return result
                            
                return None

            # First attempt with the standard function
            result = find_clade_by_name(self.tree.root, closest_node)
            
            # If that didn't work and this is an internal node ID, try a direct recursive search
            if result is None and closest_node.startswith("internal_"):
                logger.debug(f"Standard search failed for {closest_node}, trying direct ID search")
                
                # Try a direct search through all nodes
                stack = [self.tree.root]
                while stack:
                    current = stack.pop()
                    
                    # Check if the string representation of this object's ID is in the node name
                    # This is a fallback approach that sometimes works when exact ID matching fails
                    node_id_str = str(id(current))
                    if node_id_str in closest_node:
                        logger.debug(f"Found node with matching ID fragment: {node_id_str} in {closest_node}")
                        return current
                        
                    # Add all children to the stack
                    if hasattr(current, 'clades'):
                        for child in current.clades:
                            stack.append(child)
                
                # If we still didn't find it, try one more approach - find a node with approximately the same position
                logger.debug("Direct ID search also failed, trying position-based search")
                if closest_node in self.node_positions:
                    target_x, target_y = self.node_positions[closest_node]
                    
                    # Traverse all nodes to find one at a similar position
                    stack = [self.tree.root]
                    best_match = None
                    best_distance = float('inf')
                    
                    while stack:
                        current = stack.pop()
                        current_name = current.name if current.name else f"internal_{id(current)}"
                        
                        # If this node has a position, check how close it is to our target
                        if current_name in self.node_positions:
                            curr_x, curr_y = self.node_positions[current_name]
                            dist = ((curr_x - target_x)**2 + (curr_y - target_y)**2)**0.5
                            
                            if dist < best_distance:
                                best_distance = dist
                                best_match = current
                                
                        # Add children to stack
                        if hasattr(current, 'clades'):
                            for child in current.clades:
                                stack.append(child)
                    
                    # If we found a close match, use it
                    if best_match and best_distance < 5.0:  # 5.0 pixels threshold
                        logger.debug(f"Using position-based match with distance {best_distance}")
                        return best_match
            
            # Return the result of our search
            return result

        return None

            
    def swap_clades(self, clade):
        """Swap the order of child clades for the given clade"""
        if clade and len(clade.clades) > 1:
            clade.clades = list(reversed(clade.clades))
            self.update()
            return True
        return False
    
    def rotate_clade(self, clade):
        """Rotate child clades in the given clade - for clades with more than 2 children"""
        if clade and len(clade.clades) > 2:
            # Move the first clade to the end
            clade.clades = clade.clades[1:] + [clade.clades[0]]
            self.update()
            return True
        elif clade and len(clade.clades) == 2:
            return self.swap_clades(clade)
        return False
    
    def find_parent_clade(self, clade, target_clade):
        """Find the parent clade of the target clade"""
        if clade == target_clade:
            return None
        
        for child in clade.clades:
            if child == target_clade:
                return clade
                
            result = self.find_parent_clade(child, target_clade)
            if result:
                return result
                
        return None
        
    def set_branch_width(self, width):
        """Set the width of branch lines
        
        Args:
            width: Line width in pixels (1-5)
        """
        self.branch_width = max(1, min(int(width), 5))
        # Force redraw with new line width
        self.update()
    
    def draw_clade(self, painter, clade, parent_x=None, parent_y=None):
        """Draw a clade in the tree"""
        # Safety check - validate clade
        if clade is None:
            logger.error("Attempted to draw None clade")
            return
            
        # Define node and branch colors based on dark mode
        if hasattr(self, '_dark_mode') and self._dark_mode:
            branch_color = QColor(255, 255, 255)    # White branches for dark mode
            node_color = QColor(0, 120, 255)        # Blue nodes for dark mode
        else:
            branch_color = QColor(0, 0, 0)          # Black branches for light mode
            node_color = QColor(0, 80, 220)         # Darker blue nodes for light mode
        
        # Get the node name from our object ID mapping
        clade_id = id(clade)
        if hasattr(self, 'object_id_map') and clade_id in self.object_id_map:
            # Use the mapping created during layout
            clade_name = self.object_id_map[clade_id]
            logger.debug(f"Drawing: Found clade {clade_id} in object_id_map as {clade_name}")
        else:
            # Fallback to legacy naming method
            if not hasattr(clade, 'name'):
                logger.error(f"Clade {clade_id} has no name attribute")
                clade_name = f"internal_{clade_id}"
            else:
                # Check if it's a leaf node
                if not clade.clades:  # Terminal/leaf node
                    if clade.name:
                        # For leaf nodes, we need to search for the right position
                        # Try to find this leaf in the node_positions dictionary
                        leaf_key = None
                        for key in self.node_positions.keys():
                            if key.startswith(f"leaf_{clade.name}_"):
                                leaf_key = key
                                break
                        
                        if leaf_key:
                            clade_name = leaf_key
                            logger.debug(f"Drawing: Found leaf node {clade.name} as {clade_name}")
                        else:
                            # Fallback to legacy naming
                            clade_name = clade.name
                    else:
                        # For unnamed leaf nodes, try to search by pattern
                        leaf_key = None
                        for key in self.node_positions.keys():
                            if key.startswith(f"unnamed_leaf_"):
                                # Check if this position matches
                                leaf_key = key
                                break
                        
                        if leaf_key:
                            clade_name = leaf_key
                        else:
                            clade_name = "unnamed"
                else:
                    # Internal node
                    if clade.name:
                        # For named internal nodes, search for the pattern
                        internal_key = None
                        for key in self.node_positions.keys():
                            if key.startswith(f"internal_named_{clade.name}_"):
                                internal_key = key
                                break
                            
                        if internal_key:
                            clade_name = internal_key
                            logger.debug(f"Drawing: Found internal node {clade.name} as {clade_name}")
                        else:
                            # Fallback to legacy naming
                            clade_name = f"internal_{clade_id}"
                    else:
                        clade_name = f"internal_{clade_id}"
        
        # Debug output for drawing
        logger.debug(f"Drawing: Looking for position of {clade_name}")
            
        # Get the coordinates for this clade
        if clade_name not in self.node_positions:
            logger.warning(f"No position for clade: {clade_name}")
            return

        # Verify we have valid positions for the node
        try:
            x, y = self.node_positions[clade_name]
            
            # Validate coordinates are numbers
            if not (isinstance(x, (int, float)) and isinstance(y, (int, float))):
                logger.error(f"Invalid coordinates for clade {clade_name}: x={x}, y={y}")
                return
        except Exception as e:
            logger.error(f"Error getting position for clade {clade_name}: {e}")
            return

        # Convert coordinates to integers to avoid QPainter errors
        x_int = int(x)
        y_int = int(y)
        parent_x_int = int(parent_x) if parent_x is not None else None
        parent_y_int = int(parent_y) if parent_y is not None else None

        # Define marker size for all node types
        marker_size = 8
        
        # Draw node markers only for internal nodes to make them more visible and draggable
        if hasattr(clade, 'clades') and clade.clades:
            marker_rect = QRect(int(x_int - marker_size/2), int(y_int - marker_size/2), 
                                marker_size, marker_size)
            
            # Internal node - filled circle
            painter.setPen(Qt.PenStyle.NoPen)  # No border for cleaner look
            painter.setBrush(node_color)  # Blue fill from our earlier definition
            
            # Draw the node marker to make it clear what's draggable
            painter.drawEllipse(marker_rect)
        
        # Add EXTRA highlighting for hovering - brighter and larger
        # Show hover highlight for ALL nodes to make dragging more intuitive
        if clade == self.hover_node:
            # Make highlight size larger than the node
            highlight_size = marker_size + 4
            
            # Draw a highlight around the node - make it very noticeable
            # Make sure all values are integers for QRect
            hs_half = int(highlight_size/2)
            highlight_rect = QRect(x_int - hs_half, y_int - hs_half, 
                                   int(highlight_size), int(highlight_size))
            
            # Use a bright highlight color that works in both light and dark modes
            highlight_color = QColor(50, 180, 255)  # Bright blue
            painter.setPen(QPen(highlight_color, 2))
            painter.setBrush(QColor(highlight_color.red(), highlight_color.green(), highlight_color.blue(), 100))
            painter.drawEllipse(highlight_rect)

        # Draw branch from parent (if not root)
        if parent_x_int is not None and parent_y_int is not None:
            # Generate a consistent key for checking highlights
            highlight_key = clade.name if clade.name else f"internal_{id(clade)}"
            
            # Set pen color based on highlighting
            if highlight_key in self.highlighted_branches:
                color = self.highlighted_branches[highlight_key]
                pen = QPen(color, self.branch_width + 1)  # Make highlighted branches slightly thicker
                painter.setPen(pen)
            else:
                # Use the branch_color we defined earlier in paintEvent
                painter.setPen(QPen(branch_color, self.branch_width))

            # Draw horizontal line from parent's level to this node's level
            painter.drawLine(parent_x_int, parent_y_int, parent_x_int, y_int)

            # Draw horizontal line to this node - make sure x_int is different from parent_x_int
            # This ensures nodes are actually drawn at different x-coordinates
            if x_int <= parent_x_int:
                # Fix the issue with flat tree by ensuring x_int is always greater than parent_x_int
                # Calculate a proper offset based on tree depth
                # Get the current horizontal spacing factor and use it to calculate the offset
                x_offset = max(50, 100 * self.horizontal_spacing_factor)
                x_int = parent_x_int + int(x_offset)
                # Update the node position for future reference
                self.node_positions[clade_name] = (x_int, y_int)
                
            # Ensure we're actually drawing a meaningful horizontal line
            if abs(x_int - parent_x_int) < 10:
                # If nodes are too close horizontally, add minimum distance
                x_int = parent_x_int + 50  # Minimum horizontal branch length
                self.node_positions[clade_name] = (x_int, y_int)
                
            painter.drawLine(parent_x_int, y_int, x_int, y_int)

        
        # Draw bootstrap value if available and enabled
        if (parent_x_int is not None and self.show_bootstrap and hasattr(clade, 'confidence') and clade.confidence is not None):
            # Bootstrap values can be stored in different formats in different tree files
            # Handle both integer and float values
            try:
                # Try to convert to float first to handle all numeric formats
                conf_val = float(clade.confidence)
                
                # Format confidence value - use integer if it's a whole number, otherwise show decimal
                if conf_val == int(conf_val):
                    confidence = int(conf_val)
                    display_val = str(confidence)
                else:
                    # It's a fraction or percentage - keep one decimal place
                    confidence = conf_val
                    display_val = f"{confidence:.1f}"
                
                # Only display if the value is greater than 0
                if confidence > 0:
                    # Calculate position for bootstrap value
                    # Place it near the branch point
                    bs_x = (parent_x_int + x_int) // 2
                    bs_y = y_int - 10
                    
                    # Determine color threshold - if values are 0-1 scale, use 0.7 as threshold
                    # if values are 0-100 scale, use 70 as threshold
                    color_threshold = 0.7 if conf_val <= 1 else 70
                    
                    # Draw the bootstrap value with appropriate color
                    if self.dark_mode:
                        # Use color based on value (green for high confidence)
                        if confidence >= color_threshold:
                            painter.setPen(QPen(QColor(0, 255, 0), 1))  # Green
                        else:
                            painter.setPen(QPen(QColor(255, 165, 0), 1))  # Orange
                    else:
                        if confidence >= color_threshold:
                            painter.setPen(QPen(QColor(0, 100, 0), 1))  # Dark green
                        else:
                            painter.setPen(QPen(QColor(200, 100, 0), 1))  # Dark orange
                    
                    painter.drawText(bs_x, bs_y, display_val)
            except (ValueError, TypeError) as e:
                # If conversion fails, try using the value as-is
                logger.warning(f"Could not convert bootstrap value to number: {clade.confidence}, error: {e}")
                # Just use the confidence value as a string if it's not numeric
                if str(clade.confidence).strip():
                    bs_x = (parent_x_int + x_int) // 2
                    bs_y = y_int - 10
                    painter.setPen(QPen(QColor(150, 150, 150), 1))  # Gray for non-numeric values
                    painter.drawText(bs_x, bs_y, str(clade.confidence))

        # Draw label for leaf nodes
        if not clade.clades:
            # Draw full label - preserve spaces in names
            # The node name should already have full description with spaces
            # Use the original name from node_name_map if available
            if hasattr(self, 'node_name_map') and clade_name in self.node_name_map:
                full_label = self.node_name_map[clade_name]
                logger.debug(f"Drawing label using node_name_map: {full_label}")
            else:
                full_label = clade.name if clade.name else "Unnamed"
                logger.debug(f"Drawing label using clade.name: {full_label}")
            
            # Save original font for restoration later
            original_font = painter.font()
            
            # Import QFont to ensure it's available
            from PyQt6.QtGui import QFont
            
            # Get the current sequence_font
            current_font = None
            
            # Create font for leaf node labels - based on the user-selected font
            # Debug the font being used
            if hasattr(self, 'sequence_font'):
                current_font = self.sequence_font
                logger.debug(f"Drawing leaf node '{clade.name}' with sequence_font type: {type(current_font)}")
            else:
                logger.debug(f"Drawing leaf node '{clade.name}' with no sequence_font attribute")
            
            # Make sure we're using a valid QFont object
            try:
                if current_font is not None and isinstance(current_font, QFont):
                    # Extract font properties and create a brand new font to avoid any reference issues
                    family = current_font.family()
                    size = current_font.pointSize()
                    is_bold = current_font.bold()
                    is_italic = current_font.italic()
                    
                    # Create a completely new font
                    leaf_font = QFont()
                    leaf_font.setFamily(family)
                    leaf_font.setPointSize(size if size > 0 else 12)  # Default to 12pt if size is invalid
                    leaf_font.setBold(is_bold)
                    leaf_font.setItalic(is_italic)
                    
                    logger.debug(f"Created new leaf font: {leaf_font.family()}, {leaf_font.pointSize()}pt, Bold: {leaf_font.bold()}")
                else:
                    # Fallback to a default font if sequence_font is not a valid QFont
                    leaf_font = QFont("Arial", 10)
                    logger.debug("Using fallback font: Arial, 10pt")
            except Exception as e:
                # Handle any errors and use a safe fallback
                logger.error(f"Error creating font for leaf node: {e}")
                leaf_font = QFont("Arial", 10)
                logger.debug("Using fallback font after error: Arial, 10pt")
            
            # Check if this node is being hovered in delete mode
            is_delete_hover = False
            if hasattr(self, 'delete_mode') and self.delete_mode and hasattr(self, 'hovered_sequence'):
                # Convert clade_name to match the format in hovered_sequence if needed
                node_key = clade_name
                if self.hovered_sequence == node_key:
                    is_delete_hover = True
            
            # Determine text appearance based on delete mode and hover state
            if is_delete_hover:
                # In delete mode with hover, use inverted colors (white on red)
                hover_font = QFont(leaf_font)
                hover_font.setBold(True)  # Make text bold for emphasis
                painter.setFont(hover_font)
                
                # Draw a highlight rectangle with delete/danger colors
                label_rect = QRect(x_int + 5, y_int - 15, 800, 25)
                if self.dark_mode:
                    # Dark mode: Use bright red background with white text
                    highlight_color = QColor(200, 0, 0)  # Bright red background
                    text_color = QColor(255, 255, 255)  # White text
                else:
                    # Light mode: Use medium red background with white text
                    highlight_color = QColor(180, 0, 0)  # Medium red background
                    text_color = QColor(255, 255, 255)  # White text
                
                # Fill the background and set text color
                painter.fillRect(label_rect, highlight_color)
                painter.setPen(text_color)
            elif clade == self.hover_node:
                # Normal hover effect (not in delete mode)
                hover_font = QFont(leaf_font)
                hover_font.setBold(True)
                painter.setFont(hover_font)
                
                # Draw a light highlight behind the text for better visibility
                label_rect = QRect(x_int + 5, y_int - 15, 800, 25)
                highlight_color = QColor(100, 200, 255, 70)  # Light blue with transparency
                painter.fillRect(label_rect, highlight_color)
            else:
                # No hover - use the customized font
                painter.setFont(leaf_font)
            
            # Calculate font metrics for proper vertical positioning
            font_metrics = QFontMetrics(painter.font())
            font_height = font_metrics.height()
            
            # Adjust text rectangle based on font size
            y_offset = max(-10, -font_height // 2)
            height = max(30, font_height + 10)
            
            # Draw the label with potential wrapping for long names
            # Increase the width of the text rectangle to accommodate longer sequence names
            # Use a much wider text area to ensure the full name is displayed
            text_rect = QRect(x_int + 5, y_int + y_offset, 800, height)
            
            # Use correct parameters for drawText with QRect
            painter.drawText(text_rect, Qt.AlignmentFlag.AlignLeft, full_label)
            
            # Restore original pen color if we changed it for delete hover
            if is_delete_hover:
                if self.dark_mode:
                    painter.setPen(QColor(255, 255, 255))  # Reset to white for dark mode
                else:
                    painter.setPen(QColor(0, 0, 0))  # Reset to black for light mode
            
            # Restore original font
            painter.setFont(original_font)

        # Draw clade label if it exists
        if clade.name in self.clade_labels:
            painter.setFont(QFont("Arial", 10, QFont.Weight.Bold))
            painter.drawText(x_int, y_int - 10, self.clade_labels[clade.name])
            painter.setFont(QFont("Arial", 10))

        # Draw child clades with error handling
        try:
            if hasattr(clade, 'clades'):
                for child in clade.clades:
                    if child is not None:
                        self.draw_clade(painter, child, x_int, y_int)
                    else:
                        logger.error("Encountered None child clade")
            else:
                logger.warning(f"Clade {clade_name} has no clades attribute")
        except Exception as e:
            logger.error(f"Error drawing child clades for {clade_name}: {e}")

    def draw_tree(self, painter):
        """Draw the phylogenetic tree on the given painter
        
        This method handles all tree drawing operations and can be used
        both for on-screen painting and for image export.
        
        Args:
            painter: QPainter object to draw with
        """
        from PyQt6.QtGui import QFont
        
        # Fill the background based on theme
        if hasattr(self, '_dark_mode') and self._dark_mode:
            painter.fillRect(self.rect(), QColor(40, 40, 40))  # Dark background
            text_color = QColor(255, 255, 255)      # White text for dark mode
            branch_color = QColor(255, 255, 255)    # White branches for dark mode
            node_color = QColor(0, 120, 255)        # Blue nodes for dark mode
        else:
            painter.fillRect(self.rect(), QColor(255, 255, 255))  # White background
            text_color = QColor(0, 0, 0)            # Black text for light mode
            branch_color = QColor(0, 0, 0)          # Black branches for light mode
            node_color = QColor(0, 80, 220)         # Darker blue nodes for light mode

        painter.setPen(QPen(text_color, 1))
        
        # Handle case where no tree is available
        if not self.tree:
            # Set a larger font for the message
            message_font = QFont("Arial", 12)
            painter.setFont(message_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, 
                            "No tree available\nUse 'Build Tree' from the Tree menu to create a tree")
            return

        # Log current sequence font for debugging (but less frequently)
        if hasattr(self, 'sequence_font'):
            if isinstance(self.sequence_font, QFont):
                logger.debug(f"TreeCanvas.draw_tree - Using sequence font: {self.sequence_font.family()}, {self.sequence_font.pointSize()}pt")
            else:
                logger.debug(f"TreeCanvas.draw_tree - sequence_font is not a QFont: {type(self.sequence_font)}")
        else:
            logger.debug("TreeCanvas.draw_tree - No sequence_font attribute found")

        # Calculate tree layout if needed
        if not self.node_positions:
            self.calculate_tree_layout()

        # Draw the tree starting from root
        if self.tree.root:
            # Draw the tree root
            self.draw_clade(painter, self.tree.root)
        else:
            # Draw message if tree exists but root is missing
            logger.debug("Tree exists but has no root")
            message_font = QFont("Arial", 12)
            painter.setFont(message_font)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, 
                            "Tree has no root node\nPlease try rebuilding the tree")

    def paintEvent(self, event):
        """Draw the phylogenetic tree"""
        # Create painter for the widget
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        try:
            # Use the shared drawing method
            self.draw_tree(painter)
        except Exception as e:
            # Draw error message if something goes wrong
            painter.resetTransform()  # Reset transform for error message
            # Create font directly without variable
            painter.setFont(QFont("Arial", 12))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, 
                            f"Error drawing tree: {str(e)}\n\nPlease try rebuilding the tree")
            import traceback
            logger.error(f"Error in tree painting: {e}")
            logger.error(traceback.format_exc())

        finally:
            # Ensure painter is properly ended
            painter.end()

    def calculate_tree_layout(self):
        """Calculate the layout positions for the tree"""
        if not self.tree:
            logger.debug("No tree to layout")
            return

        try:
            # Log entry to this method for debugging
            logger.debug("Calculating tree layout...")
            
            # Validate the tree has the needed structure for layout
            if not hasattr(self.tree, 'root'):
                logger.error("Tree has no root - cannot calculate layout")
                return
                
            if not self.tree.root:
                logger.error("Tree root is None or empty - cannot calculate layout")
                return
                
            if not hasattr(self.tree.root, 'clades'):
                logger.error("Tree root has no clades attribute - cannot calculate layout")
                return
                
            if not self.tree.root.clades:
                logger.warning("Tree root has empty clades list - tree may not display correctly")
            
            # Calculate layout for tree
            
            # Clear existing positions and mappings to start fresh
            self.node_positions = {}
            self.node_name_map = {}  # Maps node IDs to original names
            self.object_id_map = {}  # Maps object IDs to node names

            # Get tree dimensions
            margin = 50

            # Calculate max label width to reserve space for the full names
            font_metrics = QFontMetrics(QFont("Arial", 9))  # Same font used for labels
            max_label_width = 0
            terminals = list(self.tree.get_terminals())
            # Process terminal nodes for layout
            
            for clade in terminals:
                # Use full label including spaces
                label = clade.name if clade.name else "Unnamed"
                
                # For longer labels, calculate a reasonable width
                # Most tree views truncate or wrap very long labels
                if len(label) > 50:
                    # Use a reasonable max width for very long labels
                    width = font_metrics.horizontalAdvance(label[:50] + "...") + 20
                else:
                    # For normal length labels, use the full width
                    width = font_metrics.horizontalAdvance(label) + 20
                    
                max_label_width = max(max_label_width, width)
            
            # Reserve extra space for labels to accommodate full names with spaces
            label_margin = max(max_label_width + 50, 800)  # Ensure at least 800 pixels for labels

            # Adjust available drawing space
            width = self.width() - 2 * margin - label_margin
            height = self.height() - 2 * margin
            # Calculate drawing area dimensions

            # Count leaf nodes to determine vertical spacing
            leaf_count = len(terminals)

            if leaf_count == 0:
                logger.debug("No leaf nodes found")
                return

            # Calculate vertical spacing with the scaling factor
            # Ensure we have enough space for all leaves, with at least 20 pixels per leaf
            min_spacing = 20  # Minimum pixels between leaf nodes
            calculated_spacing = (height / leaf_count)
            
            # Use the larger of calculated spacing or minimum spacing
            base_spacing = max(calculated_spacing, min_spacing)
            
            # Apply vertical spacing factor
            vertical_spacing = base_spacing * self.vertical_spacing_factor
            
            # Adjust the widget height if needed for large trees
            required_height = (vertical_spacing * leaf_count) + (2 * margin)
            
            # Cap the maximum height to prevent Qt widget size errors
            # Qt has a maximum widget dimension of 16777215 pixels, but we'll stay much lower
            # Calculate a reasonable max height based on the number of sequences
            if leaf_count > 40:
                max_allowed_height = 3000  # Stricter limit for very large trees
            elif leaf_count > 20:
                max_allowed_height = 4000  # Moderate limit for medium trees 
            else:
                max_allowed_height = 5000  # More space for small trees
            
            if required_height > self.height() and required_height < max_allowed_height:
                # Resize the widget to fit all leaves, but not beyond max allowed height
                new_height = int(min(required_height, max_allowed_height))
                logger.debug(f"Resizing tree canvas to height: {new_height}")
                
                # Only set new height if it's reasonable
                if new_height > 100 and new_height < max_allowed_height:
                    self.setMinimumHeight(new_height)
                    
                    # If parent is a scroll area, tell it to update
                    parent = self.parent()
                    if parent and hasattr(parent, 'updateTreeCanvasSize'):
                        # Update virtual size in parent
                        if hasattr(parent, 'virtual_size'):
                            current_width = parent.virtual_size.width()
                            # Keep virtual size reasonable
                            parent.virtual_size = QSize(current_width, min(new_height + 100, max_allowed_height))
                            parent.updateTreeCanvasSize()

            # Get tree depth (for horizontal scaling)
            def get_tree_depth(clade):
                if not clade.clades:
                    return 0
                return 1 + max(get_tree_depth(child) for child in clade.clades)

            tree_depth = get_tree_depth(self.tree.root)
            if tree_depth == 0:
                tree_depth = 1  # Avoid division by zero
                
            # Use tree depth for horizontal positioning

            # Calculate layout recursively
            current_leaf = [0]  # Using a list to make it mutable in nested function

            def layout_clade(clade, x_offset, depth):
                if not clade.clades:  # Terminal node (leaf)
                    # For leaf nodes (terminals), use position index as part of the node ID
                    # to handle duplicate names
                    y = margin + current_leaf[0] * vertical_spacing
                    # Apply horizontal spacing factor to leaf nodes too
                    # Ensure all leaf nodes are at the furthest depth for proper tree visualization
                    # This forces leaf nodes to be at the right edge of the tree
                    min_x_per_level = 120  # Minimum pixels between levels (same as for internal nodes)
                    
                    # For better positioning of leaf nodes, make sure they're at a consistent x-position
                    # For shallow trees, use direct level-based positioning
                    if tree_depth <= 3:
                        # For very shallow trees, force leaf nodes to the furthest level
                        x = margin + (tree_depth * min_x_per_level) * self.horizontal_spacing_factor
                    else:
                        # Position all leaves at the same rightmost position
                        # Adjust the position based on the compression level
                        if self.horizontal_spacing_factor < 0.5:
                            # When highly compressed, position leaves closer to the root
                            # Scale from 75% at 0.5 down to 55% at 0.15 (min_spacing)
                            width_percent = 0.55 + ((self.horizontal_spacing_factor - 0.15) / 0.35) * 0.2
                            x = margin + (width * width_percent) * self.horizontal_spacing_factor
                        else:
                            # Normal positioning at 85% of available width
                            x = margin + (width * 0.85) * self.horizontal_spacing_factor
                    
                    # Add position index to the node ID to make it unique
                    leaf_position = current_leaf[0]
                    current_leaf[0] += 1
                    
                    # Create a unique key that includes position information for leaf nodes
                    # This ensures nodes with the same name get different positions
                    if clade.name:
                        # Format: "leaf_{name}_{position}" - this makes it unique even with duplicate names
                        node_name = f"leaf_{clade.name}_{leaf_position}"
                        
                        # Store the original name for lookup and display purposes
                        if not hasattr(self, 'node_name_map'):
                            self.node_name_map = {}
                        # Use the full original name with all details from the clade
                        self.node_name_map[node_name] = clade.name
                    else:
                        node_name = f"unnamed_leaf_{leaf_position}"
                    
                    self.node_positions[node_name] = (x, y)
                    # Store reverse mapping from object ID to node name
                    if not hasattr(self, 'object_id_map'):
                        self.object_id_map = {}
                    self.object_id_map[id(clade)] = node_name
                    
                    logger.debug(f"Layout: Set leaf position for {node_name} to ({x}, {y})")
                    return (x, y)

                # Position internal nodes
                child_coords = []
                for child in clade.clades:
                    child_coords.append(layout_clade(child, x_offset, depth + 1))

                # Center internal nodes between their children
                min_y = min(y for _, y in child_coords)
                max_y = max(y for _, y in child_coords)
                y = (min_y + max_y) / 2
                
                # Apply horizontal spacing factor to adjust branch lengths
                # Ensure clear separation between tree levels for internal nodes
                min_x_per_level = 120  # Minimum pixels between tree levels - increased for better visibility
                
                # Calculate positions ensuring proper tree depth display
                # For shallow trees, use direct proportional positioning
                if tree_depth <= 3:
                    # For very shallow trees, use direct level-based positioning
                    x = margin + (depth * min_x_per_level) * self.horizontal_spacing_factor
                else:
                    # For deeper trees, use a combination approach
                    depth_multiplier = max(0.5, depth / max(tree_depth, 1))  # Avoid division by zero
                    x = margin + (width * depth_multiplier) * self.horizontal_spacing_factor
                
                # Ensure nodes at each depth have minimum spacing from previous level
                # Adjust minimum spacing based on compression factor
                if self.horizontal_spacing_factor < 0.5:
                    # When highly compressed, allow nodes to be much closer
                    compression_adjusted_min_x_per_level = min_x_per_level * (0.5 + self.horizontal_spacing_factor)
                    min_x = margin + (depth * compression_adjusted_min_x_per_level)
                else:
                    # Normal minimum spacing
                    min_x = margin + (depth * min_x_per_level * max(0.5, self.horizontal_spacing_factor))
                
                x = max(x, min_x)

                # Use a consistent unique identifier for internal nodes
                # This is a key issue - we need to use the same ID in both layout and drawing
                if clade.name:
                    node_name = f"internal_named_{clade.name}_{id(clade)}"
                    # Store the original name for lookup and display purposes
                    if not hasattr(self, 'node_name_map'):
                        self.node_name_map = {}
                    self.node_name_map[node_name] = clade.name
                else:
                    # Use Python's built-in id function to get a stable, unique identifier for the object
                    node_id = id(clade)
                    node_name = f"internal_{node_id}"
                
                self.node_positions[node_name] = (x, y)
                # Store reverse mapping from object ID to node name
                if not hasattr(self, 'object_id_map'):
                    self.object_id_map = {}
                self.object_id_map[id(clade)] = node_name
                
                # Debug output for important node position mapping
                logger.debug(f"Layout: Set position for {node_name} to ({x}, {y})")
                return (x, y)

            # Calculate the layout
            # Calculate node positions
            layout_clade(self.tree.root, 0, 0)
            
        except Exception as e:
            import traceback
            logger.error(f"Error in tree layout calculation: {e}")
            traceback.print_exc()
            
    def get_canvas_width(self):
        """Returns the current width of the canvas"""
        return self.width()
        
    def get_canvas_height(self):
        """Returns the current height of the canvas"""
        return self.height()
        
    def direct_delete_node(self, node_name):
        """Direct method to handle node deletion when signal chain isn't working
        
        This is a fallback method that will try to find the MainWindow and call its
        delete_tree_node method directly. It includes improved name handling for better
        matching between tree nodes and sequences.
        """
        import logging
        logger = logging.getLogger("treecraft")
        logger.debug(f"Direct delete method called for node: {node_name}")
        
        # First check if we have a valid tree - if not, exit early
        if not self.tree:
            logger.debug("No tree available for delete operation, ignoring")
            return
        
        try:
            # Clean up the node name for better matching with sequences
            # Handle potential underscore or space-separated duplicates
            if '_' in node_name:
                parts = node_name.split('_')
                if len(parts) > 1 and parts[0] == parts[1]:
                    # Handle duplicate parts like "KX897432.1_KX897432.1"
                    clean_name = parts[0]
                    logger.debug(f"Cleaned up duplicated name parts: {node_name} -> {clean_name}")
                    node_name = clean_name
                elif len(parts) > 2 and parts[0] == "leaf":
                    # Handle leaf_name_position format
                    clean_name = parts[1]  # Take just the name part
                    logger.debug(f"Extracted name from leaf format: {node_name} -> {clean_name}")
                    node_name = clean_name
            
            # Handle space-separated duplicates like "KX897432.1 KX897432.1"
            elif ' ' in node_name:
                words = node_name.split()
                # If same word appears multiple times, keep only unique parts
                if len(set(words)) < len(words):
                    clean_name = ' '.join(dict.fromkeys(words))
                    logger.debug(f"Cleaned up duplicated words in name: {node_name} -> {clean_name}")
                    node_name = clean_name
            
            # Try to find the main window by traversing the parent hierarchy
            window = self
            while window and not hasattr(window, 'delete_tree_node'):
                window = window.parent()
                
            if window and hasattr(window, 'delete_tree_node'):
                logger.debug(f"Found main window, calling delete_tree_node directly with name: {node_name}")
                
                # Get better terminal info from the tree for more accurate matching
                # BUT only use exact matches to prevent deleting similarly named sequences
                best_match = node_name
                exact_match_found = False
                
                # Make sure tree and get_terminals are available
                if self.tree and hasattr(self.tree, 'get_terminals'):
                    for terminal in self.tree.get_terminals():
                        term_name = getattr(terminal, 'name', '')
                        # Only use exact matches or node name map matches
                        if term_name == node_name:
                            best_match = term_name
                            exact_match_found = True
                            logger.debug(f"Found exact terminal match: {term_name}")
                            break
                    
                    # If we couldn't find an exact match, check if this might be a node with position suffix
                    if not exact_match_found and '_' in node_name:
                        parts = node_name.split('_')
                        if len(parts) > 1:
                            # Try matching just the first part (without position suffix)
                            base_name = parts[0]
                            for terminal in self.tree.get_terminals():
                                term_name = getattr(terminal, 'name', '')
                                if term_name == base_name:
                                    best_match = term_name
                                    exact_match_found = True
                                    logger.debug(f"Found exact base name match: {term_name}")
                                    break
                else:
                    logger.debug("No valid tree or get_terminals available")
                    return
                
                # Only proceed with deletion if we have an exact match
                if not exact_match_found:
                    logger.warning(f"No exact match found for node '{node_name}' - not deleting to avoid errors")
                    return
                
                if best_match != node_name:
                    logger.debug(f"Using better terminal match: {node_name} -> {best_match}")
                
                # Set a flag in the main window to suppress error messages
                # This is useful for RAxML trees where the name matching can be tricky
                # but the deletion still works correctly
                if hasattr(window, 'suppress_delete_errors'):
                    old_suppress = window.suppress_delete_errors
                    window.suppress_delete_errors = True
                else:
                    old_suppress = None
                    window.suppress_delete_errors = True
                
                try:
                    # Try deletion with the improved name
                    window.delete_tree_node(best_match)
                finally:
                    # Restore the previous suppress_delete_errors state
                    if old_suppress is not None:
                        window.suppress_delete_errors = old_suppress
                    else:
                        delattr(window, 'suppress_delete_errors')
                
                # Force redraw after deletion attempt
                self.update()
            else:
                logger.error(f"Could not find main window with delete_tree_node method")
        except Exception as e:
            # Catch any exceptions to prevent crashes
            logger.error(f"Error in TreeCanvas.direct_delete_node: {e}")
            import traceback
            logger.error(traceback.format_exc())
            
            # Force redraw of the canvas regardless
            self.update()
        
    def find_parent_clade(self, root_clade, target_clade):
        """Find the parent clade of the target clade in the tree
        
        Args:
            root_clade: The root clade to start searching from
            target_clade: The clade whose parent we want to find
            
        Returns:
            The parent clade if found, None otherwise
        """
        if root_clade == target_clade:
            return None
            
        # Check if any of the direct children is the target
        if hasattr(root_clade, 'clades'):
            for child in root_clade.clades:
                if child == target_clade:
                    return root_clade
                    
                # Recursive search
                parent = self.find_parent_clade(child, target_clade)
                if parent:
                    return parent
                    
        return None
        
    def swap_clades(self, clade, drag_clade=None, direction=None):
        """Swap the order of child clades for the given clade
        
        Args:
            clade: The clade whose children will be swapped
            drag_clade: The specific child clade being dragged (optional)
            direction: Direction of drag (negative = left, positive = right) (optional)
            
        Returns:
            True if the swap was successful, False otherwise
        """
        import logging
        logger = logging.getLogger("treecraft")
        logger.debug(f"Swapping clades in {clade} with {len(getattr(clade, 'clades', []))} children")
        
        if clade and hasattr(clade, 'clades') and len(clade.clades) > 1:
            # Simple approach - just reverse all children
            # This is more reliable than trying targeted swaps
            logger.debug("Reversing all children")
            clade.clades = list(reversed(clade.clades))
            
            # Force tree recalculation and update
            try:
                self.calculate_tree_layout()
            except Exception as e:
                logger.error(f"Error recalculating tree layout: {e}")
                
            # Force immediate redraw
            self.update()
            return True
            
        logger.debug(f"Could not swap: clade has {len(getattr(clade, 'clades', []))} children")
        return False
        
    def rotate_clade(self, clade, direction=None):
        """Rotate child clades in the given clade
        
        Args:
            clade: The clade whose children will be rotated
            direction: Optional direction of rotation (negative = counterclockwise, positive = clockwise)
            
        Returns:
            True if the rotation was successful, False otherwise
        """
        import logging
        logger = logging.getLogger("treecraft")
        logger.debug(f"Rotating clade with {len(getattr(clade, 'clades', []))} children")
        
        if clade and hasattr(clade, 'clades'):
            if len(clade.clades) > 2:
                # Simple approach - always rotate counter-clockwise
                # Move first clade to the end - this is more reliable
                logger.debug("Rotating counter-clockwise (standard rotation)")
                clade.clades = clade.clades[1:] + [clade.clades[0]]
                
                # Force tree recalculation and update
                try:
                    self.calculate_tree_layout()
                except Exception as e:
                    logger.error(f"Error recalculating tree layout: {e}")
                
                # Force immediate redraw
                self.update()
                return True
                
            elif len(clade.clades) == 2:
                # For binary nodes, rotation is the same as swapping
                logger.debug("Binary node - swapping children")
                return self.swap_clades(clade)
                
        logger.debug(f"Could not rotate: clade has {len(getattr(clade, 'clades', []))} children")
        return False


class PhyloCanvas(QScrollArea):
    """Scrollable canvas with zoom controls for the phylogenetic tree"""
    
    # Define signals for node interactions
    rename_node_signal = pyqtSignal(str)  # Signal to rename a node
    delete_node_signal = pyqtSignal(str)  # Signal to delete a node
    
    # Add signals that main_window.py is expecting
    node_double_clicked = pyqtSignal(str)
    node_right_clicked = pyqtSignal(str, QPoint)
    rename_signal = pyqtSignal(str)  # Signal for renaming node
    delete_signal = pyqtSignal(str)  # Signal for deleting node
    update_status = pyqtSignal(str)  # Signal for status bar updates
    
    def __init__(self, parent=None):
        super().__init__(parent)
        # Initialize private attributes first
        self._dark_mode = True
        self._sequence_font = None
        self.delete_mode = False  # Flag to track delete mode
        
        # Create the tree canvas first
        self.tree_canvas = TreeCanvas()
        self.setWidget(self.tree_canvas)
        
        # Initialize the virtual size (used for scaling)
        from PyQt6.QtCore import QSize
        self.virtual_size = QSize(1000, 800)  # Default virtual size
        
        # Set tree_canvas dark_mode directly without using the property setter
        if hasattr(self.tree_canvas, 'dark_mode'):
            self.tree_canvas.dark_mode = self._dark_mode
        
        # Connect signals for node interactions
        # Properly connect the TreeCanvas signals to our PhyloCanvas signals
        self.tree_canvas.rename_signal.connect(self.rename_signal.emit)
        self.tree_canvas.delete_signal.connect(self.delete_signal.emit)
        self.tree_canvas.rename_signal.connect(self.rename_node_signal.emit)
        self.tree_canvas.delete_signal.connect(self.delete_node_signal.emit)
        
        # Connect the new signals to forward them
        self.tree_canvas.node_double_clicked.connect(self.node_double_clicked.emit)
        self.tree_canvas.node_right_clicked.connect(self.node_right_clicked.emit)
        
        # Ensure the context menu policy is set
        self.tree_canvas.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        
        # Enable scrollbars
        self.setWidgetResizable(False)  # Important for proper scrolling
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        
        # Set minimum size
        self.setMinimumSize(400, 300)
        
        # Dark mode property
        self._dark_mode = True
        self.tree_canvas.dark_mode = True
        
        # Font property - this will be a pass-through to the tree canvas
        # Make sure we get a valid QFont object
        if hasattr(self.tree_canvas, "sequence_font") and isinstance(self.tree_canvas.sequence_font, QFont):
            self._sequence_font = self.tree_canvas.sequence_font
        else:
            # Create a default font if no valid font is available
            self._sequence_font = QFont("Arial", 9)
        
        # Virtual canvas size (we'll adjust this based on zoom and tree size)
        # To handle very large trees, we now use an even larger virtual canvas
        self.virtual_size = QSize(3000, 8000)  # Dramatically increased height and width for large trees
        self.updateTreeCanvasSize()
        
        # Set up zoom controls - only call this ONCE
        self.setup_zoom_controls()
        
        # Log signal connections for debugging
        logger.debug("PhyloCanvas initialized with tree canvas")
        logger.debug(f"Rename signal connected: {self.rename_node_signal is not None}")
        logger.debug(f"Delete signal connected: {self.delete_node_signal is not None}")
    
    def set_dark_mode(self, dark_mode):
        """Set dark mode for the tree canvas"""
        # Set our own private attribute directly to avoid property recursion
        self._dark_mode = dark_mode
        
        # Try to safely update the tree_canvas
        try:
            if hasattr(self.tree_canvas, '_dark_mode'):
                # Set private attribute directly if available
                self.tree_canvas._dark_mode = dark_mode
                self.tree_canvas.update()
            elif hasattr(self.tree_canvas, 'set_dark_mode'):
                # Use method if available
                self.tree_canvas.set_dark_mode(dark_mode)
            elif hasattr(self.tree_canvas, 'dark_mode'):
                # Last resort, use property (but might cause recursion)
                self.tree_canvas.dark_mode = dark_mode
        except Exception as e:
            logger.error(f"Error setting dark mode on tree_canvas: {e}")
            
        # Force repaint of our widget
        self.update()
        
    def updateTreeCanvasSize(self):
        """Update the tree canvas size based on zoom level"""
        scaled_width = int(self.virtual_size.width() * self.tree_canvas.scale_factor)
        scaled_height = int(self.virtual_size.height() * self.tree_canvas.scale_factor)
        self.tree_canvas.setFixedSize(scaled_width, scaled_height)
        
    def setup_zoom_controls(self):
        """Set up zoom control buttons"""
        # Create a widget to hold the zoom controls
        self.controls_widget = QWidget(self)
        controls_layout = QHBoxLayout(self.controls_widget)
        controls_layout.setContentsMargins(5, 5, 5, 5)
        controls_layout.setSpacing(2)
        
        # PyQt6 doesn't have a direct method to set tooltip delay
        # Instead, we'll style each button to provide helpful tooltips
        
        # Create delete button (proper trashcan icon)
        self.delete_btn = QPushButton("🗑️")  # Emoji with variation selector for better trash can icon
        self.delete_btn.setToolTip("Delete Mode\nClick here then click on sequences to delete them")
        self.delete_btn.setFixedSize(30, 30)
        self.delete_btn.setStyleSheet("font-size: 20px; font-weight: bold; border-radius: 4px;")  # Use default background color
        self.delete_btn.setCheckable(True)  # Make it toggleable
        self.delete_btn.clicked.connect(self.toggle_delete_mode)
        
        # Create reroot button
        self.reroot_btn = QPushButton("⚠")  # Using a root-like symbol
        self.reroot_btn.setToolTip("Reroot Tree\nClick here then click on a branch to make it the new root")
        self.reroot_btn.setFixedSize(30, 30)
        self.reroot_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        self.reroot_btn.setCheckable(True)  # Make it toggleable
        self.reroot_btn.clicked.connect(self.toggle_reroot_mode)
        
        # Create horizontal spacing control buttons with better visual design
        self.h_contract_btn = QPushButton("◀▶")
        self.h_contract_btn.setToolTip("Compress branches horizontally\nDecreases the length of branches in the tree")
        self.h_contract_btn.setFixedSize(36, 30)
        self.h_contract_btn.setStyleSheet("font-size: 12px; font-weight: bold; padding: 0px 2px;")
        self.h_contract_btn.clicked.connect(self.decrease_horizontal_spacing)
        
        self.h_expand_btn = QPushButton("◀  ▶")
        self.h_expand_btn.setToolTip("Expand branches horizontally\nIncreases the length of branches in the tree")
        self.h_expand_btn.setFixedSize(36, 30) 
        self.h_expand_btn.setStyleSheet("font-size: 12px; font-weight: bold; padding: 0px 2px;")
        self.h_expand_btn.clicked.connect(self.increase_horizontal_spacing)
        
        # Create vertical spacing control buttons with better visual design
        self.v_contract_btn = QPushButton("▲▼")
        self.v_contract_btn.setToolTip("Compress tree vertically\nReduces the space between branches")
        self.v_contract_btn.setFixedSize(36, 30)
        self.v_contract_btn.setStyleSheet("font-size: 12px; font-weight: bold; padding: 2px 0px;")
        self.v_contract_btn.clicked.connect(self.decrease_vertical_spacing)
        
        self.v_expand_btn = QPushButton("▲\n▼")
        self.v_expand_btn.setToolTip("Expand tree vertically\nIncreases the space between branches")
        self.v_expand_btn.setFixedSize(36, 30)
        self.v_expand_btn.setStyleSheet("font-size: 12px; font-weight: bold; padding: 2px 0px;")
        self.v_expand_btn.clicked.connect(self.increase_vertical_spacing)
        
        # Create 100% zoom button
        self.zoom_reset_btn = QPushButton("100%")
        self.zoom_reset_btn.setToolTip("Reset zoom level\nResets the view to standard size")
        self.zoom_reset_btn.setFixedSize(40, 30)
        self.zoom_reset_btn.setStyleSheet("font-size: 11px; font-weight: bold;")
        self.zoom_reset_btn.clicked.connect(self.reset_zoom)
        
        # Create fit to screen button (with screen icon)
        self.zoom_fit_btn = QPushButton("⤧")
        self.zoom_fit_btn.setToolTip("Fit tree to window\nResizes the tree to fit the entire view")
        self.zoom_fit_btn.setFixedSize(30, 30)
        self.zoom_fit_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        self.zoom_fit_btn.clicked.connect(self.fit_to_screen)
        
        # Create zoom in button
        self.zoom_in_btn = QPushButton("+")
        self.zoom_in_btn.setToolTip("Zoom in\nMakes the tree larger")
        self.zoom_in_btn.setFixedSize(30, 30)
        self.zoom_in_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        self.zoom_in_btn.clicked.connect(self.zoom_in)
        
        # Create zoom out button
        self.zoom_out_btn = QPushButton("-")
        self.zoom_out_btn.setToolTip("Zoom out\nMakes the tree smaller")
        self.zoom_out_btn.setFixedSize(30, 30)
        self.zoom_out_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        self.zoom_out_btn.clicked.connect(self.zoom_out)
        
        # Add buttons to layout in the requested order (with delete and reroot first)
        controls_layout.addWidget(self.delete_btn)
        controls_layout.addWidget(self.reroot_btn)
        controls_layout.addWidget(self.h_contract_btn)
        controls_layout.addWidget(self.h_expand_btn)
        controls_layout.addWidget(self.v_contract_btn)
        controls_layout.addWidget(self.v_expand_btn)
        controls_layout.addWidget(self.zoom_fit_btn)
        controls_layout.addWidget(self.zoom_reset_btn)
        controls_layout.addWidget(self.zoom_out_btn)
        controls_layout.addWidget(self.zoom_in_btn)
        
        # Position the controls at the bottom right
        self.controls_widget.setFixedSize(controls_layout.sizeHint())
        self.update_controls_position()
        
        # Make sure it stays on top
        self.controls_widget.raise_()
        
    def update_controls_position(self):
        """Update the position of zoom controls when the window is resized"""
        self.controls_widget.move(
            self.width() - self.controls_widget.width() - 10,
            self.height() - self.controls_widget.height() - 10
        )
    
    def resizeEvent(self, event):
        """Handle resize events to reposition controls"""
        super().resizeEvent(event)
        self.update_controls_position()
        
    def get_canvas_width(self):
        """Get the minimum width needed to display the tree"""
        if not self.tree or not self.node_positions:
            return self.width()
            
        # Calculate tree width based on node positions
        margin = 50
        max_label_width = 0
        
        # Calculate max label width to reserve space
        font_metrics = QFontMetrics(self.sequence_font)
        for clade in self.tree.get_terminals():
            label = clade.name if hasattr(clade, 'name') and clade.name else "Unnamed"
            width = font_metrics.horizontalAdvance(label) + 10  # Add padding
            max_label_width = max(max_label_width, width)
            
        # Find rightmost x position
        max_x = 0
        for pos in self.node_positions.values():
            max_x = max(max_x, pos[0])
            
        # Total width = rightmost node + max label width + right margin
        return max_x + max_label_width + margin
        
    def get_canvas_height(self):
        """Get the minimum height needed to display the tree"""
        if not self.tree or not self.node_positions:
            return self.height()
            
        # Find the bottom-most y position
        margin = 50
        max_y = 0
        for pos in self.node_positions.values():
            max_y = max(max_y, pos[1])
            
        # Total height = bottom-most node + bottom margin
        return max_y + margin
        
    def zoom_in(self):
        """Zoom in the tree view"""
        if self.tree_canvas.scale_factor < self.tree_canvas.max_scale:
            self.tree_canvas.scale_factor *= 1.2
            self.updateTreeCanvasSize()
            # Force a layout recalculation
            self.tree_canvas.calculate_tree_layout()
            logger.debug(f"Zoomed in to {self.tree_canvas.scale_factor:.2f}x")
            
    def zoom_out(self):
        """Zoom out the tree view"""
        if self.tree_canvas.scale_factor > self.tree_canvas.min_scale:
            self.tree_canvas.scale_factor *= 0.8
            self.updateTreeCanvasSize()
            # Force a layout recalculation
            self.tree_canvas.calculate_tree_layout()
            logger.debug(f"Zoomed out to {self.tree_canvas.scale_factor:.2f}x")
            
    def reset_zoom(self):
        """Reset zoom to 100%"""
        self.tree_canvas.scale_factor = 1.0
        self.updateTreeCanvasSize()
        
        # Center the view more precisely using float division
        h_center = (self.tree_canvas.width() - self.viewport().width()) / 2.0
        v_center = (self.tree_canvas.height() - self.viewport().height()) / 2.0
        
        self.horizontalScrollBar().setValue(int(h_center))
        self.verticalScrollBar().setValue(int(v_center))
        
        logger.debug("Zoom reset to 1.0x")
        
    def export_image(self, filepath, options=None):
        """Export tree as an image file
        
        Args:
            filepath: Path to save the image
            options: Dictionary with export options:
                - background_color: "current", "black", or "white"
                - high_resolution: Boolean for 2x resolution
                - ultra_resolution: Boolean for 4x resolution
        """
        from PyQt6.QtGui import QPixmap, QPainter
        
        if not self.tree_canvas.tree:
            logger.error("No tree to export")
            return False
            
        if options is None:
            options = {
                "background_color": "current",
                "high_resolution": False,
                "ultra_resolution": False
            }
            
        # Determine resolution scale factor
        scale = 1.0
        if options["high_resolution"]:
            scale = 2.0
        elif options["ultra_resolution"]:
            scale = 4.0
            
        # Save current state
        original_scale = self.tree_canvas.scale_factor
        original_dark_mode = self.tree_canvas.dark_mode
        
        # Set background color according to options
        if options["background_color"] != "current":
            self.tree_canvas.dark_mode = (options["background_color"] == "black")
        
        # Calculate optimal size for the exported image
        self.tree_canvas.calculate_tree_layout()
        tree_width = max(1000, int(self.tree_canvas.width() * scale))
        tree_height = max(800, int(self.tree_canvas.height() * scale))
        
        # Create pixmap at the desired resolution
        pixmap = QPixmap(tree_width, tree_height)
        pixmap.fill(QColor(0, 0, 0, 0))  # Transparent background
        
        # Create painter for the pixmap
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        painter.scale(scale, scale)
        
        # Draw the tree onto the pixmap
        self.tree_canvas.draw_tree(painter)
        painter.end()
        
        # Save the image
        file_format = filepath.split('.')[-1].lower()
        success = pixmap.save(filepath, file_format)
        
        # Restore original settings
        self.tree_canvas.scale_factor = original_scale
        self.tree_canvas.dark_mode = original_dark_mode
        self.tree_canvas.calculate_tree_layout()
        self.tree_canvas.update()
        
        if success:
            logger.info(f"Tree exported to {filepath}")
        else:
            logger.error(f"Failed to export tree to {filepath}")
            
        return success
        
    def fit_to_screen(self):
        """Fit the tree to the available viewport space"""
        if not self.tree_canvas.tree:
            return
        
        # Get the viewport dimensions - this is the actual visible area
        viewport_width = self.viewport().width()
        viewport_height = self.viewport().height()
        
        logger.debug(f"Viewport size: {viewport_width}x{viewport_height}")
        
        # Calculate base tree size by rendering at 1.0 scale first
        self.tree_canvas.scale_factor = 1.0
        self.updateTreeCanvasSize()
        
        # Force layout calculation
        self.tree_canvas.calculate_tree_layout()
        
        # Calculate the actual tree bounds from node positions
        if hasattr(self.tree_canvas, 'node_positions') and self.tree_canvas.node_positions:
            # Find the min/max coordinates of all nodes to get actual tree bounds
            min_x = min(x for x, y in self.tree_canvas.node_positions.values())
            max_x = max(x for x, y in self.tree_canvas.node_positions.values())
            min_y = min(y for x, y in self.tree_canvas.node_positions.values())
            max_y = max(y for x, y in self.tree_canvas.node_positions.values())
            
            # Account for labels on leaves - add extra space
            # Find label widths for leaf nodes
            label_padding = 0
            font_metrics = self.tree_canvas.fontMetrics()
            for clade in self.tree_canvas.tree.get_terminals():
                if clade.name:
                    label_width = font_metrics.horizontalAdvance(clade.name) + 10
                    label_padding = max(label_padding, label_width)
            
            # Calculate tree dimensions with dynamic padding based on tree size
            tree_width = max_x - min_x
            tree_height = max_y - min_y
            
            # Use minimal padding to maximize tree size when fitting to screen
            h_padding = max(20, label_padding, int(tree_width * 0.05))
            v_padding = max(20, int(tree_height * 0.05))
            
            # Adjust bounds with padding
            total_width = tree_width + h_padding
            total_height = tree_height + v_padding
            
            logger.debug(f"Actual tree size: {total_width}x{total_height} with padding")
            
            # Calculate required scale
            width_scale = viewport_width / total_width
            height_scale = viewport_height / total_height
            
            # Use the smaller scale factor to ensure the tree fits
            # Use 98% of available space to get close to the edges as requested
            ideal_scale = min(width_scale, height_scale) * 0.98
            
            # Ensure reasonable bounds and handle very small trees
            ideal_scale = max(0.3, min(self.tree_canvas.max_scale, ideal_scale))
            
            logger.debug(f"Calculated ideal scale: {ideal_scale:.2f}")
            
            # Apply scale
            self.tree_canvas.scale_factor = ideal_scale
        else:
            # Fallback if node positions aren't available
            terminal_count = len(list(self.tree_canvas.tree.get_terminals()))
            
            # Use a simple heuristic based on terminal count
            if terminal_count <= 5:
                # For very small trees, use a larger scale
                self.tree_canvas.scale_factor = 0.5
            elif terminal_count <= 10:
                self.tree_canvas.scale_factor = 0.7
            elif terminal_count <= 50:
                self.tree_canvas.scale_factor = 0.9
            else:
                # For large trees, use a smaller scale
                self.tree_canvas.scale_factor = min(1.0, 50 / terminal_count)
        
        self.updateTreeCanvasSize()
        
        # Force another layout calculation with the new scale
        self.tree_canvas.calculate_tree_layout()
        
        # Center the view more precisely
        hbar = self.horizontalScrollBar()
        vbar = self.verticalScrollBar()
        
        if hbar.maximum() > 0:
            # Center horizontally using float division for more precision
            hbar.setValue(int(hbar.maximum() / 2.0))
        if vbar.maximum() > 0:
            # Center vertically using float division for more precision
            vbar.setValue(int(vbar.maximum() / 2.0))
        
    def handle_wheel_zoom(self, event):
        """Handle mouse wheel events for zooming"""
        # Make sure we have a tree to zoom
        if not self.tree_canvas.tree:
            return
            
        # We shouldn't actually need this method anymore since
        # the TreeCanvas now handles zooming directly
        logger.debug("PhyloCanvas.handle_wheel_zoom called but should be handled by TreeCanvas")
        
        # Pass the event to the tree canvas
        self.tree_canvas.wheelEvent(event)
            
    def set_tree(self, tree):
        """Set the tree to display"""
        try:
            logger.debug(f"Setting tree in PhyloCanvas: {type(tree)}")
            
            if tree is None:
                logger.warning("Cannot set tree - tree is None")
                # Show message in the tree canvas
                self.tree_canvas.tree = None
                # Clear any cached data
                if hasattr(self.tree_canvas, 'node_name_map'):
                    self.tree_canvas.node_name_map = {}
                if hasattr(self.tree_canvas, 'object_id_map'):
                    self.tree_canvas.object_id_map = {}
                self.tree_canvas.update()
                return
                
            # Check if the tree has a root and terminals
            if hasattr(tree, 'root'):
                logger.debug(f"Tree root: {tree.root}")
                terminals = list(tree.get_terminals())
                logger.debug(f"Tree has {len(terminals)} terminal nodes")
                
                # Log some terminal names for debugging
                for i, term in enumerate(terminals[:3]):
                    logger.debug(f"  Terminal {i}: {term.name}")
                    
                if len(terminals) < 2:
                    logger.warning(f"Tree has only {len(terminals)} terminal nodes - may not display properly")
            else:
                logger.warning("Tree has no root attribute - may not be valid")
            
            # IMPORTANT: Make sure the tree is valid before setting it
            # This is a critical step to avoid render issues
            from Bio.Phylo.BaseTree import Tree
            if not isinstance(tree, Tree):
                logger.error(f"Object is not a Bio.Phylo.BaseTree.Tree: {type(tree)}")
                # Attempt to convert or wrap it?
                if hasattr(tree, 'as_phyloxml') and callable(getattr(tree, 'as_phyloxml', None)):
                    logger.debug("Attempting to convert tree to PhyloXML format")
                    try:
                        tree = tree.as_phyloxml()
                    except Exception as e:
                        logger.error(f"Failed to convert tree: {e}")
            
            # Pass the tree to the TreeCanvas
            self.tree_canvas.set_tree(tree)
            
            # Force layouts to be recalculated
            if hasattr(self.tree_canvas, 'calculate_tree_layout'):
                self.tree_canvas.calculate_tree_layout()
            
            # Resize the tree canvas to be sure it's large enough
            self.tree_canvas.setMinimumSize(800, 600)
            self.updateTreeCanvasSize()
            
            # Auto-fit to screen when a new tree is set
            if tree:
                logger.debug("Auto-fitting tree to screen")
                self.fit_to_screen()
                
            # Force a redraw of the entire widget to ensure changes are visible
            self.tree_canvas.update()
            self.viewport().update()
            self.update()
            
            # Add a delayed second update for reliability 
            from PyQt6.QtWidgets import QApplication
            QApplication.processEvents()
            self.tree_canvas.update()
            
        except Exception as e:
            import traceback
            logger.error(f"Error setting tree in PhyloCanvas: {str(e)}")
            logger.error(traceback.format_exc())
            # Try to display the error in the tree canvas
            self.tree_canvas.tree = None
            self.tree_canvas.update()
        
    def highlight_branch(self, clade_name, color=QColor(255, 0, 0)):
        """Highlight a branch with the given color"""
        self.tree_canvas.highlight_branch(clade_name, color)
        
    def add_clade_label(self, clade_name, label):
        """Add a label to a clade"""
        self.tree_canvas.add_clade_label(clade_name, label)
        
    def clear_highlights(self):
        """Clear all branch highlights"""
        self.tree_canvas.clear_highlights()
    
    def increase_vertical_spacing(self):
        """Increase the vertical spacing between branches"""
        result = self.tree_canvas.increase_vertical_spacing()
        # Update canvas size after spacing change
        if result:
            self.updateTreeCanvasSize()
        return result
        
    def decrease_vertical_spacing(self):
        """Decrease the vertical spacing between branches"""
        result = self.tree_canvas.decrease_vertical_spacing()
        # Update canvas size after spacing change
        if result:
            self.updateTreeCanvasSize()
        return result
        
    def increase_horizontal_spacing(self):
        """Increase the horizontal branch lengths"""
        result = self.tree_canvas.increase_horizontal_spacing()
        # Update canvas size after spacing change
        if result:
            self.updateTreeCanvasSize()
        return result
        
    def decrease_horizontal_spacing(self):
        """Decrease the horizontal branch lengths"""
        result = self.tree_canvas.decrease_horizontal_spacing()
        # Update canvas size after spacing change
        if result:
            self.updateTreeCanvasSize()
        return result
        
    def enter_reroot_mode(self):
        """Enter reroot mode to allow the user to reroot the tree"""
        return self.tree_canvas.enter_reroot_mode()
        
    def exit_reroot_mode(self):
        """Exit reroot mode"""
        result = self.tree_canvas.exit_reroot_mode()
        # Update button state
        self.reroot_btn.setChecked(False)
        return result
        
    def enter_delete_mode(self):
        """Enter delete mode to allow the user to delete sequences"""
        import logging
        logger = logging.getLogger("treecraft")
        
        if not self.tree_canvas.tree:
            logger.warning("Cannot enter delete mode - no tree loaded")
            return False
        
        # Set the delete mode flag
        self.delete_mode = True
        logger.debug("PhyloCanvas delete_mode set to True")
        
        # Create and set a trash can cursor for delete mode
        from PyQt6.QtGui import QCursor, QPixmap
        
        # First try to create a custom trash can cursor
        try:
            # Create a small pixmap for the cursor
            cursor_size = 32
            pixmap = QPixmap(cursor_size, cursor_size)
            pixmap.fill(Qt.GlobalColor.transparent)
            
            # Draw the trash can on the pixmap
            from PyQt6.QtGui import QPainter, QColor, QPen, QFont
            painter = QPainter(pixmap)
            painter.setFont(QFont("Arial", 20))
            painter.setPen(QPen(QColor(255, 0, 0)))
            painter.drawText(0, 25, "🗑️")  # Draw trash can emoji
            painter.end()
            
            # Create and set the cursor
            cursor = QCursor(pixmap, 8, 8)  # Set hotspot near top-left
            self.tree_canvas.setCursor(cursor)
            logger.debug("Set custom trash can cursor")
        except Exception as e:
            # Fall back to cross cursor if custom cursor fails
            logger.debug(f"Failed to create custom cursor: {e}, using crosshair instead")
            self.tree_canvas.setCursor(Qt.CursorShape.CrossCursor)
        
        # Update status message to inform user
        self.update_status.emit("Delete Mode: Click on sequences to delete them")
        
        # Connect mouse events in tree canvas if needed
        self.tree_canvas.delete_mode = True
        logger.debug("TreeCanvas delete_mode set to True")
        
        # Force the view to update
        self.tree_canvas.update()
        self.update()
        
        # Print a very explicit message to confirm we entered delete mode
        logger.info("🔥 DELETE MODE ACTIVATED 🔥")
        
        return True
        
    def exit_delete_mode(self):
        """Exit delete mode"""
        import logging
        logger = logging.getLogger("treecraft")
        
        # Reset the delete mode flag
        self.delete_mode = False
        logger.debug("PhyloCanvas delete_mode set to False")
        
        # Reset cursor to default
        self.tree_canvas.setCursor(Qt.CursorShape.ArrowCursor)
        
        # Update status message
        self.update_status.emit("Delete mode disabled")
        
        # Disconnect mouse events in tree canvas if needed
        self.tree_canvas.delete_mode = False
        
        # Clear any hovered sequence highlighting
        if hasattr(self.tree_canvas, 'hovered_sequence'):
            self.tree_canvas.hovered_sequence = None
        
        logger.debug("TreeCanvas delete_mode set to False")
        
        # Update button state
        self.delete_btn.setChecked(False)
        
        # Force the view to update
        self.tree_canvas.update()
        self.update()
        
        # Print a very explicit message to confirm we exited delete mode
        logger.info("❌ DELETE MODE DEACTIVATED ❌")
        
        return True
        
    def set_branch_width(self, width):
        """Set the branch line width
        
        Args:
            width: Width in pixels (1-5)
        """
        # Make sure width is a valid integer
        width = max(1, min(int(width), 5))
        
        # Update tree canvas
        if hasattr(self.tree_canvas, 'set_branch_width'):
            self.tree_canvas.set_branch_width(width)
        elif hasattr(self.tree_canvas, 'branch_width'):
            self.tree_canvas.branch_width = width
            self.tree_canvas.update()
        
    def direct_delete_node(self, node_name):
        """Direct method to delete a node from the tree
        
        This method will try to find the MainWindow and call its delete_tree_node method directly.
        It's used as a fallback when the signal chain isn't working.
        """
        import logging
        logger = logging.getLogger("treecraft")
        logger.debug(f"PhyloCanvas direct_delete_node called for: {node_name}")
        
        try:
            # Try to find the main window by traversing the parent hierarchy
            window = self
            while window and not hasattr(window, 'delete_tree_node'):
                window = window.parent()
                
            if window and hasattr(window, 'delete_tree_node'):
                logger.debug(f"Found MainWindow, calling delete_tree_node directly")
                # Call the MainWindow method directly
                window.delete_tree_node(node_name)
            else:
                # If we can't find the MainWindow, log an error
                logger.error(f"Could not find MainWindow with delete_tree_node method")
        except Exception as e:
            # Catch any exceptions to prevent crashes
            logger.error(f"Error in direct_delete_node: {str(e)}")
            
            # Force update of the tree view to ensure it refreshes
            self.update()
    
    def toggle_delete_mode(self, checked):
        """Toggle delete mode on/off based on button state"""
        if not self.tree_canvas.tree:
            # Can't enter delete mode without a tree
            self.delete_btn.setChecked(False)
            return
            
        # Make sure we don't have both delete mode and reroot mode enabled
        if checked and self.reroot_btn.isChecked():
            self.reroot_btn.setChecked(False)
            self.exit_reroot_mode()
            
        if checked:
            # Entering delete mode
            success = self.enter_delete_mode()
            if not success:
                # Reset button if entering delete mode failed
                self.delete_btn.setChecked(False)
        else:
            # Exiting delete mode
            self.exit_delete_mode()
    
    def toggle_reroot_mode(self, checked):
        """Toggle reroot mode on/off based on button state"""
        if not self.tree_canvas.tree:
            # Can't enter reroot mode without a tree
            self.reroot_btn.setChecked(False)
            return
            
        # Make sure we don't have both delete mode and reroot mode enabled
        if checked and self.delete_btn.isChecked():
            self.delete_btn.setChecked(False)
            self.exit_delete_mode()
            
        if checked:
            # Entering reroot mode
            success = self.enter_reroot_mode()
            if not success:
                # Reset button if entering reroot mode failed
                self.reroot_btn.setChecked(False)
        else:
            # Exiting reroot mode
            self.exit_reroot_mode()
        
    @property
    def dark_mode(self):
        return self._dark_mode
        
    @dark_mode.setter
    def dark_mode(self, value):
        # Guard against self setting and store old value
        old_value = self._dark_mode
        self._dark_mode = value
        
        # Only update tree_canvas if value changed and tree_canvas exists
        if old_value != value and hasattr(self, 'tree_canvas'):
            try:
                # Access private attribute directly to avoid recursion
                if hasattr(self.tree_canvas, '_dark_mode'):
                    self.tree_canvas._dark_mode = value
                # Fall back to set_dark_mode method if available
                elif hasattr(self.tree_canvas, 'set_dark_mode'):
                    self.tree_canvas.set_dark_mode(value)
                # Last resort - try the property
                else:
                    self.tree_canvas.dark_mode = value
                
                # Force update
                self.tree_canvas.update()
            except Exception as e:
                logger.error(f"Error setting dark mode on tree_canvas: {e}")
                
        # Update self
        self.update()
        
    @property
    def sequence_font(self):
        return self._sequence_font
        
    @sequence_font.setter
    def sequence_font(self, font):
        # Log for debugging
        logger.debug(f"Setting sequence_font in PhyloCanvas to: {font}")
        self._sequence_font = font
        
        # Only update tree_canvas if it exists
        if hasattr(self, 'tree_canvas'):
            # Direct attribute setting to prevent potential recursion
            if hasattr(self.tree_canvas, '_sequence_font'):
                self.tree_canvas._sequence_font = font
            else:
                self.tree_canvas.sequence_font = font
            
            # Force a recalculation of the tree layout
            if hasattr(self.tree_canvas, 'calculate_tree_layout'):
                self.tree_canvas.calculate_tree_layout()
                
            # Force a complete repaint
            self.tree_canvas.update()
        
        # Update the scroll area
        if hasattr(self, 'viewport'):
            self.viewport().update()
        self.update()