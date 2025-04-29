from PyQt6.QtWidgets import QListWidget, QListWidgetItem, QMenu, QInputDialog, QMessageBox
from PyQt6.QtCore import Qt, pyqtSignal

class SequenceListWidget(QListWidget):
    """Widget to display and manage sequence list"""
    
    # Define signals
    rename_signal = pyqtSignal(str)  # Signal when renaming a sequence (passes seq_id)
    delete_signal = pyqtSignal(str)  # Signal when deleting a sequence (passes seq_id)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequences = {}  # Dictionary of sequence_id -> sequence
        self.descriptions = {}  # Dictionary of sequence_id -> description
        
        # Set horizontal scroll mode to ScrollPerPixel for smoother scrolling
        self.setHorizontalScrollMode(QListWidget.ScrollMode.ScrollPerPixel)
        
        # Enable horizontal scrollbar if needed
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        
        # Enable multi-selection mode
        self.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        
        # Enable context menu
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)
        
        # Connect double-click signal
        self.itemDoubleClicked.connect(self.on_item_double_clicked)
        
    def show_context_menu(self, position):
        """Show a context menu for the selected item"""
        # Only show if an item is selected
        if not self.selectedItems():
            return
            
        # Create context menu
        context_menu = QMenu(self)
        
        # Add actions
        rename_action = context_menu.addAction("Rename Sequence")
        delete_action = context_menu.addAction("Delete Sequence")
        
        # Show the menu and get the selected action
        action = context_menu.exec(self.mapToGlobal(position))
        
        # Handle the selected action
        if action == rename_action:
            self.handle_rename()
        elif action == delete_action:
            self.handle_delete()
            
    def on_item_double_clicked(self, item):
        """Handle double click on an item"""
        self.handle_rename()
        
    def handle_rename(self):
        """Rename the selected sequence"""
        # Get the selected item
        if not self.selectedItems() or len(self.selectedItems()) != 1:
            return
            
        import logging
        logger = logging.getLogger("treecraft.sequence_list")
            
        item = self.selectedItems()[0]
        item_idx = self.row(item)
        description = item.text()
        
        # Try to get sequence ID directly from item data
        seq_id = item.data(Qt.ItemDataRole.UserRole)
        logger.debug(f"Item UserRole data: {seq_id}")
        
        # If we have a valid sequence ID from item data
        if seq_id and seq_id in self.sequences:
            logger.debug(f"Using sequence ID from item data: {seq_id}")
            if self.rename_signal:
                self.rename_signal.emit(seq_id)
            return
        
        # Find the sequence ID for this item - need to handle duplicate descriptions
        sequence_id = None
        
        # First, try to match using both description and list index
        matching_ids = []
        for id, desc in self.descriptions.items():
            if desc == description:
                matching_ids.append(id)
        
        # If we have multiple matches, try to determine which one by position
        if len(matching_ids) > 1:
            logger.debug(f"Multiple sequence IDs match description '{description}': {matching_ids}")
            logger.debug(f"Selected item index: {item_idx}")
            
            # The index in the list should correspond to the position in matching_ids
            if item_idx < len(matching_ids):
                sequence_id = matching_ids[item_idx]
                logger.debug(f"Selected sequence ID by position: {sequence_id}")
            else:
                # Fallback - just use the first one
                sequence_id = matching_ids[0]
                logger.debug(f"Fallback to first matching ID: {sequence_id}")
        elif matching_ids:
            # Just one match
            sequence_id = matching_ids[0]
        
        if sequence_id and self.rename_signal:
            logger.debug(f"Emitting rename signal for sequence ID: {sequence_id}")
            self.rename_signal.emit(sequence_id)
            
    def handle_delete(self):
        """Delete the selected sequence"""
        # Get the selected item
        if not self.selectedItems() or len(self.selectedItems()) != 1:
            return
            
        import logging
        logger = logging.getLogger("treecraft.sequence_list")
        
        item = self.selectedItems()[0]
        description = item.text()
        
        # Try to get sequence ID directly from item data
        seq_id = item.data(Qt.ItemDataRole.UserRole)
        logger.debug(f"Item UserRole data for delete: {seq_id}")
        
        # If we have a valid sequence ID from item data
        if seq_id and seq_id in self.sequences:
            logger.debug(f"Using sequence ID from item data for delete: {seq_id}")
            if self.delete_signal:
                self.delete_signal.emit(seq_id)
            return
        
        # Find the sequence ID for this item
        sequence_id = None
        for id, desc in self.descriptions.items():
            if desc == description:
                sequence_id = id
                break
                
        if sequence_id and self.delete_signal:
            logger.debug(f"Emitting delete signal for sequence ID: {sequence_id}")
            self.delete_signal.emit(sequence_id)
    def rebuild_sequence_list(self):
        """Completely rebuild the list display from internal data"""
        # Remember the current selection if possible
        current_row = self.currentRow()
        current_text = self.currentItem().text() if self.currentItem() else None
        
        # Clear the list widget
        self.clear()
        
        # Rebuild from our internal data
        for seq_id, desc in self.descriptions.items():
            item = QListWidgetItem(desc)
            item.setData(Qt.ItemDataRole.UserRole, seq_id)  # Store sequence ID in item data
            self.addItem(item)
        
        # Try to restore selection
        if current_text:
            # Find the item with the same text
            for i in range(self.count()):
                if self.item(i).text() == current_text:
                    self.setCurrentRow(i)
                    break
        elif current_row >= 0 and current_row < self.count():
            # Fall back to row number if text doesn't match
            self.setCurrentRow(current_row)

    def update_sequence(self, old_id, new_id, new_sequence, new_description):
        """Update a sequence with new data and refresh the display"""
        # Find the item in the list
        item_to_update = None
        item_index = -1
        
        # Look for the item with the old description
        old_description = self.descriptions.get(old_id, old_id)
        for i in range(self.count()):
            if self.item(i).text() == old_description:
                item_to_update = self.item(i)
                item_index = i
                break
        
        if item_to_update is None:
            print(f"Warning: Could not find item with description '{old_description}' to update")
            return False
        
        # Update the internal data
        if old_id != new_id:
            # ID has changed
            self.sequences[new_id] = new_sequence
            self.descriptions[new_id] = new_description
            
            # Delete old entries
            if old_id in self.sequences:
                del self.sequences[old_id]
            if old_id in self.descriptions:
                del self.descriptions[old_id]
        else:
            # Just update the existing entries
            self.sequences[old_id] = new_sequence
            self.descriptions[old_id] = new_description
        
        # Force update the display
        # The key is to remove and re-add the item
        self.takeItem(item_index)
        new_item = QListWidgetItem(new_description)
        new_item.setData(Qt.ItemDataRole.UserRole, new_id)  # Store sequence ID in item data
        self.insertItem(item_index, new_item)
        self.setCurrentRow(item_index)
        
        return True


    def refresh(self):
        """Refresh the list display from the internal data"""
        # Remember the current selection if any
        selected_items = self.selectedItems()
        selected_texts = [item.text() for item in selected_items]
        
        # Clear and rebuild the list
        self.clear()
        
        # Add all sequences back to the list
        for seq_id, description in self.descriptions.items():
            item = QListWidgetItem(description)
            item.setData(Qt.ItemDataRole.UserRole, seq_id)  # Store sequence ID in item data
            self.addItem(item)
        
        # Restore selection if possible
        for i in range(self.count()):
            item = self.item(i)
            if item and item.text() in selected_texts:
                item.setSelected(True)

    def add_sequence(self, sequence_id, sequence, description=None):
        """Add a sequence to the list"""
        self.sequences[sequence_id] = sequence
        
        # Store the description - fix duplicate display names but preserve full content
        import logging
        logger = logging.getLogger("treecraft.sequence_list")
        
        if description is None:
            # If no description provided, use the sequence ID
            description = sequence_id
            logger.debug(f"No description provided, using ID as description: {description}")
        elif description.startswith(sequence_id) and description.count(sequence_id) > 1:
            # Check for duplicate IDs in the description (e.g., "HQ604103.1_ HQ604103.1_")
            # This can happen when loading files with repeated IDs in the description
            logger.debug(f"Found duplicate ID in description: '{description}', fixing...")
            
            # Remove duplicated IDs from description
            # Split by spaces and filter out duplicates while preserving order
            parts = description.split()
            seen = set()
            unique_parts = []
            for part in parts:
                if part not in seen:
                    seen.add(part)
                    unique_parts.append(part)
            description = " ".join(unique_parts)
            logger.debug(f"Fixed description: '{description}'")
            
        # Ensure description is never shorter than sequence_id
        if len(description) < len(sequence_id):
            logger.debug(f"Description shorter than ID, using ID instead: {sequence_id}")
            description = sequence_id
            
        # Make sure we're using the full description
        logger.debug(f"Final description for {sequence_id}: '{description}'")
        self.descriptions[sequence_id] = description
        
        # Add to the list with the fixed description
        item = QListWidgetItem(description)
        item.setData(Qt.ItemDataRole.UserRole, sequence_id)  # Store sequence ID in item data
        self.addItem(item)
    
    def remove_selected_sequences(self):
        """Remove selected sequences from the list"""
        for item in self.selectedItems():
            # Find the sequence ID for this item
            description = item.text()
            sequence_id = None
            
            # Find the ID that matches this description
            for id, desc in self.descriptions.items():
                if desc == description:
                    sequence_id = id
                    break
            
            if sequence_id:
                del self.sequences[sequence_id]
                del self.descriptions[sequence_id]
            
            self.takeItem(self.row(item))
    
    def clear_sequences(self):
        """Clear all sequences from the list"""
        self.sequences = {}
        self.descriptions = {}
        self.clear()
        
    def remove_sequence(self, sequence_id):
        """Remove a specific sequence by ID"""
        if sequence_id in self.sequences:
            # Find the item in the list widget
            description = self.descriptions.get(sequence_id, sequence_id)
            for i in range(self.count()):
                if self.item(i).text() == description:
                    self.takeItem(i)
                    break
                    
            # Remove from our dictionaries
            del self.sequences[sequence_id]
            del self.descriptions[sequence_id]
            return True
        return False
    
    def get_sequences(self):
        """Return all sequences as a dictionary"""
        return self.sequences
    
    def get_selected_ids(self):
        """Return IDs of selected sequences"""
        selected_ids = []
        for item in self.selectedItems():
            description = item.text()
            # Find the ID that matches this description
            for id, desc in self.descriptions.items():
                if desc == description:
                    selected_ids.append(id)
                    break
        return selected_ids
