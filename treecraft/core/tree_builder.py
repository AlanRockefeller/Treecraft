from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import tempfile
import shutil
import os
import logging
from PyQt6.QtWidgets import QMessageBox
from treecraft.utils.external_tools import find_external_tool, run_raxml, run_mrbayes, find_tree_file

logger = logging.getLogger("treecraft.tree_builder")


def build_tree(alignment):
    """Build a phylogenetic tree using UPGMA method"""
    if len(alignment) < 2:
        QMessageBox.warning(None, "Warning", "Need at least 2 sequences to build a tree")
        return None

    try:
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)

        # Build tree
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        return tree
    except Exception as e:
        QMessageBox.critical(None, "Error", f"Failed to build tree: {str(e)}")
        return None



def build_tree_with_params(alignment, params, progress_callback=None):
    """Build a phylogenetic tree with bootstrap support if requested"""
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio import AlignIO, Phylo
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from io import StringIO
    import tempfile
    import os
    import random
    import logging
    
    logger = logging.getLogger("treecraft.tree_builder")
    
    try:
        # Check if sequences are aligned (all same length)
        seq_lengths = set(len(str(record.seq)) for record in alignment)
        if len(seq_lengths) != 1:
            raise ValueError("Sequences must be aligned (all same length) before building a tree")
        
        # Current parameters
        method = params.get("method", "upgma")
        distance_model = params.get("distance_model", "identity")
        bootstrap = params.get("bootstrap", False)
        bootstrap_replicates = params.get("bootstrap_replicates", 100)
        
        # Create map of record IDs to full descriptions - crucial for preserving spaces
        id_to_desc = {}
        for record in alignment:
            full_desc = record.description if record.description else record.id
            id_to_desc[record.id] = full_desc
        
        # Create a temporary file for the alignment
        fd, temp_path = tempfile.mkstemp(suffix=".fasta")
        os.close(fd)
        
        try:
            # Create a modified alignment without spaces in IDs
            # This is to work around Biopython's limitations with spaces in tree nodes
            safe_alignment = []
            id_map = {}  # Maps safe IDs to original IDs
            
            for i, record in enumerate(alignment):
                # Create a safe ID without spaces for tree building
                safe_id = f"seq_{i}"
                
                # Map the safe ID to the original record's full description
                id_map[safe_id] = record.description if record.description else record.id
                
                # Create new record with safe ID
                safe_record = SeqRecord(
                    record.seq,
                    id=safe_id,
                    description=""
                )
                safe_alignment.append(safe_record)
                
            # Write modified alignment to temp file
            with open(temp_path, "w") as f:
                for record in safe_alignment:
                    f.write(f">{record.id}\n{record.seq}\n")
            
            # Read as alignment
            align = AlignIO.read(temp_path, "fasta")
            
            # Build the main tree
            calculator = DistanceCalculator(distance_model.lower())
            constructor = DistanceTreeConstructor()
            
            # Calculate distance matrix
            dm = calculator.get_distance(align)
            
            # Build tree based on method
            if method == "upgma":
                tree = constructor.upgma(dm)
            elif method == "nj":
                tree = constructor.nj(dm)
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            # Replace safe IDs with full descriptions in the tree
            for terminal in tree.get_terminals():
                if terminal.name in id_map:
                    terminal.name = id_map[terminal.name]
                    logger.debug(f"Restored full name: {terminal.name}")
            
            # Perform bootstrap analysis if requested
            if bootstrap and bootstrap_replicates > 0:
                # Manually implement bootstrapping
                from Bio.Phylo.Consensus import _count_clades, get_support
                
                # Generate bootstrap trees
                bootstrap_trees = []
                
                # Set sequence length
                seq_len = len(align[0])
                
                for i in range(bootstrap_replicates):
                    if progress_callback:
                        progress_callback(int(i / bootstrap_replicates * 100))
                    
                    # Generate bootstrap alignment
                    # Sample columns with replacement
                    indices = [random.randrange(seq_len) for _ in range(seq_len)]
                    
                    # Create bootstrapped alignment
                    boot_align = align[:, 0:0]  # Empty alignment with same structure
                    for idx in indices:
                        boot_align = boot_align + align[:, idx:idx+1]
                    
                    try:
                        # Calculate distance matrix
                        boot_dm = calculator.get_distance(boot_align)
                        
                        # Build bootstrap tree
                        if method == "upgma":
                            boot_tree = constructor.upgma(boot_dm)
                        else:
                            boot_tree = constructor.nj(boot_dm)
                        
                        # Replace IDs with full descriptions in bootstrap tree
                        for terminal in boot_tree.get_terminals():
                            if terminal.name in id_map:
                                terminal.name = id_map[terminal.name]
                                
                        bootstrap_trees.append(boot_tree)
                    except Exception as e:
                        print(f"Error in bootstrap replicate {i}: {e}")
                        continue
                
                # Calculate branch support
                tree = get_support(tree, bootstrap_trees)
                
                if progress_callback:
                    progress_callback(100)
            
            return tree
            
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.remove(temp_path)
                
    except Exception as e:
        print(f"Error building tree: {e}")
        return None


def ensure_aligned(records):
    """Ensure sequences are aligned using a simple method"""
    from Bio import pairwise2
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # For very simple cases, just check if sequences are the same length
    lengths = [len(record.seq) for record in records]
    if len(set(lengths)) == 1:
        # All sequences are the same length, might be aligned already
        return records
    
    # For more complex cases, we need to align
    try:
        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        
        # Use the first sequence as reference
        ref_seq = records[0].seq
        
        # Align each sequence to reference
        aligned_records = [records[0]]  # Start with reference
        
        for record in records[1:]:
            alignment = aligner.align(ref_seq, record.seq)[0]
            aligned_seq = ""
            ref_idx = 0
            seq_idx = 0
            
            for op in alignment.path:
                if op[0] == 0 and op[1] == 0:  # Match or mismatch
                    aligned_seq += record.seq[seq_idx]
                    ref_idx += 1
                    seq_idx += 1
                elif op[0] == 1 and op[1] == 0:  # Gap in reference
                    aligned_seq += record.seq[seq_idx]
                    seq_idx += 1
                elif op[0] == 0 and op[1] == 1:  # Gap in sequence
                    aligned_seq += "-"
                    ref_idx += 1
            
            # Create new SeqRecord with aligned sequence
            aligned_record = SeqRecord(
                Seq(aligned_seq),
                id=record.id,
                description=record.description
            )
            aligned_records.append(aligned_record)
        
        return aligned_records
    except Exception as e:
        print(f"Failed to align sequences: {str(e)}")
        return records  # Return original records if alignment fails


def build_upgma_tree(alignment, distance_model="identity", bootstrap=False, bootstrap_replicates=100, progress_callback=None):
    """Build a phylogenetic tree using UPGMA method"""
    try:
        # Calculate distance matrix
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        
        calculator = DistanceCalculator(distance_model)
        dm = calculator.get_distance(alignment)
        
        # Build tree
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)
        
        return tree
    except Exception as e:
        print(f"Failed to build UPGMA tree: {str(e)}")
        return None

def build_nj_tree(alignment, distance_model="identity", bootstrap=False, bootstrap_replicates=100, progress_callback=None):
    """Build a phylogenetic tree using Neighbor Joining algorithm"""
    try:
        # Calculate distance matrix
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        
        logger.info(f"Building NJ tree for {len(alignment)} sequences with {distance_model} model")
        logger.info(f"First sequence: {alignment[0].id} - Length: {len(alignment[0].seq)}")
        
        calculator = DistanceCalculator(distance_model)
        dm = calculator.get_distance(alignment)
        
        logger.info(f"Distance matrix calculated - size: {len(dm)}")
        
        # Build NJ tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        
        logger.info(f"NJ tree built successfully - tree type: {type(tree)}")
        if hasattr(tree, 'root'):
            logger.info(f"Tree has root: {tree.root}")
            logger.info(f"Terminal nodes: {len(list(tree.get_terminals()))}")
        else:
            logger.warning("Tree has no root attribute!")
            
        return tree
    except Exception as e:
        import traceback
        logger.error(f"Failed to build NJ tree: {str(e)}")
        traceback.print_exc()
        return None


def build_ml_tree(alignment, progress=None):
    """Build a phylogenetic tree using RAxML or MrBayes"""
    from Bio import Phylo, SeqIO
    from Bio.SeqRecord import SeqRecord
    import logging
    import re
    
    logger = logging.getLogger("treecraft.tree_builder")
    
    if len(alignment) < 4:
        QMessageBox.warning(None, "Warning", "Need at least 4 sequences for ML tree construction")
        return None
    
    # Check if RAxML or MrBayes is installed
    raxml_path = find_external_tool("raxml")
    mrbayes_path = find_external_tool("mb")
    
    if not raxml_path and not mrbayes_path:
        reply = QMessageBox.question(
            None, 
            "External Tool Required",
            "RAxML or MrBayes is required for ML tree construction but not found.\n\n"
            "Would you like to download and install one of these tools?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            # Provide instructions for installation
            QMessageBox.information(
                None,
                "Installation Instructions",
                "Please install one of the following:\n\n"
                "1. RAxML: https://github.com/stamatak/standard-RAxML\n"
                "2. MrBayes: https://nbisweden.github.io/MrBayes/\n\n"
                "After installation, ensure the executable is in your PATH."
            )
        return None
    
    # Use safe sequence IDs for tree building and map them back to full names
    safe_alignment = []
    id_map = {}  # Maps safe IDs to original full descriptions
    parentheses_map = {}  # Maps escaped parentheses back to original ones
    
    for i, record in enumerate(alignment):
        # Create a safe ID without spaces for tree building
        safe_id = f"seq_{i}"
        
        # Get the original description
        original_description = record.description if record.description else record.id
        
        # Map the safe ID to the original record's full description
        id_map[safe_id] = original_description
        
        # If the description contains parentheses, store a mapping with escaped parentheses
        if '(' in original_description or ')' in original_description:
            # Replace parentheses with escaped versions
            escaped_description = original_description.replace('(', '\\(').replace(')', '\\)')
            # Add mapping to restore the original
            parentheses_map[escaped_description] = original_description
            # Update id_map to store the escaped version
            id_map[safe_id] = escaped_description
            logger.debug(f"Escaped parentheses in name: '{original_description}' -> '{escaped_description}'")
        
        # Create new record with safe ID
        safe_record = SeqRecord(
            record.seq,
            id=safe_id,
            description=""
        )
        safe_alignment.append(safe_record)
    
    # Save alignment to temporary file
    temp_dir = tempfile.mkdtemp()
    temp_fasta = os.path.join(temp_dir, "alignment.fasta")
    
    try:
        # Write sequences with safe IDs
        with open(temp_fasta, "w") as f:
            for record in safe_alignment:
                f.write(f">{record.id}\n{record.seq}\n")
        
        # Run external tool (RAxML or MrBayes)
        if raxml_path:
            run_raxml(temp_dir, temp_fasta, progress)
        else:
            run_mrbayes(temp_dir, temp_fasta, progress)
            
        # Load resulting tree
        tree_file = find_tree_file(temp_dir)
        if tree_file:
            # Read tree with safe IDs
            tree = Phylo.read(tree_file, "newick")
            
            # Replace safe IDs with full descriptions in the tree
            for terminal in tree.get_terminals():
                if terminal.name in id_map:
                    original_name = id_map[terminal.name]
                    logger.debug(f"Restoring name from '{terminal.name}' to '{original_name}'")
                    terminal.name = original_name
                    logger.debug(f"Restored full name in ML tree: {terminal.name}")
                else:
                    # For RAxML trees, they may include location suffixes with underscores
                    # Try to find a match by checking if the name starts with a known prefix
                    found = False
                    if "_" in terminal.name:
                        base_name = terminal.name.split("_")[0]
                        for safe_id, full_name in id_map.items():
                            # Try different matching strategies to find the full name
                            if (safe_id.startswith(base_name) or 
                                terminal.name.startswith(safe_id) or
                                base_name in full_name or
                                any(part in full_name for part in terminal.name.split("_"))):
                                logger.info(f"Found match: {terminal.name} -> {full_name}")
                                terminal.name = full_name
                                found = True
                                break
                    
                    # Extended search if not found with prefix
                    if not found:
                        # Look for fragments of the terminal name in full names
                        for safe_id, full_name in id_map.items():
                            # Check if any word in terminal name appears in the full name
                            if terminal.name in full_name or any(part in full_name for part in terminal.name.split("_")):
                                logger.info(f"Found partial match by content: {terminal.name} -> {full_name}")
                                terminal.name = full_name
                                found = True
                                break
                    
                    if not found:
                        logger.warning(f"Could not find mapping for terminal: {terminal.name}")
            
            # Now restore any escaped parentheses to their original form
            for terminal in tree.get_terminals():
                if terminal.name in parentheses_map:
                    original_with_parentheses = parentheses_map[terminal.name]
                    logger.debug(f"Restoring parentheses: '{terminal.name}' -> '{original_with_parentheses}'")
                    terminal.name = original_with_parentheses
                # Handle case where RAxML kept the escaped parentheses (backslashes may be present)
                elif '\\(' in terminal.name or '\\)' in terminal.name:
                    # Create a version with normal parentheses for comparison
                    unescaped = terminal.name.replace('\\(', '(').replace('\\)', ')')
                    if unescaped in parentheses_map.values():
                        logger.debug(f"Restoring from escaped sequence: '{terminal.name}' -> '{unescaped}'")
                        terminal.name = unescaped
            
            return tree
        else:
            QMessageBox.warning(None, "Warning", "No tree file generated")
            return None
            
    except Exception as e:
        QMessageBox.critical(None, "Error", f"Failed to build ML tree: {str(e)}")
        return None
    finally:
        # Clean up
        try:
            shutil.rmtree(temp_dir)
        except:
            pass



