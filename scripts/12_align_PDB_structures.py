import os
import pymol
from pymol import cmd


def load_pymol_session_align(session_input, session_output, structures):
    """
    Loads a PyMOL session, adds all specified structures, and saves it to a new location.

    :param session_input: Path to the .pse.gz session file.
    :param session_output: Path to save the modified PyMOL session.
    :param structures: List of full-path structure filenames to be loaded.
    """

    # Reinitialize PyMOL to clear previous states
    cmd.reinitialize()

    # Load pymol session
    cmd.load(session_input)

    # Get reference structure to align with
    reference = [i for i in cmd.get_names() if 'alphafold2' in i and 'pocket' not in i]
    if len(reference) != 1:
        print("Something weird is going on. Please check.")
    else:
        reference = reference[0]

    # Load each structure file
    for structure in structures:

        # Load structure
        cmd.load(structure)

        # Align with reference
        structure_name = os.path.basename(structure).replace('.ent', '').replace('.pdb', '')
        cmd.align(structure_name, reference)

    # Hide all surface representations
    cmd.hide("surface", "all")

    # Center
    cmd.zoom(reference)

    # Save the session in the new directory
    cmd.save(session_output)


# Define file paths 
root = os.path.dirname(os.path.abspath(__file__))
pymol_sessions = os.path.join(root, "..", "processed", "pymol_sessions")
pymol_sessions_pdbe = os.path.join(root, "..", "processed", "pymol_sessions_pdbe")
pdb_structures = os.path.join(root, "..", "data", "structures", "pdbe_database")
os.makedirs(pymol_sessions_pdbe, exist_ok=True)

# For all proteins having > 0 PDBe structures (both apo and holo)
for uni in sorted(os.listdir(pdb_structures)):

    # Enumerate PDB structures
    structures = sorted(os.listdir(os.path.join(pdb_structures, uni, uni + "_archive-PDB")))
    structures = [i for i in structures if '.ent' in i and '.txt' not in i]
    structures = [os.path.join(pdb_structures, uni, uni + "_archive-PDB", i) for i in structures]

    load_pymol_session_align(os.path.join(pymol_sessions, uni + ".pse.gz"), 
                       os.path.join(pymol_sessions_pdbe, uni + ".pse.gz"),
                       structures)
    

    print(f"Created PyMol session for {uni}")
