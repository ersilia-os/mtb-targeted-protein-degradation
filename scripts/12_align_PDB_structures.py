import os
import pymol
from pymol import cmd
import numpy as np
import pandas as pd


def load_pymol_session_align(session_input, session_output, structures):
    """
    Loads a PyMOL session, adds all specified structures, and saves it to a new location.

    :param session_input: Path to the .pse.gz session file.
    :param session_output: Path to save the modified PyMOL session.
    :param structures: List of full-path structure filenames to be loaded.
    """

    # Storing PDBe annotations
    report = []
    uni = session_input.split("/")[-1].split(".pse.gz")[0]

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

    # Get all pocket coordinates
    pocket_to_coords = {}
    pockets = sorted([i for i in cmd.get_names("all") if 'pocket' in i])
    for pocket in pockets:
        model = cmd.get_model(pocket)
        if len(model.atom) != 1:
            print(f"Warning: {obj} has {len(model.atom)} atoms (expected 1)")
        atom = model.atom[0]
        pocket_to_coords[pocket] = [atom.coord[0], atom.coord[1], atom.coord[2]] 

    # Load each structure file
    for structure in structures:

        # Load structure
        cmd.load(structure)

        # Align with reference
        structure_name = os.path.basename(structure).replace('.ent', '').replace('.pdb', '').replace('.cif', '')
        cmd.align(structure_name, reference)

        standard_amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY','HIS', 'ILE', 'LEU', 
                                'LYS', 'MET', 'PHE', 'PRO', 'SER','THR', 'TRP', 'TYR', 'VAL', 'SEC', 'PYL'}

        # Extract HETATMs and HETRESs
        model = cmd.get_model(structure_name)
        hetatms = [atom for atom in model.atom if atom.hetatm and atom.resn.upper() not in set(['HOH']) and atom.resn.upper() not in standard_amino_acids]
        hetres = set([atom.resn.upper() for atom in hetatms])

        # For each HETRES, check distances vs all pockets
        for hr in sorted(hetres):

            # Get all HETRES atoms
            atoms = [atom for atom in hetatms if atom.resn.upper() == hr]
            atom_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in atoms])

            # For each pocket, check MIN distance
            for pocket in sorted(pocket_to_coords):
                min_distance = float(round(min(np.array([np.linalg.norm(i - pocket_to_coords[pocket]) for i in atom_coords])), 3))
                report.append([uni, structure_name, pocket, hr, min_distance])                

    # Hide all surface representations
    cmd.hide("surface", "all")

    # Center
    cmd.zoom(reference)

    # Save the session in the new directory
    cmd.save(session_output)

    return report


# Define file paths 
root = os.path.dirname(os.path.abspath(__file__))
pymol_sessions = os.path.join(root, "..", "output", "pymol_sessions")
pymol_sessions_pdbe = os.path.join(root, "..", "output", "pymol_sessions_pdbe")
pdb_structures = os.path.join(root, "..", "data", "structures", "pdbe_database")
os.makedirs(pymol_sessions_pdbe, exist_ok=True)

REPORT_ALL = []

# For all proteins having > 0 PDBe structures (both apo and holo)
for uni in sorted(os.listdir(pdb_structures)):

    # Enumerate PDB structures
    structures = sorted(os.listdir(os.path.join(pdb_structures, uni, uni + "_archive-PDB")))
    structures = [i for i in structures if ('.ent' in i or '.cif' in i) and '.txt' not in i]
    structures = [os.path.join(pdb_structures, uni, uni + "_archive-PDB", i) for i in structures]

    report = load_pymol_session_align(os.path.join(pymol_sessions, uni + ".pse.gz"), 
                       os.path.join(pymol_sessions_pdbe, uni + ".pse.gz"),
                       structures)
    
    if len(report) > 0:
        REPORT_ALL.extend(report)
    print(f"Created PyMol session for {uni}")

# Store PDBe report annotations
df = pd.DataFrame(REPORT_ALL, columns=["Uniprot ID", "PDB Structure", "Pocket", "HET RES", "Min Distance"])
df.to_csv(os.path.join(root, "..", "output" ,"pdbe_annotation_report.csv"), index=False, sep="\t")