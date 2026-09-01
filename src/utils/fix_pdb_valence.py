"""
Removes PDB2PQR-mislabeled hydrogens that break RDKit's PDB valence sanitization -- a rare edge
case (confirmed on 4/38 human aaRS structures during script 94's receptor prep) where a hydrogen
ends up geometrically closer to, and RDKit-perceived as bonded to, the wrong heavy atom than the
one its own PDB atom name says it belongs to (e.g. an "HA" ending up bonded to "CB" instead of
"CA"). Confirmed empirically before writing this: the affected residues' heavy-atom geometry and
pLDDT are both normal in the source AlphaFold structure (CA-CB bond length ~1.54A, pLDDT 81) -- this
is a PDB2PQR hydrogen-placement artifact, not an input-structure problem. PDB2PQR's own log flags
these exact residues as "Unable to debump".

Run inside the "unidock_tools" env (needs RDKit).
"""
import argparse

from rdkit import Chem


def fix_valence_problems(input_pdb, output_pdb):
    """Removes any bonded hydrogen whose own PDB atom name doesn't match the valence-exceeded
    heavy atom it ended up perceived as bonded to. Returns the number of atoms removed (0 if the
    structure had no valence problems -- safe/idempotent no-op copy in that case)."""
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False, sanitize=False)
    problems = Chem.DetectChemistryProblems(mol)

    atoms_to_remove = set()
    for p in problems:
        if p.GetType() != "AtomValenceException":
            continue
        atom = mol.GetAtomWithIdx(p.GetAtomIdx())
        heavy_name = atom.GetPDBResidueInfo().GetName().strip()
        heavy_suffix = heavy_name[1:]
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "H":
                continue
            h_name = nbr.GetPDBResidueInfo().GetName().strip()
            if not h_name.startswith("H" + heavy_suffix):
                atoms_to_remove.add(nbr.GetIdx())

    with open(input_pdb) as f:
        lines = f.readlines()
    atom_line_idxs = [i for i, line in enumerate(lines) if line.startswith("ATOM") or line.startswith("HETATM")]
    remove_line_idxs = {atom_line_idxs[i] for i in atoms_to_remove}
    kept = [line for i, line in enumerate(lines) if i not in remove_line_idxs]

    with open(output_pdb, "w") as f:
        f.writelines(kept)

    return len(atoms_to_remove)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    n_removed = fix_valence_problems(args.input, args.output)
    print(f"Removed {n_removed} mislabeled hydrogen(s).")
