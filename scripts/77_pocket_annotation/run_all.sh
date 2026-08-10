#!/bin/bash
# Runs the full pocket/domain annotation pipeline for all 21 tRNA synthetases.
# Steps 06, 07, 08, 11 use PyMOL (adda4tb conda env); the rest use plain
# python3 (pandas/requests/scipy, no PyMOL needed).
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

eval "$(conda shell.bash hook)"

ACS=$(python3 -c "
import pandas as pd
df = pd.read_csv('$REPO_ROOT/data/mtb_trna_synthetases_bosch_2021_fig5_annotated.csv')
print(' '.join(df['uniprot_ac']))
")

echo "=== Step 1: fetch InterPro data for all 21 proteins ==="
for ac in $ACS; do
    python3 "$SCRIPT_DIR/01_fetch_interpro.py" "$ac"
done

echo "=== Step 2-3: build + categorize annotation tables ==="
for ac in $ACS; do
    python3 "$SCRIPT_DIR/02_build_annotation_table.py" "$ac"
    python3 "$SCRIPT_DIR/03_categorize.py" "$ac"
done

echo "=== Step 4: query PDB cross-references ==="
python3 "$SCRIPT_DIR/04_query_pdb_xrefs.py"

echo "=== Step 5: download experimental PDB structures ==="
python3 "$SCRIPT_DIR/05_download_pdb_structures.py"

conda activate adda4tb

echo "=== Step 6: identify chains + numbering offsets ==="
python3 "$SCRIPT_DIR/06_identify_chains.py"

echo "=== Step 7: pocket-local alignment + direct-PDB ligand evidence ==="
python3 "$SCRIPT_DIR/07_align_and_extract_ligands.py"

echo "=== Step 8: AlphaFill ligand evidence ==="
python3 "$SCRIPT_DIR/08_extract_alphafill_evidence.py"

conda deactivate

echo "=== Step 9: 3D pocket clustering ==="
python3 "$SCRIPT_DIR/09_cluster_pockets.py"

echo "=== Step 10: assemble final table -> output/77_pocket_annotation/pocket_detection_interpro_updated.csv ==="
python3 "$SCRIPT_DIR/10_assemble_final_table.py"

conda activate adda4tb

echo "=== Step 11: build per-protein PyMOL sessions ==="
python3 "$SCRIPT_DIR/11_build_pymol_sessions.py"

echo "Done."
