"""Shared constants and helpers for docking hit analysis scripts (48_docking_hits_raw, 49_docking_hits_selective)."""
import os
import sys

import numpy as np
import pandas as pd
from rdkit.Chem import MolFromSmiles, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.QED import qed as calc_qed
from tqdm import tqdm

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

LIBRARIES = {
    "DL":   os.path.join(ROOT, "output", "unidock_docking",        "docking_results"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking_2", "docking_results"),
}

REFERENCE_POCKET_CSV = os.path.join(ROOT, "output", "reference_pocket.csv")

SMILES_PATHS = {
    "DL":   os.path.join(ROOT, "data", "Enamine", "Enamine_Hit_Locator_Library_100K_Set_plated_100160cmpds_20250629.smiles"),
    "REAL": os.path.join(ROOT, "output", "unidock_REAL_docking", "inference_10B", "clustered_compounds.csv"),
}

SMILES_SEPS = {
    "DL":   "\t",
    "REAL": ",",
}

SMILES_COLS = {
    "DL":   ("Catalog ID", "SMILES"),
    "REAL": ("id",         "smiles"),
}


def lookup_smiles(compound_ids, lib):
    """Stream the library SMILES file in chunks, returning {id: smiles} for the requested ids."""
    path = SMILES_PATHS[lib]
    if not os.path.isfile(path):
        print(f"  Warning: SMILES file not found at {path}, omitting SMILES column.")
        return {}
    id_col, smi_col = SMILES_COLS[lib]
    sep = SMILES_SEPS[lib]
    needed = set(compound_ids)
    found = {}
    for chunk in pd.read_csv(path, sep=sep, usecols=[id_col, smi_col], chunksize=200_000):
        match = chunk[chunk[id_col].isin(needed)]
        found.update(zip(match[id_col], match[smi_col]))
        if len(found) == len(needed):
            break
    return found


def load_gene_map():
    path = os.path.join(ROOT, "data", "mtb_trna_synthetases_bosch_2021_fig5.csv")
    df = pd.read_csv(path)
    return {row["gene_name_in_bosch_2021"]: row["uniprot_ac"] for _, row in df.iterrows()}


def load_reference_pockets():
    if not os.path.isfile(REFERENCE_POCKET_CSV):
        print(f"Error: {REFERENCE_POCKET_CSV} not found.")
        print("Run script 47 and manually populate output/reference_pocket.csv")
        print("  (columns: gene_name, pocket_name)")
        sys.exit(1)
    df = pd.read_csv(REFERENCE_POCKET_CSV)
    return dict(zip(df["gene_name"], df["pocket_name"]))


def load_scores(pocket_path):
    """Return a pd.Series of scores indexed by compound id."""
    df = pd.read_csv(pocket_path)
    return df.set_index("compound")["score"]


def build_matrix(pocket_map, results_dir, label=""):
    """
    Build a raw-score DataFrame from a dict {column_label: pocket_name}.
    Rows = compounds, columns = column_labels. Missing cells → NaN.
    """
    series = {}
    for col_label, pocket_name in tqdm(pocket_map.items(), desc=label, unit="pocket"):
        report = os.path.join(results_dir, pocket_name, "report.csv")
        if not os.path.isfile(report):
            print(f"  Warning: report not found for pocket '{pocket_name}', skipping.")
            continue
        series[col_label] = load_scores(report)
    if not series:
        return pd.DataFrame()
    return pd.concat(series, axis=1)


# --- Physicochemical profiling ---

PROP_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED"]
DISCRETE_PROPS = {"HBD", "HBA", "RotBonds", "AromaticRings"}


def sample_background_smiles(lib, exclude_ids, n, seed):
    """Return {id: smiles} for n randomly sampled compounds not in exclude_ids."""
    path = SMILES_PATHS[lib]
    if not os.path.isfile(path):
        print(f"  Warning: SMILES file not found at {path}, cannot sample background.")
        return {}
    id_col, smi_col = SMILES_COLS[lib]
    sep = SMILES_SEPS[lib]
    pool = {}
    for chunk in pd.read_csv(path, sep=sep, usecols=[id_col, smi_col], chunksize=200_000):
        cands = chunk[~chunk[id_col].isin(exclude_ids)]
        pool.update(zip(cands[id_col], cands[smi_col]))
    ids = list(pool.keys())
    sampled = np.random.default_rng(seed).choice(ids, size=min(n, len(ids)), replace=False)
    return {cid: pool[cid] for cid in sampled}


def compute_properties(smiles_dict):
    """Compute physicochemical properties for a {id: smiles} dict."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)
    records = []
    for cid, smi in tqdm(smiles_dict.items(), desc="Computing properties", unit="mol"):
        if not isinstance(smi, str):
            continue
        mol = MolFromSmiles(smi)
        if mol is None:
            continue
        records.append({
            "id": cid,
            "MW":            rdMolDescriptors.CalcExactMolWt(mol),
            "cLogP":         rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
            "TPSA":          rdMolDescriptors.CalcTPSA(mol),
            "HBD":           rdMolDescriptors.CalcNumHBD(mol),
            "HBA":           rdMolDescriptors.CalcNumHBA(mol),
            "RotBonds":      rdMolDescriptors.CalcNumRotatableBonds(mol),
            "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "QED":           calc_qed(mol),
            "is_pains":      catalog.HasMatch(mol),
        })
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).set_index("id")


def plot_profiling(sel_props, bg_props, out_path):
    """KDE (continuous) + overlaid bars (discrete) + PAINS bar."""
    from scipy.stats import gaussian_kde
    import stylia

    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()

    XLIMS = {
        "MW":            (200, 500),
        "cLogP":         (-4, 8),
        "TPSA":          (0, 200),
        "QED":           (0, 1),
        "HBD":           (0, 7),
        "HBA":           (0, 13),
        "RotBonds":      (0, 13),
        "AromaticRings": (0, 8),
    }

    pains_sel = 100 * sel_props["is_pains"].sum() / len(sel_props) if len(sel_props) > 0 else 0.0
    pains_bg  = 100 * bg_props["is_pains"].sum()  / len(bg_props)  if len(bg_props)  > 0 else 0.0
    print(f"  Selected  : {len(sel_props):,} compounds, {pains_sel:.1f}% PAINS")
    print(f"  Background: {len(bg_props):,} compounds, {pains_bg:.1f}% PAINS")

    fig, axs = stylia.create_figure(3, 3, width=1.3, height=0.8)

    for prop in PROP_COLUMNS:
        ax = axs.next()
        bg_data  = bg_props[prop].dropna().values
        sel_data = sel_props[prop].dropna().values
        lo, hi = XLIMS[prop]
        if prop in DISCRETE_PROPS:
            vals = sorted(set(bg_data.astype(int)) | set(sel_data.astype(int)))
            bg_freq  = pd.Series(bg_data.astype(int)).value_counts(normalize=True).reindex(vals, fill_value=0)
            sel_freq = pd.Series(sel_data.astype(int)).value_counts(normalize=True).reindex(vals, fill_value=0)
            ax.bar(vals, bg_freq.values,  color=nc.gray,   alpha=1,   label="Background")
            ax.bar(vals, sel_freq.values, color=nc.purple, alpha=0.5, label="Selected")
            ax.set_xlim(lo - 0.5, hi + 0.5)
            stylia.label(ax, xlabel=prop, ylabel="Frequency")
        else:
            for data, color, label in [
                (bg_data,  nc.gray,   f"Background (n={len(bg_props):,})"),
                (sel_data, nc.purple, f"Selected (n={len(sel_props):,})"),
            ]:
                if len(data) < 2:
                    continue
                kde = gaussian_kde(data)
                x = np.linspace(lo, hi, 300)
                ax.plot(x, kde(x), color=color, label=label)
            ax.set_xlim(lo, hi)
            stylia.label(ax, xlabel=prop, ylabel="Density")
            if prop == "cLogP":
                ax.legend(loc="upper left")

    ax = axs.next()
    ax.bar([0, 1], [pains_bg, pains_sel], color=[nc.gray, nc.purple], alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Background", "Selected"])
    ax.set_ylim(0, 10)
    stylia.label(ax, xlabel="", ylabel="PAINS (%)")

    stylia.save_figure(out_path)
