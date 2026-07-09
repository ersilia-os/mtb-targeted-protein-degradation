"""Shared constants and helpers for docking hit analysis scripts (48_docking_hits_raw, 49_docking_hits_selective)."""
import os
import sys
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*converter.*already registered.*")

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

# Enamine REAL 10M round-1 screening (distinct from LIBRARIES["REAL"], which is round 2 / 10B).
# Each pocket's input_ligands file lists the prioritized ("active") compounds first, followed
# by the rest of the screened 10M library (the "negative set") — see notebooks/46_docking_exploration_I.ipynb.
REAL_ROUND1_RESULTS_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking", "docking_results")
REAL_ROUND1_INPUT_LIGANDS_DIR = os.path.join(ROOT, "output", "unidock_REAL_docking", "input_ligands")
REAL_ROUND1_N_ACTIVES = 100_000

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


def load_real_negative_scores(pocket_name):
    """Return a pd.Series of docking scores for the Enamine REAL 10M round-1 screening's
    negative set (compounds not among the top REAL_ROUND1_N_ACTIVES prioritized for this pocket).
    Empty Series if the round-1 files aren't available for this pocket."""
    input_ligands_path = os.path.join(REAL_ROUND1_INPUT_LIGANDS_DIR, f"input_ligands_{pocket_name}.txt")
    report_path = os.path.join(REAL_ROUND1_RESULTS_DIR, pocket_name, "report.csv")
    if not os.path.isfile(input_ligands_path) or not os.path.isfile(report_path):
        return pd.Series(dtype=float)
    with open(input_ligands_path) as f:
        compound_ids = [line.strip().replace(".sdf", "").split("/")[-1] for line in f]
    negative_ids = set(compound_ids[REAL_ROUND1_N_ACTIVES:])
    scores = load_scores(report_path)
    return scores[scores.index.isin(negative_ids)]


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


# --- UpSet plots ---

def save_upset(gene_top, top_n, output_dir, lib, trna_tag):
    """UpSet plot of top-N hit overlap across target genes."""
    import warnings
    import matplotlib.pyplot as plt
    from upsetplot import from_contents, plot as upset_plot
    warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

    top_sets = {g: set(ids[:top_n]) for g, ids in gene_top.items()}
    data = from_contents(top_sets)
    upset_plot(data)
    plt.suptitle(f"{lib} — top {top_n:,}")
    fname = f"{trna_tag}_{lib}_upset_top{top_n}.png"
    path = os.path.join(output_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {path}")


# --- Score boxplots ---

def simple_boxplot(ax, x, y, c, width):
    lw = 0.3
    ax.boxplot(y, positions=[x], widths=width, patch_artist=True, whis=[1, 99], showfliers=False,
               boxprops=dict(facecolor=c, color='black', lw=lw),
               medianprops=dict(color='black', lw=lw),
               whiskerprops=dict(color='black', lw=lw),
               capprops=dict(color='none', lw=lw))


def _safe_import_stylia():
    """Imports stylia, first ensuring matplotlib's cache dir exists - stylia's own __init__.py
    unconditionally runs shutil.rmtree(mpl.get_cachedir()), which raises FileNotFoundError if that
    directory doesn't already exist (e.g. right after a previous stylia import just deleted it, or
    on a machine/user account where matplotlib has never built its cache yet)."""
    import matplotlib as mpl
    os.makedirs(mpl.get_cachedir(), exist_ok=True)
    import stylia
    return stylia


def plot_score_boxplots(scores1, sel_ids, genes, out_path, sel_label="Selected", sel_color="purple",
                         dl_scores=None, real_negative_scores=None, xtick_rotation=0):
    """Boxplots of raw docking scores per gene: pre-screened vs selected, optionally alongside
    two external reference distributions for the same reference pockets:
      dl_scores            — gene-column DataFrame (same shape as scores1) of Enamine DL scores
      real_negative_scores — {gene: pd.Series} of Enamine REAL 10M round-1 negative-set scores
      xtick_rotation       — degrees to rotate x-tick labels (0 = horizontal, unrotated); use e.g.
                             45 when "genes" are actually long pocket names (script 53's NON-CAT
                             per-pocket columns), left at 0 for short gene-name labels (script 52).
    """
    import matplotlib.patches as mpatches
    stylia = _safe_import_stylia()

    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()
    c_sel = getattr(nc, sel_color)

    fig, axs = stylia.create_figure(1, 1, height=0.4, width=0.7)
    ax = axs.next()

    n_bg = len(scores1)
    n_sel = len(sel_ids)
    n_dl = len(dl_scores) if dl_scores is not None else 0
    n_rn = len(next(iter(real_negative_scores.values()))) if real_negative_scores else 0

    width = 0.25
    gap = 0.3
    group_spacing = 1.5
    offsets = [0, gap, 2 * gap, 3 * gap]  # DL, REAL negative set, pre-screened, selected

    ticks, tick_labels = [], []
    for i, gene in enumerate(genes):
        if gene not in scores1.columns:
            continue
        base = group_spacing * i
        bg_scores  = scores1[gene].dropna().values
        sel_scores = scores1.loc[scores1.index.isin(sel_ids), gene].dropna().values
        if dl_scores is not None and gene in dl_scores.columns:
            simple_boxplot(ax, base + offsets[0], dl_scores[gene].dropna().values, nc.mint, width)
        if real_negative_scores is not None and gene in real_negative_scores:
            simple_boxplot(ax, base + offsets[1], real_negative_scores[gene].dropna().values, nc.blue, width)
        simple_boxplot(ax, base + offsets[2], bg_scores,  nc.yellow, width)
        simple_boxplot(ax, base + offsets[3], sel_scores, c_sel,   width)
        ticks.append(base + sum(offsets) / len(offsets))
        tick_labels.append(gene)

    ax.set_xticks(ticks)
    if xtick_rotation:
        ax.set_xticklabels(tick_labels, rotation=xtick_rotation, ha="right")
    else:
        ax.set_xticklabels(tick_labels)
    legend_handles = []
    if dl_scores is not None:
        legend_handles.append(mpatches.Patch(facecolor=nc.mint, label=f"Enamine DL (n={n_dl:,})"))
    if real_negative_scores is not None:
        legend_handles.append(mpatches.Patch(facecolor=nc.blue, label=f"Enamine REAL negative set (n={n_rn:,})"))
    legend_handles.append(mpatches.Patch(facecolor=nc.yellow, label=f"Pre-screened (n={n_bg:,})"))
    legend_handles.append(mpatches.Patch(facecolor=c_sel,   label=f"{sel_label} (n={n_sel:,})"))
    ax.legend(handles=legend_handles)
    stylia.label(ax, xlabel="", ylabel="Docking score")
    stylia.save_figure(out_path)


# --- Physicochemical profiling ---

PROP_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings", "QED"]
DISCRETE_PROPS = {"HBD", "HBA", "RotBonds", "AromaticRings"}


def sample_prescreened_smiles(lib, exclude_ids, n, seed):
    """Return {id: smiles} for n randomly sampled compounds not in exclude_ids."""
    path = SMILES_PATHS[lib]
    if not os.path.isfile(path):
        print(f"  Warning: SMILES file not found at {path}, cannot sample pre-screened compounds.")
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


def plot_tsne(sel_smiles, bg_smiles, out_path, seed, sel_color="purple"):
    """
    t-SNE projection of ECFP4 (Morgan, radius=2, 2048 bits) fingerprints for a selected compound
    set vs. a background set - background plotted first (gray) so the selected set (colored)
    draws on top and stays visible. {id: smiles} dicts in, one scatter plot out.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from sklearn.manifold import TSNE
    stylia = _safe_import_stylia()

    # Format: slide | Style: ersilia - change with stylia.set_format() / stylia.set_style()
    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()
    c_sel = getattr(nc, sel_color)

    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    def to_fp(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return morgan_gen.GetFingerprintAsNumPy(mol)

    groups, fps = [], []
    for group, smiles_dict in [("bg", bg_smiles), ("sel", sel_smiles)]:
        for smi in smiles_dict.values():
            fp = to_fp(smi)
            if fp is None:
                continue
            fps.append(fp)
            groups.append(group)
    groups = np.array(groups)
    coords = TSNE(n_components=2, random_state=seed, init="random").fit_transform(np.array(fps))

    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    ax = axs.next()
    is_bg = groups == "bg"
    ax.scatter(coords[is_bg, 0], coords[is_bg, 1], color=nc.gray, label=f"Pre-screened (n={is_bg.sum():,})")
    ax.scatter(coords[~is_bg, 0], coords[~is_bg, 1], color=c_sel, label=f"Selected (n={(~is_bg).sum():,})")
    ax.legend()
    stylia.label(ax, xlabel="t-SNE 1", ylabel="t-SNE 2")
    stylia.save_figure(out_path)


def plot_profiling(sel_props, bg_props, out_path, sel_color="purple"):
    """KDE (continuous) + overlaid bars (discrete) + PAINS bar."""
    from scipy.stats import gaussian_kde
    stylia = _safe_import_stylia()

    stylia.set_format("slide")
    stylia.set_style("ersilia")
    nc = stylia.NamedColors()
    c_sel = getattr(nc, sel_color)

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
    print(f"  Pre-screened: {len(bg_props):,} compounds, {pains_bg:.1f}% PAINS")

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
            ax.bar(vals, bg_freq.values,  color=nc.gray,   alpha=1,   label="Pre-screened")
            ax.bar(vals, sel_freq.values, color=c_sel, alpha=0.5, label="Selected")
            ax.set_xlim(lo - 0.5, hi + 0.5)
            stylia.label(ax, xlabel=prop, ylabel="Frequency")
        else:
            for data, color, label in [
                (bg_data,  nc.gray,   f"Pre-screened (n={len(bg_props):,})"),
                (sel_data, c_sel, f"Selected (n={len(sel_props):,})"),
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
    ax.bar([0, 1], [pains_bg, pains_sel], color=[nc.gray, c_sel], alpha=0.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre-screened", "Selected"])
    ax.set_ylim(0, 10)
    stylia.label(ax, xlabel="", ylabel="PAINS (%)")

    stylia.save_figure(out_path)
