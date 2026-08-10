import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
final_output_dir = output_dir

CATEGORY_COLUMNS = {
    "Catalytic Domain (ATP/ligase)": "Catalytic_Domain_support",
    "Anticodon Binding Domain": "Anticodon_Binding_support",
    "Editing Domain": "Editing_Domain_support",
    "tRNA Binding Domain": "tRNA_Binding_support",
    "RNA Binding Domain": "RNA_Binding_support",
    "Non-catalytic (structural homolog)": "NonCatalytic_support",
}


def parse_pocket_residues(residues_str):
    residues = set()
    for tok in residues_str.split():
        chain, resn = tok.split("_")
        residues.add(int(resn))
    return residues


def gene_map():
    bosch = pd.read_csv(os.path.join(repo_root, "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
    return dict(zip(bosch["uniprot_ac"], bosch["gene_name_in_bosch_2021"]))


def load_support(ac):
    path = os.path.join(output_dir, "{}_residue_support.csv".format(ac))
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, keep_default_na=False)
    m = {}
    for _, row in df.iterrows():
        m.setdefault(row["category"], {})[row["residue"]] = row["entry_support_count"]
    return m


def pocket_labels_and_support(pocket_residues, support_by_category):
    labels = []
    support_col_values = {}
    for category, col in CATEGORY_COLUMNS.items():
        cat_support = support_by_category.get(category, {})
        hits = [cat_support[r] for r in pocket_residues if r in cat_support]
        if hits:
            labels.append(category)
            support_col_values[col] = min(hits)  # weakest link: least-supported residue in the pocket
        else:
            support_col_values[col] = ""
    return labels, support_col_values


def agg_ligand_evidence(sub_df, cols):
    if sub_df.empty:
        return {c: "" for c in ["ligands"] + cols}
    out = {}
    out["ligands"] = "|".join(sub_df["ligand_resn"].astype(str))
    for c in cols:
        out[c] = "|".join(str(v) for v in sub_df[c])
    return out


if __name__ == "__main__":
    pockets_df = pd.read_csv(os.path.join(repo_root, "output", "pocket_detection_data.csv"))
    genes = gene_map()

    direct_ev = pd.read_csv(os.path.join(output_dir, "direct_pdb_ligand_evidence.csv"), keep_default_na=False)
    af_ev = pd.read_csv(os.path.join(output_dir, "alphafill_ligand_evidence.csv"), keep_default_na=False)
    clusters = pd.read_csv(os.path.join(output_dir, "pocket_clusters.csv"), keep_default_na=False)

    support_cache = {}
    rows = []
    for _, prow in pockets_df.iterrows():
        ac = prow["Uniprot AC"]
        file_name = prow["File name"]
        pocket_number = prow["Pocket number"]

        if ac not in support_cache:
            support_cache[ac] = load_support(ac)
        pocket_residues = parse_pocket_residues(prow["Pocket residues (chain_resn)"])
        labels, support_vals = pocket_labels_and_support(pocket_residues, support_cache[ac])

        d_sub = direct_ev[(direct_ev["uniprot_ac"] == ac) & (direct_ev["pocket_file_name"] == file_name) & (direct_ev["pocket_number"] == pocket_number)]
        a_sub = af_ev[(af_ev["uniprot_ac"] == ac) & (af_ev["pocket_file_name"] == file_name) & (af_ev["pocket_number"] == pocket_number)]

        d_agg = agg_ligand_evidence(d_sub, ["pdb_id", "strength", "distance"])
        a_agg = agg_ligand_evidence(a_sub, ["source_pdb", "strength", "distance", "local_rmsd", "homolog_identity"])

        cl = clusters[(clusters["Uniprot AC"] == ac) & (clusters["File name"] == file_name) & (clusters["Pocket number"] == pocket_number)]
        cluster_id = cl["spatial_cluster_id"].iloc[0] if len(cl) else ""
        cluster_size = cl["n_models_in_cluster"].iloc[0] if len(cl) else ""

        has_direct_strong = (d_sub["strength"] == "strong").any()
        has_af_strong = (a_sub["strength"] == "strong").any()
        has_any_weak = (d_sub["strength"] == "weak").any() or (a_sub["strength"] == "weak").any()
        has_catalytic_label = "Catalytic Domain (ATP/ligase)" in labels

        if not has_catalytic_label:
            confidence = 0
        elif has_direct_strong:
            confidence = 4
        elif has_af_strong:
            confidence = 3
        else:
            confidence = 1 + int(has_any_weak)

        row = {
            "Uniprot AC": ac, "Gene": genes.get(ac, ""), "File name": file_name,
            "Prediction type": prow["Prediction type"], "Pocket number": pocket_number,
            "Pocket score": prow["Pocket score"], "Pocket probability": prow["Pocket probability"],
            "curated_labels": "|".join(labels),
        }
        row.update(support_vals)
        row["has_direct_ligand_evidence"] = "yes" if not d_sub.empty else "no"
        row["direct_pdb_ligands"] = d_agg["ligands"]
        row["direct_pdb_ligand_sources"] = d_agg["pdb_id"]
        row["direct_pdb_ligand_strength"] = d_agg["strength"]
        row["direct_pdb_ligand_distance"] = d_agg["distance"]
        row["alphafill_ligands"] = a_agg["ligands"]
        row["alphafill_sources"] = a_agg["source_pdb"]
        row["alphafill_ligand_strength"] = a_agg["strength"]
        row["alphafill_ligand_distance"] = a_agg["distance"]
        row["alphafill_local_rmsd"] = a_agg["local_rmsd"]
        row["alphafill_identity"] = a_agg["homolog_identity"]
        row["spatial_cluster_id"] = cluster_id
        row["n_models_in_cluster"] = cluster_size
        row["catalytic_confidence"] = confidence

        rows.append(row)

    final_df = pd.DataFrame(rows)
    out_path = os.path.join(final_output_dir, "pocket_detection_interpro_updated.csv")
    final_df.to_csv(out_path, index=False)
    print("Saved final table with {} rows -> {}".format(len(final_df), out_path))
    print()
    print("catalytic_confidence distribution:")
    print(final_df["catalytic_confidence"].value_counts().sort_index())
    print()
    print("Proteins with zero pockets reaching confidence >= 3:")
    max_conf = final_df.groupby("Uniprot AC")["catalytic_confidence"].max()
    zero_high = max_conf[max_conf < 3]
    genes_rev = genes
    for ac, m in zero_high.items():
        print("  {} ({}) -> max confidence {}".format(ac, genes_rev.get(ac, "?"), m))
