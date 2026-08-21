import os
import sys
import re
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.join(root, "..", "..")
output_dir = os.path.join(repo_root, "output", "77_pocket_annotation")
os.makedirs(output_dir, exist_ok=True)

UNIPROT_AC = sys.argv[1] if len(sys.argv) > 1 else "P9WFU3"  # default: pheS

# Explicit GO-id -> category map, ordered most- to least-specific.
# GO:0000166 (generic "nucleotide binding") is deliberately excluded: it is
# shared by structural domains (e.g. IPR010978 "tRNA-binding arm") and is not
# a reliable catalytic signal on its own.
# "Catalytic" also matches any GO term whose *name* fits the generic
# "<amino acid>-tRNA ligase/synthetase activity" pattern, since InterPro
# mints a distinct GO id per amino acid (e.g. GO:0004826 for phenylalanine)
# and enumerating all 20 by id would be brittle.
GO_ID_CATEGORY = {
    "GO:0004812": "Catalytic Domain (ATP/ligase)",   # aminoacyl-tRNA ligase activity
    "GO:0005524": "Catalytic Domain (ATP/ligase)",   # ATP binding
    "GO:0043039": "Catalytic Domain (ATP/ligase)",   # tRNA aminoacylation
    "GO:0006418": "Catalytic Domain (ATP/ligase)",   # tRNA aminoacylation for protein translation
    "GO:0043771": "Anticodon Binding Domain",        # anticodon binding
    "GO:0002161": "Editing Domain",                  # aminoacyl-tRNA editing activity
    "GO:0000049": "tRNA Binding Domain",              # tRNA binding
    "GO:0003723": "RNA Binding Domain",               # RNA binding
    # gatA/gatB (GatCAB amidotransferase complex, P9WQA1/P9WN61) use a
    # different catalytic mechanism than the other 19 targets (transamidate
    # a mischarged Glu-tRNA(Gln)/Asp-tRNA(Asn) using glutamine as the
    # amido-N-donor, rather than the amino-acid+ATP->aminoacyl-adenylate
    # ligase mechanism) — a different GO vocabulary, but genuinely the same
    # kind of thing: this protein's core catalytic active site.
    "GO:0050567": "Catalytic Domain (ATP/ligase)",   # glutaminyl-tRNA synthase (glutamine-hydrolyzing) activity (gatA)
    "GO:0016884": "Catalytic Domain (ATP/ligase)",   # carbon-nitrogen ligase activity, with glutamine as amido-N-donor (gatB)
    "GO:0016874": "Catalytic Domain (ATP/ligase)",   # ligase activity, attached to gatB's specific catalytic-domain entry (IPR006075)
}

CATALYTIC_NAME_PATTERN = re.compile(r"-tRNA (ligase|synthetase) activity", re.IGNORECASE)

CATEGORY_PRIORITY = [
    "Catalytic Domain (ATP/ligase)",
    "Anticodon Binding Domain",
    "Editing Domain",
    "tRNA Binding Domain",
    "RNA Binding Domain",
]

# Curated (uniprot_ac, accession) -> category overrides, checked before
# everything else. For cases verified against the literature where an
# entry's GO terms/name would otherwise mislabel it. Scoped per-protein
# (not per-accession) because the same InterPro entry can be the real
# catalytic module in one protein and a non-catalytic paralogous fold in
# another (see IPR045864 below).
#
# pheT (P9WFU1): bacterial PheRS is an alpha2beta2 heterotetramer. The
# alpha subunit (pheS) holds the sole catalytic module (CAM, domains A1/A2).
# The beta subunit's "core domain" is a structural paralog of that fold
# (CLM, "catalytic-like module", domains B6/B7) explicitly described as
# "catalytically not active", and its B5 domain is an unrelated DNA-binding
# -like (HTH) domain borrowed from BirA. Source: Klipcan et al. 2010,
# "Structural Aspects of Phenylalanylation and Quality Control in Three
# Major Forms of Phenylalanyl-tRNA Synthetase", J. Amino Acids
# (https://onlinelibrary.wiley.com/doi/10.4061/2010/983503); B5 domain,
# Wikipedia (https://en.wikipedia.org/wiki/B5_protein_domain).
#
# NOTE: IPR004532 (1-830) and IPR045060 (158-723) are whole/near-whole-chain
# "beta subunit" identity tags, not localized domain calls — the same
# "family entry spans nearly the whole protein" issue seen elsewhere in this
# project. They carry the whole-complex's GO terms (phenylalanine-tRNA
# ligase activity, ATP binding), so WITHOUT an explicit override here they
# fall to "Catalytic Domain (ATP/ligase)" via the normal GO-term path —
# reintroducing, across the whole 830-residue chain, the exact "pheT looks
# catalytic" mislabeling this investigation set out to fix. They must NOT be
# folded into "Non-catalytic (structural homolog)" either: that would make
# that category's span swallow the genuinely localized, literature-verified
# Anticodon Binding (735-831), Editing (211-407), and tRNA Binding (44-155)
# regions below. Both issues were caught via structural validation against
# PDB 7K98 (M. tuberculosis PheRS-tRNA complex) — explicitly override to
# "Other/Unclassified".
CURATED_OVERRIDES = {
    ("P9WFU1", "IPR004532"): "Other/Unclassified",                 # whole-chain identity tag, not locally catalytic
    ("P9WFU1", "IPR045060"): "Other/Unclassified",                 # whole-chain identity tag, not locally catalytic
    ("P9WFU1", "IPR045864"): "Non-catalytic (structural homolog)",  # CLM fold (B6/B7), catalytically inactive
    ("P9WFU1", "IPR041616"): "Non-catalytic (structural homolog)",  # CLM "core domain" (B6/B7), catalytically inactive
    ("P9WFU1", "IPR005147"): "Other/Unclassified",                 # B5 domain, DNA-binding-like HTH, unrelated to catalysis
}

# Strong name-based overrides, checked BEFORE GO terms. These patterns are
# narrow and highly specific enough that the entry's own name is more
# trustworthy than its GO terms, which InterPro often attaches at the
# whole-protein level rather than to this particular sub-domain (e.g. pheT's
# B3/B4 domain — the literature-confirmed editing domain of bacterial PheRS,
# see Kotik-Kogan et al. 2005 Structure — carries the parent protein's
# "phenylalanine-tRNA ligase activity" GO term even though it isn't the
# catalytic site).
STRONG_NAME_OVERRIDE_PATTERNS = [
    (re.compile(r"tRNA.binding arm", re.IGNORECASE), "tRNA Binding Domain"),
    (re.compile(r"anticodon.(binding|recognition)", re.IGNORECASE), "Anticodon Binding Domain"),
    (re.compile(r"\bediting\b", re.IGNORECASE), "Editing Domain"),
    (re.compile(r"B3.{0,2}/.{0,2}B?4", re.IGNORECASE), "Editing Domain"),  # bacterial PheRS beta-subunit editing domain (B3/B4)
]

# Broad, last-resort name fallback: only used when an entry has neither a
# strong name override nor any informative GO term at all.
# NOTE: requires "synthetase"/"ligase" to co-occur with "tRNA"/"transfer RNA"
# specifically — a bare "synthetase" substring match previously mislabeled
# gatB's IPR014746 ("Glutamine synthetase/guanido kinase, catalytic domain")
# as a tRNA-ligase-type Catalytic Domain, when glutamine synthetase is an
# unrelated enzyme family (EC 6.3.1.2) that merely shares a structural fold
# class — a coincidental name match, not evidence of aaRS-type activity.
NAME_FALLBACK_PATTERNS = [
    (re.compile(r"[a-z]_core\b", re.IGNORECASE), "Catalytic Domain (ATP/ligase)"),
    (re.compile(r"(tRNA|transfer RNA)[\s-]*(synthetase|ligase)|aminoacyl.{0,3}tRNA", re.IGNORECASE), "Catalytic Domain (ATP/ligase)"),
]


def categorize_by_strong_name_override(name):
    for pattern, category in STRONG_NAME_OVERRIDE_PATTERNS:
        if pattern.search(name):
            return category
    return None


def categorize_go_terms(go_terms_str):
    if not go_terms_str:
        return None, []
    hits = []
    for term in go_terms_str.split(";"):
        go_id, go_name = term.split("|", 1)
        if go_id in GO_ID_CATEGORY:
            hits.append((go_id, go_name, GO_ID_CATEGORY[go_id]))
        elif CATALYTIC_NAME_PATTERN.search(go_name):
            hits.append((go_id, go_name, "Catalytic Domain (ATP/ligase)"))
    if not hits:
        return None, []
    categories_hit = {h[2] for h in hits}
    for cat in CATEGORY_PRIORITY:
        if cat in categories_hit:
            return cat, hits
    return None, hits


def categorize_by_name_fallback(name):
    for pattern, category in NAME_FALLBACK_PATTERNS:
        if pattern.search(name):
            return category
    return None


def categorize(uniprot_ac, accession, name, go_terms_str):
    if (uniprot_ac, accession) in CURATED_OVERRIDES:
        return CURATED_OVERRIDES[(uniprot_ac, accession)], [], "curated override"
    cat = categorize_by_strong_name_override(name)
    if cat is not None:
        return cat, [], "name override"
    cat, hits = categorize_go_terms(go_terms_str)
    if cat is not None:
        return cat, hits, "GO term"
    cat = categorize_by_name_fallback(name)
    if cat is not None:
        return cat, hits, "name fallback"
    return "Other/Unclassified", hits, "none"


# aaRS Class I / Class II, detected from InterPro entry names (Eriani et al. 1990 classification,
# as discussed in Klipcan et al. 2010, https://doi.org/10.4061/2010/983503). Class II entries are
# unambiguous ("class II" never appears as a substring of anything else), but a naive "class I"
# substring match would false-positive on IPR010978, "Class I and II aminoacyl-tRNA synthetase,
# tRNA-binding arm" - a domain shared by both classes, present on many of these proteins regardless
# of their own class. Excluded via negative lookahead. gatA/gatB (not classical aaRSs - see the
# GatCAB amidotransferase note above) correctly fall out as "N/A" since neither has any
# class-mentioning InterPro entry at all; no hardcoded per-protein override needed.
AARS_CLASS_II_PATTERN = re.compile(r"\bclass[- ]?II\b", re.IGNORECASE)
AARS_CLASS_I_PATTERN = re.compile(r"\bclass[- ]?I\b(?!\s*and\s*II)", re.IGNORECASE)


def detect_aars_class(entry_names):
    is_class_i = any(AARS_CLASS_I_PATTERN.search(n) for n in entry_names)
    is_class_ii = any(AARS_CLASS_II_PATTERN.search(n) for n in entry_names)
    if is_class_i and is_class_ii:
        return "Ambiguous (both class I and II entries found)"
    if is_class_i:
        return "Class I"
    if is_class_ii:
        return "Class II"
    return "N/A (not a classical aaRS)"


if __name__ == "__main__":
    in_path = os.path.join(output_dir, "{}_annotation_table.csv".format(UNIPROT_AC))
    df = pd.read_csv(in_path, keep_default_na=False)

    categories = []
    matched_go = []
    category_sources = []
    for _, row in df.iterrows():
        cat, hits, source = categorize(UNIPROT_AC, row["accession"], row["name"], row["go_terms"])
        categories.append(cat)
        matched_go.append(";".join("{}|{}".format(h[0], h[1]) for h in hits))
        category_sources.append(source)

    df["category"] = categories
    df["category_matched_go_terms"] = matched_go
    df["category_source"] = category_sources

    aars_class = detect_aars_class(df["name"].tolist())
    df["aars_class"] = aars_class

    out_path = os.path.join(output_dir, "{}_annotation_table_categorized.csv".format(UNIPROT_AC))
    df.to_csv(out_path, index=False)
    print("Categorized {} entries -> {}".format(len(df), out_path))
    print("aaRS class: {}".format(aars_class))
    print(df[["accession", "name", "type", "specificity", "start", "end", "category"]].to_string(index=False))

    # entry_support_count: per (residue, category), how many distinct entries
    # in this table assign that category to that residue. A residue backed by
    # many independent domain-type entries is more trustworthy than one
    # backed only by a single whole-chain family entry (see e.g. pheS's
    # IPR045864/IPR022911 vs. its 6 narrower catalytic entries).
    support_rows = []
    for _, row in df.iterrows():
        if row["category"] == "Other/Unclassified":
            continue
        for r in range(row["start"], row["end"] + 1):
            support_rows.append((r, row["category"], row["accession"]))
    support_df = pd.DataFrame(support_rows, columns=["residue", "category", "accession"])
    support_counts = (
        support_df.groupby(["residue", "category"])["accession"]
        .nunique()
        .reset_index()
        .rename(columns={"accession": "entry_support_count"})
    )
    support_out_path = os.path.join(output_dir, "{}_residue_support.csv".format(UNIPROT_AC))
    support_counts.to_csv(support_out_path, index=False)
    print("Residue support counts ({} residue-category pairs) -> {}".format(len(support_counts), support_out_path))
