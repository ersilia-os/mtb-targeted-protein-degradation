RANDOM_SEED = 42

# aaRS Class I/II mapping (same mapping figure_1_plot.py/ExtendedDataFigure2.py already used,
# centralized here as a project-wide constant) - gatA/gatB are transamidases (GatCAB complex),
# not aaRS ligases, so they have no Class I/II designation and are intentionally excluded here.
AARS_CLASS_LABELS = {
    "alaS": "Class II", "argS": "Class I", "aspS": "Class II", "cysS1": "Class I",
    "gltS": "Class I", "glyS": "Class II", "hisS": "Class II", "ileS": "Class I",
    "leuS": "Class I", "lysS": "Class II", "metS": "Class I", "pheS": "Class II",
    "pheT": "Class II", "proS": "Class II", "serS": "Class II", "thrS": "Class II",
    "trpS": "Class I", "tyrS": "Class I", "valS": "Class I",
}
