"""Shared ligand weak/strong classification rules for direct-PDB and AlphaFill
evidence extraction. Default is weak; "strong" is reserved for evidence that
specifically implicates the catalytic pocket of THIS protein's own reaction:
the free amino acid substrate, the ATP/ADP/AMP cofactor (or a non-hydrolyzable
analog), a per-protein reaction-intermediate mimic (aminoacyl-sulfamoyl-
adenosine / aminoacyl-AMP class), or a literature-recognized validated
inhibitor chemotype curated per protein.
"""

BUFFER_COMPONENTS = {
    "HOH", "SO4", "GOL", "MG", "CL", "NA", "K", "ZN", "CA", "MN", "PGE", "PGO", "PGR",
    "P6G", "ACT", "DMS", "EPE", "TRS", "PEG", "EDO", "BME", "IOD", "FMT", "MPD", "IMD",
}

# Each synthetase's own free amino acid substrate (3-letter PDB code).
SUBSTRATE_BY_AC = {
    "P9WFU3": "PHE",  # pheS
    "P9WFW3": "ASP",  # aspS
    "P9WFV1": "LEU",  # leuS
    "P9WFV9": "GLU",  # gltS
    "P9WFV7": "GLY",  # glyS
    "P9WFU5": "MET",  # metS
    "P9WFU9": "LYS",  # lysS
    "P9WFT3": "TRP",  # trpS
    "P9WFT1": "TYR",  # tyrS
    "P9WFT5": "THR",  # thrS
    "P9WFV3": "ILE",  # ileS
    "P9WFS9": "VAL",  # valS
    "P9WFT9": "PRO",  # proS
    "P9WFT7": "SER",  # serS
    "P9WFW1": "CYS",  # cysS1
    "P9WFV5": "HIS",  # hisS
    "P9WFW5": "ARG",  # argS
    "P9WFW7": "ALA",  # alaS
}

# ATP and its non-hydrolyzable analogs (the reaction cofactor for all 19
# classical aaRS's, and the phosphoryl donor in gatB's transamidation step).
COFACTOR_CODES = {"ATP", "ADP", "AMP", "APC"}  # APC = AMPCPP, non-hydrolyzable ATP analog

# Per-protein reaction-intermediate mimics: aminoacyl-sulfamoyl-adenosine or
# aminoacyl-AMP analogs, chemically the closest possible mimic of the actual
# catalytic transition state. Confirmed via RCSB chemical component lookup;
# only mapped for the protein whose own amino acid the compound carries (a
# homology-transplanted mimic for a DIFFERENT protein's substrate, e.g. YSA
# turning up on pheS via AlphaFill, does not count as strong for pheS).
INTERMEDIATE_ANALOG_BY_AC = {
    "P9WFU3": {"W5Y"},          # pheS: 5'-O-(L-phenylalanylsulfamoyl)adenosine
    "P9WFW3": {"DSZ"},          # aspS: 5'-O-(L-alpha-aspartylsulfamoyl)adenosine
    "P9WFV7": {"G5A"},          # glyS: 5'-O-(glycylsulfamoyl)adenosine
    "P9WFU9": {"KAA"},          # lysS: 5'-O-[(L-lysylamino)sulfonyl]adenosine
    "P9WFV1": {"LSS"},          # leuS: 5'-O-(L-leucylsulfamoyl)adenosine
    "P9WFT3": {"TYM"},          # trpS: tryptophanyl-5'-AMP
    "P9WFT1": {"YSA"},          # tyrS: 5'-O-[N-(L-tyrosyl)sulfamoyl]adenosine
    "P9WFU5": {"ME8"},          # metS: methionyl-adenylate analog
    "P9WFW7": {"A5A"},          # alaS: 5'-O-(N-(L-alanyl)-sulfamoyl)adenosine
    "P9WFT7": {"SSA"},          # serS: 5'-O-(N-(L-seryl)-sulfamoyl)adenosine
    "P9WFT5": {"TSB"},          # thrS: 5'-O-(N-(L-threonyl)-sulfamoyl)adenosine
}

# Literature-recognized validated inhibitor chemotypes, curated per protein
# (not a generic screening-fragment default -- those stay weak).
VALIDATED_INHIBITOR_BY_AC = {
    # leuS: benzoxaborole-modified adenosine class, related to the
    # AN3365/GSK2251052 antibacterial LeuRS editing-site inhibitor series.
    "P9WFV1": {"38Y", "A2H", "A52", "81D"},
}


def classify_ligand(resn, uniprot_ac):
    if resn in BUFFER_COMPONENTS:
        return "weak"
    if SUBSTRATE_BY_AC.get(uniprot_ac) == resn:
        return "strong"
    if resn in COFACTOR_CODES:
        return "strong"
    if resn in INTERMEDIATE_ANALOG_BY_AC.get(uniprot_ac, set()):
        return "strong"
    if resn in VALIDATED_INHIBITOR_BY_AC.get(uniprot_ac, set()):
        return "strong"
    return "weak"
