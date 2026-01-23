# to be run in herbert
from scipy.spatial import distance
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

root = os.path.dirname(os.path.abspath(__file__))

perc = "1"

# Defining path to Enamine screening results and path to output
PATH_TO_RESULTS = "/aloy/home/acomajuncosa/Ersilia/gcadda4tb-enamine-real-screening/results"
PATH_TO_OUTPUT = os.path.join(root, "..", "processed", 'unidock_REAL_docking', 'inference_10B')
os.makedirs(os.path.join(PATH_TO_OUTPUT, "shared_compounds"), exist_ok=True)

# Load chunk info
CHUNKS = pd.read_csv("/aloy/home/acomajuncosa/Ersilia/gcadda4tb-enamine-real-screening/data/chunks/chunks.csv", header=None)[0].tolist()

# Load pocket info
POCKETS = sorted(set(sorted([i.split("_ind_")[0] for i in sorted(os.listdir(os.path.join(PATH_TO_RESULTS, "Enamine_REAL_LeadLike_000")))])))