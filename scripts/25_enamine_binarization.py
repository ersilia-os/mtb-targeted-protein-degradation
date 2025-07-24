from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
import os

root = os.path.dirname(os.path.abspath(__file__))

# Get input data
PATH_TO_REPORTS = os.path.join(root, "..", "processed", "unidock_docking", "docking_results")
PATH_TO_OUTPUT = os.path.join(root, "..", "processed", "unidock_docking", "binarized_reports")
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)


# For each pocket structure
for st in tqdm(sorted(os.listdir(PATH_TO_REPORTS))):

    # Read report
    report = pd.read_csv(os.path.join(PATH_TO_REPORTS, st, "report.csv"))

    # Define percentiles
    perc_01 = np.percentile(report['score'], 0.1)
    perc_05 = np.percentile(report['score'], 0.5)
    perc_1 = np.percentile(report['score'], 1)

    # Create three extra columns
    report["bin_01"] = (report["score"] < perc_01).astype(int)
    report["bin_05"] = (report["score"] < perc_05).astype(int)
    report["bin_1"] = (report["score"] < perc_1).astype(int)

    # Save report
    report.to_csv(os.path.join(PATH_TO_OUTPUT, f"report_bin_{st}.csv"), index=False)