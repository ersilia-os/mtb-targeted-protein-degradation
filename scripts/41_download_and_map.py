from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
from google_auth_httplib2 import AuthorizedHttp
from googleapiclient.errors import HttpError
from collections import Counter
from tqdm import tqdm
import pandas as pd
import httplib2
import time
import os
import io

def download_file(outfile):
    root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"
    file = os.path.basename(outfile)
    service_file = os.path.join(root, "service.json")
    GDRIVE_FOLDER_ID = "1FBELagBf9hlKVgvkaZ8YF60jKRAmsHPo"
    folder_id = GDRIVE_FOLDER_ID
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.readonly"])
    for attempt in range(10):
        try:
            http = httplib2.Http(timeout=600)
            authed_http = AuthorizedHttp(creds, http=http)
            service = build("drive", "v3", http=authed_http)
            query = f"name='{file}' and '{folder_id}' in parents and trashed=false"
            results = service.files().list(q=query, fields="files(id)", supportsAllDrives=True, includeItemsFromAllDrives=True).execute()
            break
        except (HttpError, OSError):
            time.sleep(5)
            if attempt == 9:
                raise
    files = results.get("files", [])
    if not files:
        raise FileNotFoundError(f"'{file}' not found in folder {folder_id}. Consider checking available chunks in data/chunks/chunks.csv")
    if len(files) > 1:
        raise RuntimeError(f"Multiple files named '{file}' are found...")
    file_id = files[0]["id"]
    request = service.files().get_media(fileId=file_id, supportsAllDrives=True)
    with io.FileIO(outfile, "wb") as fh:
        downloader = MediaIoBaseDownload(fh, request,chunksize=100 * 1024 * 1024)  # 100 MB chunks
        done = False
        retries = 0
        while not done:
            try:
                status, done = downloader.next_chunk()
                if status:
                    print(f"Download {int(status.progress() * 100)}%\n")
            except (HttpError, OSError):
                retries += 1
                print(f"Error found when downloading file. Trying again...[{retries}/10]")
                if retries >= 10:
                    raise
                time.sleep(10)

def get_splits(path):
    splits = os.listdir(path)
    return [i.replace(".csv", "") for i in sorted(splits)]

root = "/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation"
PATH_TO_OUTPUT = os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B")

# Get splits (994)
SPLITS = get_splits(os.path.join(PATH_TO_OUTPUT, 'A_proteins'))


A_pockets = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "A_pockets.csv"))
B_pockets = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "B_pockets.csv"))
A_proteins = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "A_proteins.csv"))
B_proteins = pd.read_csv(os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "B_proteins.csv"))

# Check that all pockets have same number of cpds associated
assert set(Counter(B_pockets['pocket']).values()) == set([1000])
assert set(Counter(B_proteins['pocket']).values()) == set([13000])

# Drop pocket column
B_pockets = B_pockets.drop(columns=['pocket'])
B_proteins = B_proteins.drop(columns=['pocket'])

# Include label
A_pockets['label'] = "A_pockets" 
B_pockets['label'] = "B_pockets" 
A_proteins['label'] = "A_proteins" 
B_proteins['label'] = "B_proteins"

# Merge all, sort and drop duplicates
COMPOUNDS = pd.concat([A_pockets, B_pockets, A_proteins, B_proteins]).sort_values(by=['split', 'index']).drop_duplicates(subset=['split', 'index']).reset_index(drop=True)

print(f"Number of selected compounds {len(A_pockets) + len(B_pockets) + len(A_proteins) + len(B_proteins)}")
print(f"Number of unique compounds {len(COMPOUNDS)}")
print(f"Intersection between sets OR duplication within sets: {(len(A_pockets) + len(B_pockets) + len(A_proteins) + len(B_proteins)) - len(COMPOUNDS)}")

# Create output directory
PATH_TO_OUTPUT = os.path.join(root, "processed", "unidock_REAL_docking", "inference_10B", "selected_compounds")
os.makedirs(PATH_TO_OUTPUT, exist_ok=True)

for split in tqdm(SPLITS):

    PATH_TO_IDs = os.path.join(root, "tmp", f"{split}_SMILES_IDs.tsv.zip")

    # Download splits id file
    if os.path.exists(PATH_TO_IDs) == False:
        download_file(PATH_TO_IDs)

    # Identify subset
    df = COMPOUNDS[COMPOUNDS['split'] == split].reset_index(drop=True)
    inds = df['index'].to_numpy()

    # Load split info
    smiles_ids = pd.read_csv(PATH_TO_IDs, sep='\t')
    smiles_ids = smiles_ids.iloc[inds].reset_index()
    assert (smiles_ids['index'] == df['index']).all

    # Concatanate horizontally
    df = pd.concat([df, smiles_ids], axis=1)

    # Save
    df.to_csv(os.path.join(PATH_TO_OUTPUT, f"{split}.csv"), index=False)