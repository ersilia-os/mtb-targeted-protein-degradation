"""Shared helpers for Enamine REAL 10B surrogate-screen selectivity scripts (57, 58).

The screen (ersilia-os/gcadda4tb-enamine-real-screening) never persists per-compound
probabilities: only the indices clearing the 99th percentile threshold per pocket per chunk
(ind_1.npz) are saved, alongside the scalar threshold itself.
"""
import io
import os
import time

import httplib2
import numpy as np
import pandas as pd
from google.oauth2.service_account import Credentials
from google_auth_httplib2 import AuthorizedHttp
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaIoBaseDownload

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
PATH_TO_POCKET_DATA = os.path.join(ROOT, "output", "pocket_detection_data.csv")
GDRIVE_FOLDER_ID = "1FBELagBf9hlKVgvkaZ8YF60jKRAmsHPo"


def download_file(outfile):
    service_file = os.path.join(ROOT, "service.json")
    file = os.path.basename(outfile)
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.readonly"])
    for attempt in range(10):
        try:
            http = httplib2.Http(timeout=600)
            authed_http = AuthorizedHttp(creds, http=http)
            service = build("drive", "v3", http=authed_http)
            query = f"name='{file}' and '{GDRIVE_FOLDER_ID}' in parents and trashed=false"
            results = service.files().list(q=query, fields="files(id)", supportsAllDrives=True, includeItemsFromAllDrives=True).execute()
            break
        except (HttpError, OSError):
            time.sleep(5)
            if attempt == 9:
                raise
    files = results.get("files", [])
    if not files:
        raise FileNotFoundError(f"'{file}' not found in folder {GDRIVE_FOLDER_ID}.")
    if len(files) > 1:
        raise RuntimeError(f"Multiple files named '{file}' are found...")
    file_id = files[0]["id"]
    request = service.files().get_media(fileId=file_id, supportsAllDrives=True)
    with io.FileIO(outfile, "wb") as fh:
        downloader = MediaIoBaseDownload(fh, request, chunksize=100 * 1024 * 1024)
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


def get_pocket_to_ac():
    """{pocket_name: Uniprot AC} for all 276 pockets, from output/pocket_detection_data.csv."""
    df = pd.read_csv(PATH_TO_POCKET_DATA)
    pocket_names = df["File name"].str.replace(".pdb", "", regex=False) + "_pocket_" + df["Pocket number"].astype(str)
    return dict(zip(pocket_names, df["Uniprot AC"]))


def load_ind1(tar, chunk, pocket):
    """Read a pocket's ind_1.npz directly from the chunk's tar archive, no disk extraction."""
    member = tar.extractfile(f"{chunk}/{pocket}_ind_1.npz")
    npz = np.load(io.BytesIO(member.read()))
    return npz["ind_1"], float(npz["thr"])
