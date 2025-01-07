import requests
import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))


def download_alphafill_data(entry_id, output_dir):
    """
    Downloads AlphaFill data and JSON metadata for the given entry ID.
    
    Args:
        entry_id (str): The entry ID to download data for.
        output_dir (str): Directory to save the downloaded files.
    """
    # Base URLs
    base_url = f"https://alphafill.eu/v1/aff/{entry_id}"
    json_url = f"{base_url}/json"

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # File paths
    cif_file_path = os.path.join(output_dir, f"{entry_id}.cif")
    json_file_path = os.path.join(output_dir, f"{entry_id}.json")

    # Download AlphaFill data
    try:
        print(f"Downloading AlphaFill file for entry ID: {entry_id}")
        aff_response = requests.get(base_url, stream=True)
        aff_response.raise_for_status()
        with open(cif_file_path, "wb") as aff_file:
            for chunk in aff_response.iter_content(chunk_size=8192):
                aff_file.write(chunk)
        print(f"AlphaFill file saved to: {cif_file_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download AlphaFill file: {e}")
        return

    # Download JSON metadata
    try:
        print(f"Downloading JSON metadata for entry ID: {entry_id}")
        json_response = requests.get(json_url)
        json_response.raise_for_status()
        with open(json_file_path, "w") as json_file:
            json_file.write(json_response.text)
        print(f"JSON metadata saved to: {json_file_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download JSON metadata: {e}")


df = pd.read_csv(os.path.join(root, "..", "data", "mtb_trna_synthetases_bosch_2021_fig5_annotated.csv"))
for uniprot_ac in list(df["uniprot_ac"]):
    dirname = os.path.join(root, "..", "data", "structures", "alphafill_database", uniprot_ac)
    if not os.path.exists(dirname):
        os.makedirs(dirname, exist_ok=True)
        download_alphafill_data(uniprot_ac, dirname)