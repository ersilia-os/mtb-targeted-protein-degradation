import tarfile
import os

root = os.path.dirname(os.path.abspath(__file__))
UNIDOCK_PATH = os.path.join(root, "..", "processed", "unidock_docking", "docking_results")

# For each pocket
for pocket in sorted(os.listdir(UNIDOCK_PATH)):
    
    # Read log file inside tar.gz
    with tarfile.open(os.path.join(UNIDOCK_PATH, pocket, "logs.tar.gz"), "r:gz") as tar:
            for member in tar.getnames():
                logs = tar.extractfile("logs.log").readlines()

    # Iterate over lines
    for line in logs:
         line = line.decode('utf-8')
         if "error" in line.lower() or ("warning" in line.lower() and "in add_to_output_container" not in line):
              print(pocket, line)
              break