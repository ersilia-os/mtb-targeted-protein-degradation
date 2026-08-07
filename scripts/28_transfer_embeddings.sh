#!/bin/bash

LOCAL_DIR="/home/acomajuncosa/Documents_GPU/mtb-targeted-protein-degradation/output/enamine_REAL_characterization/embeddings/"  # results are gathered from norrsken-gpu-wsl
REMOTE_USER="acomajuncosa"
REMOTE_HOST="10.7.108.17"   # dante IP
REMOTE_DIR="/aloy/home/acomajuncosa/Ersilia/mtb/output/enamine_REAL_characterization/embeddings"

rsync -av --partial --progress \
  -e "ssh -T -c aes128-gcm@openssh.com -o Compression=no" \
  "$LOCAL_DIR" \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}"
