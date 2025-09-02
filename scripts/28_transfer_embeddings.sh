#!/bin/bash

LOCAL_DIR="../processed/enamine_REAL_characterization/embeddings/"
REMOTE_USER="acomajuncosa"
REMOTE_HOST="10.7.108.17"   # dante IP
REMOTE_DIR="/aloy/home/acomajuncosa/Ersilia/mtb/processed/enamine_REAL_characterization/embeddings"

rsync -av --partial --progress \
  -e "ssh -T -c aes128-gcm@openssh.com -o Compression=no" \
  "$LOCAL_DIR" \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}"
