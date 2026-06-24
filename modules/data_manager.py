import os

from pathlib import Path

import requests
import shutil

import sqlite3

import modules.lib as lib

class DataManager:
    def __init__(self, db_path):
        self.db_path = db_path

    # - Directories & Network functions - #
    def download_chromosome(self, chromosome_id):        
        chrom_clean = str(chromosome_id).lower().replace("chr", "").strip()
        url = lib.CHROMOSOME_BASE_URL + f"chr{chrom_clean}.fna.gz"
        
        os.makedirs(lib.TEMP_DIR, exist_ok=True)
        target_path = os.path.join(lib.TEMP_DIR, f"chr{chrom_clean}.fasta.gz")
        
        lib.log(f"Downloading chromosome #{chromosome_id}...")

        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            # Safe binary download stream (Template from the internet)
            chunk_size = 1024 * 64 # 64 kb
            with open(target_path, "wb") as f_bin:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk: f_bin.write(chunk)

            lib.log(f"Chromosome #{chromosome_id} downloaded into {lib.TEMP_DIR}.")
            return lib.gzip_open(target_path)
            
        except Exception as e:
            lib.log(f"Download error: {e}")
            return None

    def purge_temp(self):
        if os.path.exists(lib.TEMP_DIR):
            try:
                for root, dirs, files in os.walk(lib.TEMP_DIR):
                    for f in files: os.remove(os.path.join(root, f))
                    for d in dirs: shutil.rmtree(os.path.join(root, d))

                lib.log("Temporary data buffer purged successfully.")

            except Exception as e: lib.log(f"Execution handling error during directory cleanup: {e}")