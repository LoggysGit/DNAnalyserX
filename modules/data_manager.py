import os

from pathlib import Path

import requests
import shutil

import gzip

import sqlite3

import modules.lib as lib

class DataManager:
    def __init__(self, db_path="genome_data.db"):
        self.db_path = db_path
        self.init_db()

    # - Database functions - #
    def get_db_connection(self): return sqlite3.connect(self.db_path)

    def init_db(self):
        with self.get_db_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("PRAGMA foreign_keys = ON;")

            # 1. Create Diseases table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS diseases (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL UNIQUE,
                    description TEXT,
                    significance TEXT,
                    pmid TEXT
                );
            """)

            # 2. Create Known Variants table
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS known_variants (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    ref_base TEXT NOT NULL,
                    alt_base TEXT NOT NULL,
                    disease_id INTEGER,
                    FOREIGN KEY (disease_id) REFERENCES diseases(id) ON DELETE SET NULL
                );
            """)

            cursor.execute("CREATE INDEX IF NOT EXISTS idx_coords ON known_variants(chromosome, position);")
            
            cursor.execute("SELECT COUNT(*) FROM diseases;")
            if cursor.fetchone()[0] == 0: self._seed_data(cursor)
                
            conn.commit()

    def _seed_data(self, cursor): pass

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
            
            #total_size_bytes = int(response.headers.get('content-length', 0))
            downloaded_bytes = 0
            chunk_size = 1024 * 64 
            
            # Safe binary download stream (Template from the internet)
            with open(target_path, "wb") as f_bin:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f_bin.write(chunk)
                        downloaded_bytes += len(chunk)

            lib.log(f"Chromosome #{chromosome_id} downloaded.")
            return gzip.open(target_path, "rt")
            
        except Exception as e:
            lib.log(f"Network processing system error: {e}")
            return None

    def purge_temp(self):
        if os.path.exists(lib.TEMP_DIR):
            try:
                for root, dirs, files in os.walk(lib.TEMP_DIR):
                    for f in files: os.remove(os.path.join(root, f))
                    for d in dirs: shutil.rmtree(os.path.join(root, d))
                lib.log("Temporary data buffer purged successfully.")
            except Exception as e: lib.log(f"Execution handling error during directory cleanup: {e}")