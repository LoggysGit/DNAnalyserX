import os

from pathlib import Path
from datetime import datetime, timedelta

import gzip
import requests

import shutil
import sqlite3

import modules.lib as lib

class DiseaseDatabase:
    def __init__(self, db_path, gui_cmds):
        self.db_path = db_path
        self.gui_cmd_buffer = gui_cmds
        self._init_db()

    def _init_db(self):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()

            # Main DB
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS clinvar_mutations (
                    allele_id INTEGER PRIMARY KEY,
                    mutation_type TEXT,
                    chromosome TEXT,
                    position INTEGER,
                    ref_allele TEXT,
                    alt_allele TEXT,
                    clinical_significance TEXT,
                    disease_name TEXT
                );
            """)
            
            # Composite search index
            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_variant_search 
                ON clinvar_mutations (chromosome, position, ref_allele, alt_allele);
            """)

            # Metadata DB
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS db_metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT
                );
            """)
            
            conn.commit()

    def insert_mutation_batch(self, batch_data):
        query = """
            INSERT OR REPLACE INTO clinvar_mutations (
                allele_id, mutation_type, chromosome, position, 
                ref_allele, alt_allele, clinical_significance, disease_name
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?);
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.executemany(query, batch_data)
            conn.commit()

    def find_mutation(self, chromosome, position, ref, alt):
        query = """
            SELECT clinical_significance, disease_name 
            FROM clinvar_mutations 
            WHERE chromosome = ? AND position = ? AND ref_allele = ? AND alt_allele = ?
        """
        chrom_clean = str(chromosome).lower().replace("chr", "").strip()
        
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query, (chrom_clean, position, ref, alt))
            return cursor.fetchone()
        
    # - Metadata handle functions - #
    def get_last_update_date(self):
        query = "SELECT value FROM db_metadata WHERE key = 'last_update'"
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            row = cursor.fetchone()
            return row[0] if row else None

    def sys_set_last_update_now(self):
        now_str = datetime.now().date().isoformat() #'YYYY-MM-DD'
        query = "INSERT OR REPLACE INTO db_metadata (key, value) VALUES ('last_update', ?)"
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query, (now_str,))
            conn.commit()

class DataManager:
    def __init__(self, db_path):
        self.db_path = db_path
        self.disease_database = DiseaseDatabase(db_path)

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

    def update_disease_database(self):
        os.makedirs(lib.TEMP_DIR, exist_ok=True)
        target_path = os.path.join(lib.TEMP_DIR, "variant_summary.txt.gz")

        self.gui_cmd_buffer.put(["DB_UPDATE", None])
        
        lib.log("Starting ClinVar database download...")
        # Download
        try:
            response = requests.get(lib.DISEASE_BASE_URL, stream=True, timeout=60)
            response.raise_for_status()
            with open(target_path, "wb") as f_bin:
                for chunk in response.iter_content(chunk_size=1024 * 64):
                    if chunk: f_bin.write(chunk)
            lib.log("ClinVar archive downloaded successfully.")
        except Exception as e:
            lib.log(f"Failed to download ClinVar update archive: {e}")
            return

        lib.log("Beginning parsing pipeline into SQLite...")
        batch = []
        batch_size = 50000
        records_processed = 0

        # Parse
        try:
            with gzip.open(target_path, "rt", encoding="utf-8") as f_text:
                # Skip column descriptions row
                _ = f_text.readline()  
                
                for line in f_text:
                    tokens = line.split("\t")
                    if len(tokens) < 34: continue
                    
                    # FIlter old gene formats
                    assembly = tokens[16].strip()
                    if assembly != lib.GENOME_VER: continue
                        
                    allele_id = int(tokens[0])
                    mut_type = tokens[1].strip()
                    chromosome = tokens[18].strip().lower().replace("chr", "")
                    
                    # Skip unmapped coordinate tracking locations
                    try: position = int(tokens[30])
                    except ValueError: continue
                        
                    ref_allele = tokens[31].strip().upper()
                    alt_allele = tokens[32].strip().upper()
                    clinical_sig = tokens[6].strip()
                    disease_name = tokens[13].strip()

                    batch.append((
                        allele_id, mut_type, chromosome, position, 
                        ref_allele, alt_allele, clinical_sig, disease_name
                    ))

                    if len(batch) >= batch_size:
                        self.disease_database.insert_mutation_batch(batch)
                        records_processed += len(batch)
                        lib.log(f"Parsed and integrated {records_processed} target entries...")
                        batch = []

                # Insert into DB
                if batch:
                    self.disease_database.insert_mutation_batch(batch)
                    records_processed += len(batch)

            lib.log(f"Database generation completed. Integrated {records_processed} variants.")
        except Exception as e: lib.log(f"Processing error during database compilation pipeline: {e}")
        # Clear source file
        finally: self.purge_temp()
        
        self.gui_cmd_buffer.put(["DONE", None])

    def handle_disease_db_update(self):
        last_update_str = self.disease_database.get_last_update_date()

        # Check time treshhold
        if last_update_str:
            last_update = datetime.fromisoformat(last_update_str).date()
            if datetime.now().date() - last_update < timedelta(days=7): return

        # Update DB
        lib.log("Weekly update threshold reached. Initializing pipeline...")
        self.update_disease_database()

        # Save update timestamp
        self.disease_database.sys_set_last_update_now()