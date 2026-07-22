""" Data manager - controls files and DB"""

import os
import re
import time

from pathlib import Path
from datetime import datetime, timedelta

import gzip
import shutil
import sqlite3

import requests

from modules import lib

class Database:
    """ Database class"""

    def __init__(self, db_path, gui_cmds):
        self.db_path = db_path
        self.gui_cmd_buffer = gui_cmds

        self.init_db()

    def init_db(self):
        """ Creates and inits a DB"""
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()

            cursor.execute("""
                CREATE TABLE IF NOT EXISTS Mutations (
                    hgvs TEXT,
                    chromosome TEXT,
                    gene TEXT,
                    position INTEGER,
                    ref_allele TEXT,
                    alt_allele TEXT,
                    clnvs TEXT,
                    clinical_significance TEXT,
                    disease_name TEXT,
                    UNIQUE(gene, hgvs, alt_allele)
                );
            """)

            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_lookup 
                ON Mutations (gene, hgvs, alt_allele);
            """)

            cursor.execute("""
                CREATE TABLE IF NOT EXISTS Metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT
                );
            """)

            db.commit()

            lib.log("Disease database initialized.")

    def insert_mutation_batch(self, batch_data):
        """ Inserts data in DB """
        query = """
            INSERT OR REPLACE INTO Mutations (
                hgvs, chromosome, gene, position, ref_allele, 
                alt_allele, clnvs, clinical_significance, disease_name
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
        """
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.executemany(query, batch_data)
            db.commit()

    def find_mutation(self, hgvs_id, gene, alt) -> list:
        """ Filters and returns a mutation"""
        query = """
            SELECT chromosome, position, ref_allele, alt_allele, clnvs, clinical_significance, disease_name 
            FROM Mutations WHERE gene = ? AND hgvs = ? AND alt_allele = ?
        """
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(query, (str(gene).strip(), str(hgvs_id).strip(), str(alt).strip()))

            result = cursor.fetchone()
            return result if result else [None, None, None, alt, None, "Unknown", "Unknown"]

    # = Metadata control = #
    def get_last_update_date(self) -> str:
        """ Returns last update date value"""
        query = "SELECT value FROM Metadata WHERE key = 'last_update'"
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(query)

            row = cursor.fetchone()
            return row[0] if row else None

    def sys_set_last_update_now(self):
        """ Saves last update date in Metadata Table """
        now_str = datetime.now().date().isoformat()
        query = "INSERT OR REPLACE INTO Metadata (key, value) VALUES ('last_update', ?)"
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(query, (now_str,))
            db.commit()

class DataManager:
    """ Main data manager"""

    def __init__(self, db_path, gui_cmds):
        self.db_path = db_path
        self.gui_cmd_buffer = gui_cmds

        self.disease_database = Database(db_path, gui_cmds)

    def purge_temp(self):
        """ Purges temp directory """
        if os.path.exists(lib.TEMP_DIR):
            try:
                for root, dirs, files in os.walk(lib.TEMP_DIR):
                    for f in files:
                        os.remove(os.path.join(root, f))
                    for d in dirs:
                        shutil.rmtree(os.path.join(root, d))

                lib.log("Temporary data buffer purged successfully.")

            except Exception as e:
                lib.log(f"Cleanup error: {e}")

    def save_mutations_to_vcf(self, export_path, mutation_list, reference_name=""):
        """ Saves mutation list into VCF file """
        path = Path(export_path)

        header_lines = [
            "##fileformat=VCFv4.2",
            f"##fileDate={time.strftime('%Y%m%d')}",
            "##source=DNAnalyserX_Engine",
            f"##reference={reference_name}",
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">',
            '##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS protein notation">',
            '##INFO=<ID=SIG,Number=1,Type=String,Description="Clinical significance from ClinVar">',
            '##INFO=<ID=DIS,Number=1,Type=String,Description="Associated disease name">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        ]

        def safe_str(value, default="."):
            if value is None:
                return default
            text = str(value).strip()
            return text if text else default

        def sanitize_str(value, default="."):
            text = safe_str(value, default)
            return text.replace(" ", "_").replace(";", "_").replace("=", "_").replace(",", "_")

        try:
            rows = []
            for item in mutation_list:
                chrom = safe_str(item.get("chr"))
                gene = safe_str(item.get("gene"))
                hgvs = safe_str(item.get("hgvs"))
                ref = safe_str(item.get("ref"))
                alt = safe_str(item.get("alt"))
                clnvs = sanitize_str(item.get("clnvs", "UNKNOWN"))
                significance = sanitize_str(item.get("clnsign", "-"))
                disease_name = sanitize_str(item.get("name", "Not_found"))

                try:
                    pos = int(item.get("position"))
                    if pos <= 0:
                        raise ValueError

                except (TypeError, ValueError):
                    lib.log(f"Skipping VCF row: invalid POS "
                            f"'{item.get('position')}' for {gene} {hgvs}")
                    skipped += 1
                    continue

                if ref == "." or ref == "":
                    ref = "N"
                if alt == "." or alt == "":
                    alt = "N"

                info_block = (f"SVTYPE={clnvs};SIG={significance};"
                              f"DIS={disease_name};GENE={gene};HGVS={hgvs}")
                rows.append((chrom, pos, ref, alt, info_block))

            # Sort rows by position
            rows.sort(key=lambda r: (r[0], r[1]))

            with open(path, "w", encoding="utf-8") as f_out:
                for header in header_lines:
                    f_out.write(header + "\n")

                for chrom, pos, ref, alt, info_block in rows:
                    vcf_row = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info_block}\n"
                    f_out.write(vcf_row)

                lib.log(f"VCF data exported in {path}.")
                return True

        except Exception as e:
            lib.log(f"VCF export critical failure: {e}")
            return False

    def update_disease_database(self):
        """ Updates disease DB """
        os.makedirs(lib.TEMP_DIR, exist_ok=True)
        target_path = os.path.join(lib.TEMP_DIR, "disease_variant_summary.txt.gz")

        self.gui_cmd_buffer.put(["DB_UPDATE", None])

        # Download file
        lib.log("Starting ClinVar database download...")
        try:
            response = requests.get(lib.DISEASE_BASE_URL, stream=True, timeout=60)
            response.raise_for_status()
            with open(target_path, "wb") as f_bin:
                for chunk in response.iter_content(chunk_size=1024 * 64):
                    if chunk:
                        f_bin.write(chunk)
            lib.log("ClinVar database downloaded successfully.")

        except Exception as e:
            lib.log(f"Failed to download ClinVar database: {e}")
            return

        # Parsing downloaded
        lib.log("Parsing into SQLite DB...")
        batch = []
        batch_size = 50000
        records_processed = 0

        aa_map = {
            "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
            "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
            "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
            "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
            "Ter": "*", "Xaa": "X"
        }

        hgvs_regex = re.compile(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*|fs)?\)")

        try:
            with gzip.open(target_path, "rt", encoding="utf-8") as f_text:
                _ = f_text.readline()

                for line in f_text:
                    tokens = line.split("\t")
                    if len(tokens) < 34:
                        continue

                    # Only actual genome version
                    if tokens[16].strip() != lib.GENOME_VER:
                        continue

                    # Check gene data
                    gene_symbol = tokens[4].strip()
                    if not gene_symbol or gene_symbol == "-":
                        continue

                    variant_name = tokens[2].strip()
                    hgvs_match = hgvs_regex.search(variant_name)
                    if not hgvs_match:
                        continue

                    ref_aa_3, pos_num, alt_aa_3 = hgvs_match.groups()

                    ref_aa_1 = aa_map.get(ref_aa_3)
                    if not ref_aa_1:
                        continue

                    if alt_aa_3 == "fs":
                        alt_aa_1 = "fs"
                    elif alt_aa_3:
                        alt_aa_1 = aa_map.get(alt_aa_3, alt_aa_3)
                    else: alt_aa_1 = ""

                    hgvs_key = f"p.{ref_aa_1}{pos_num}{alt_aa_1}"

                    chromosome = tokens[18].strip().lower().replace("chr", "")

                    try:
                        position = int(tokens[19])
                    except ValueError:
                        continue

                    ref_allele = tokens[32].strip().upper()
                    alt_allele = tokens[33].strip().upper()

                    clinical_sig = tokens[6].strip()
                    clnvs = tokens[1].strip()
                    disease_name = tokens[13].strip()

                    batch.append((
                        hgvs_key, chromosome, gene_symbol, position,
                        ref_allele, alt_allele, clnvs, clinical_sig, disease_name
                    ))

                    if len(batch) >= batch_size:
                        self.disease_database.insert_mutation_batch(batch)
                        records_processed += len(batch)
                        lib.log(f"Parsed {records_processed} variants...")
                        batch = []

                if batch:
                    self.disease_database.insert_mutation_batch(batch)
                    records_processed += len(batch)

            lib.log(f"Database generation completed. Integrated {records_processed} variants.")
            self.disease_database.sys_set_last_update_now()

        except Exception as e:
            lib.log(f"An error occured while parsing DB: {e}")
        finally:
            self.purge_temp()

        self.gui_cmd_buffer.put(["DONE", None])

    def handle_disease_db_update(self):
        """ Handles DB update """
        last_update_str = self.disease_database.get_last_update_date()

        # Check time treshhold
        if last_update_str:
            last_update = datetime.fromisoformat(last_update_str).date()
            if datetime.now().date() - last_update < timedelta(days=7):
                return

        # Update DB
        lib.log("Weekly DB update threshold reached. Initializing update...")
        self.update_disease_database()

        # Save update timestamp
        self.disease_database.sys_set_last_update_now()
