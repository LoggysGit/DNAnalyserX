import os

from pathlib import Path
from datetime import datetime as dt
import gzip

import modules.lib as lib

ENV = "win" if os.name == "nt" else "linux"

DATA_DIR = Path("data")
TEMP_DIR = DATA_DIR / "temp"
LOGS_FILE_DIR = "logs.log"

DB_PATH = DATA_DIR / "disease_database.db"

MAX_NUCL_LENGTH = 50000

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2

GENOME_VER = "GRCh38"
CHROMOSOME_BASE_URL = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_{GENOME_VER}.p14/GCF_000001405.40_{GENOME_VER}.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
DISEASE_BASE_URL = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

def log(data):
    timestamp = dt.now().strftime("%d.%m.%Y-%H:%M:%S:%f")[:-3]
    with open(LOGS_FILE_DIR, 'a', encoding="utf-8") as f: f.write(f"{timestamp} {data}\n")
    print(data)

def open_file(path):
    try:
        lib.log(f"Opening file: {path}")
        with open(path, 'r', encoding="utf-8") as f: return f.readlines()
    except Exception as e: 
        lib.log(f"Error opening {path} file: {e}")
        return []

def gzip_open(path):
    try:
        lib.log(f"Opening GZ file: {path}")
        with gzip.open(path, "rt", encoding="utf-8") as f: return f.readlines()
    except Exception as e:
        lib.log(f"Error opening GZ {path} file: {e}")
        return []