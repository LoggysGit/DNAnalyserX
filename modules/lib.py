import os

from pathlib import Path
from datetime import datetime as dt

import configparser
import gzip

import modules.lib as lib

IS_DEBUG = False

# --- Paths --- #
ENV = "win" if os.name == "nt" else "linux"

DATA_DIR = Path("data")
TEMP_DIR = DATA_DIR / "temp"

LOGS_FILE_DIR = "logs.log"
CONFIG_DIR = "config.cfg"

DB_PATH = DATA_DIR / "disease_database.db"

# --- Constants --- #
config = configparser.ConfigParser()
config.read(CONFIG_DIR)

GENOME_VER = config.get('GENERAL', 'GENOME_VER', fallback='GRCh38')
CHROMOSOME_BASE_URL = config.get('GENERAL', 'CHROMOSOME_BASE_URL')
DISEASE_BASE_URL = config.get('GENERAL', 'DISEASE_BASE_URL')

MAX_NUCL_LENGTH = config.getint('LIMITS', 'MAX_NUCLEOTIDE_LENGTH', fallback=100000)
MAX_INDEL_SIZE = config.getint('LIMITS', 'MAX_INDEL_SIZE', fallback=1000)
START_POS_PADDING = config.getint('LIMITS', 'START_NUCL_POS_PADDING', fallback=0)

MATCH_SCORE = config.getint('SCORING', 'SW_MATCH_SCORE', fallback=1)
MISMATCH_SCORE = config.getint('SCORING', 'SW_MISMATCH_SCORE', fallback=-3)
GAP_OPEN_SCORE = config.getint('SCORING', 'SW_GAP_OPEN_SCORE', fallback=-5)
GAP_EXT_SCORE = config.getint('SCORING', 'SW_GAP_EXTEND_SCORE', fallback=-2)

# --- Functions --- #

def log(data):
    timestamp = dt.now().strftime("%d.%m.%Y-%H:%M:%S:%f")[:-3]
    with open(LOGS_FILE_DIR, 'a', encoding="utf-8") as f: f.write(f"{timestamp} | {data}\n")
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
        with gzip.open(path, "rt", encoding="utf-8") as f: lines = f.readlines()
            
        if not lines: return ""
            
        sequence_lines = [line.strip() for line in lines if not line.startswith(">")]
        return "".join(sequence_lines)
        
    except Exception as e:
        lib.log(f"Error opening GZ {path} file: {e}")
        return ""
    
def dbg(str):
    if IS_DEBUG: print(str)