import os

from pathlib import Path
from datetime import datetime as dt

import configparser

IS_DEBUG = True

# --- Paths --- #
ENV = "win" if os.name == "nt" else "linux"

DATA_DIR = Path("data")
TEMP_DIR = DATA_DIR / "temp"

LOGS_FILE_DIR = "logs.log"
CONFIG_DIR = "config.cfg"

APP_INFO_DIR = DATA_DIR / "app_info.md"
GENES_CACHE_DIR = DATA_DIR / "genes"
DB_PATH = DATA_DIR / "database.db"

# --- Constants --- #
config = configparser.ConfigParser()
config.read(CONFIG_DIR)

GENOME_VER = config.get('GENERAL', 'GENOME_VER', fallback='GRCh38')
CHROMOSOME_BASE_URL = config.get('GENERAL', 'CHROMOSOME_BASE_URL')
DISEASE_BASE_URL = config.get('GENERAL', 'DISEASE_BASE_URL')

MAX_NUCL_LENGTH = config.getint('GENERAL', 'MAX_NUCLEOTIDE_LENGTH', fallback=100000)

# --- Functions --- #

def log(data):
    timestamp = dt.now().strftime("%d.%m.%Y-%H:%M:%S:%f")[:-3]
    with open(LOGS_FILE_DIR, 'a', encoding="utf-8") as f: f.write(f"{timestamp} | {data}\n")
    print(data)
    
def dbg(str):
    if IS_DEBUG: print(f"> {str}")