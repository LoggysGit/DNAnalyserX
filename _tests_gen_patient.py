import os
import gzip
import time
import random
import sqlite3
from pathlib import Path
import modules.lib as lib
import modules.data_manager as dm

def generate_smart_test_dataset(db_path, chromosome_id, start_pos, reference_sequence, length=2000, max_mutations=5, output_dir="data"):
    """
    Accepts a full reference genome string, extracts the exact window coordinate segment,
    fetches pathogenic variants from the local ClinVar database, injects them,
    and outputs a compressed patient FASTA file matching the exact target length.
    """
    chrom_clean = str(chromosome_id).lower().replace("chr", "").strip()
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    patient_file_path = output_path / f"_pat_chr{chrom_clean}_pos{start_pos}.fasta.gz"
    
    # CRUCIAL FIX: Slice the full chromosome sequence to extract ONLY the requested window
    # We slice from start_pos to start_pos + length
    end_pos = start_pos + length
    window_reference = reference_sequence[start_pos-1 : end_pos-1]
    
    if len(window_reference) < length:
        print(f"Error: Extracted window reference sequence length ({len(window_reference)}) is shorter than requested length ({length}).")
        return None
        
    # Convert the sliced target window reference to a mutable list
    ref_list = list(window_reference)
    
    # Connect to the local database and fetch real clinical mutations in this coordinate window
    query = """
        SELECT position, ref_allele, alt_allele, mutation_type, clinical_significance, disease_name
        FROM clinvar_mutations
        WHERE chromosome = ? AND position >= ? AND position <= ?
        LIMIT 50;
    """
    
    real_mutations = []
    try:
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query, (chrom_clean, start_pos, end_pos))
            real_mutations = cursor.fetchall()
            print(f"Found {len(real_mutations)} real ClinVar mutations in this coordinate segment.")
    except Exception as e:
        print(f"Database extraction failed: {e}")
        return None

    if not real_mutations:
        print("CRITICAL: No real clinical variants found in this specific database block. Alignment test canceled.")
        return None

    # Limit and sample variants to avoid full sequence structural destruction
    sample_size = min(len(real_mutations), max_mutations)
    selected_mutations = random.sample(real_mutations, sample_size)
    
    # Sort by position descending to mutate from right to left, preserving lower index integrity
    selected_mutations.sort(key=lambda x: x[0], reverse=True)

    injected_stats = {"SNP": 0, "Deletion": 0, "Insertion": 0}
    
    for mut in selected_mutations:
        pos, db_ref, db_alt, mut_type, significance, dis_name = mut
        
        # Calculate local tracking index relative to the window start point
        local_idx = pos - start_pos
        if local_idx < 0 or local_idx >= length:
            continue
            
        # Overwrite slice sequence to match ClinVar expected reference baseline before mutation deployment
        if db_ref != "na" and db_ref != ".":
            for offset, char in enumerate(db_ref):
                if local_idx + offset < len(ref_list):
                    ref_list[local_idx + offset] = char

        # Deploy variations based on structural database categories
        if "single nucleotide variant" in mut_type or (len(db_ref) == 1 and len(db_alt) == 1):
            ref_list[local_idx] = db_alt
            injected_stats["SNP"] += 1
            print(f"Injected SNP at pos {pos}: {db_ref} -> {db_alt} ({significance}-{dis_name})")
            
        elif "deletion" in mut_type or db_alt == ".":
            del_len = len(db_ref) if db_ref != "na" else 1
            del ref_list[local_idx : local_idx + del_len]
            injected_stats["Deletion"] += 1
            print(f"Injected Deletion at pos {pos}: Removed {db_ref} ({significance}-{dis_name})")
            
        elif "insertion" in mut_type or db_ref == ".":
            ins_seq = list(db_alt) if db_alt != "na" else ['A']
            ref_list[local_idx:local_idx] = ins_seq
            injected_stats["Insertion"] += 1
            print(f"Injected Insertion at pos {pos}: Added {db_alt} ({significance}-{dis_name})")

    # Reconstruct final patient nucleotide string structure
    final_patient_string = "".join(ref_list)

    # Save target file with gzip compression encoding standard
    with gzip.open(patient_file_path, "wt") as f_out:
        f_out.write(f">patient_chr{chrom_clean}_smart position={start_pos} length={length}\n")
        f_out.write(final_patient_string + "\n")

    print("\n" + "="*50)
    print(f"SUCCESS: Validated smart dataset generated from real reference genome window.")
    print(f"Output File: {patient_file_path}")
    print(f"Total Injected: {injected_stats['SNP']} SNPs, {injected_stats['Deletion']} Del, {injected_stats['Insertion']} Ins")
    print("="*50)

    # Verify injections
    with gzip.open(patient_file_path, "rt") as f:
        lines = f.readlines()
        seq = "".join(l.strip() for l in lines if not l.startswith(">"))

    print("\n--- INJECTION VERIFICATION ---")
    print(f"Window ref [0:5]: {window_reference[0:5]}")
    print(f"Saved patient[0:5]: {seq[0:5]}")
    for mut in selected_mutations:
        local_idx = mut[0] - start_pos
        if 0 <= local_idx < len(seq):
            status = "OK" if seq[local_idx] == mut[2] else "FAIL"
            print(f"[{status}] pos={mut[0]} local={local_idx} got='{seq[local_idx]}' expected='{mut[2]}'")

    return patient_file_path

data_man = dm.DataManager(lib.DB_PATH, None)
if __name__ == "__main__":
    DB_FILE = lib.DATA_DIR / "disease_database.db"
    
    chrid = "12"
    # Downloads or loads the massive raw chromosome sequence string
    FULL_CHROMOSOME_SEQUENCE = data_man.download_chromosome(chrid)
    print(f"Downloaded file len: {len(FULL_CHROMOSOME_SEQUENCE)}")

    print(f"D-> reference_sequence[210001]: {FULL_CHROMOSOME_SEQUENCE[210001]}")  # что здесь?
    print(f"D-> reference_sequence[210000]: {FULL_CHROMOSOME_SEQUENCE[210000]}")  # что здесь?
    print(f"D-> reference_sequence[209999]: {FULL_CHROMOSOME_SEQUENCE[209999]}")  # и здесь?

    #print(FULL_CHROMOSOME_SEQUENCE[0:140])
    
    generate_smart_test_dataset(
        db_path=DB_FILE, 
        chromosome_id=chrid, 
        start_pos=32879004, 
        reference_sequence=FULL_CHROMOSOME_SEQUENCE,
        length=2000, 
        max_mutations=15
    )

    data_man.purge_temp()