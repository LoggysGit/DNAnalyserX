import os
import gzip
import time
import random
import sqlite3
from pathlib import Path
import modules.lib as lib

def generate_smart_test_dataset(db_path, chromosome_id, start_pos, length=2000, max_mutations=5, output_dir="data"):
    """
    Fetches real pathogenic variants from the local ClinVar database,
    constructs a baseline reference sequence, injects these real-world mutations 
    at exact positions with correct Ref/Alt pairs, and saves a compressed patient FASTA.
    """
    chrom_clean = str(chromosome_id).lower().replace("chr", "").strip()
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    patient_file_path = output_path / f"patient_chr{chrom_clean}_smart_mutated.fasta.gz"
    
    # 1. Generate a mock base reference string using the target coordinate scale
    # (Since we have local boundaries, we fill it with predictable mock bases)
    random.seed(int(start_pos + length))
    bases_pool = ['A', 'C', 'T', 'G']
    ref_list = [random.choice(bases_pool) for _ in range(length * 2)]
    
    # 2. Connect to the local database and fetch real mutations in this coordinate window
    end_pos = start_pos + length
    query = """
        SELECT position, ref_allele, alt_allele, mutation_type, clinical_significance
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

    # If the database region is empty, we fall back to a couple of random variations
    if not real_mutations:
        print("No real variants found in this specific block. Injecting simulated variants.")
        real_mutations = [
            (start_pos + 100, 'T', 'A', 'single nucleotide variant', 'Pathogenic_Mock'),
            (start_pos + 250, 'A', 'G', 'single nucleotide variant', 'Likely_Pathogenic_Mock')
        ]

    # Select a specific limited sample size to avoid completely destroying the alignment sequence
    sample_size = min(len(real_mutations), max_mutations)
    selected_mutations = random.sample(real_mutations, sample_size)
    # Sort by position descending to mutate from right to left without breaking left-side indices
    selected_mutations.sort(key=lambda x: x[0], reverse=True)

    # 3. Inject the real-world variants directly into the array structure
    injected_stats = {"SNP": 0, "Deletion": 0, "Insertion": 0}
    
    for mut in selected_mutations:
        pos, db_ref, db_alt, mut_type, significance = mut
        
        # Calculate local index relative to our window start position
        local_idx = pos - start_pos
        if local_idx < 0 or local_idx >= length:
            continue
            
        # Overwrite reference tracking sequence to match ClinVar requirements before mutation
        if db_ref != "na" and db_ref != ".":
            for offset, char in enumerate(db_ref):
                if local_idx + offset < len(ref_list):
                    ref_list[local_idx + offset] = char

        # Process variations deployment logic
        if "single nucleotide variant" in mut_type or len(db_ref) == 1 and len(db_alt) == 1:
            ref_list[local_idx] = db_alt
            injected_stats["SNP"] += 1
            print(f"Injected SNP at pos {pos}: {db_ref} -> {db_alt} ({significance})")
            
        elif "deletion" in mut_type or db_alt == ".":
            del_len = len(db_ref) if db_ref != "na" else 1
            # Remove items from list slice block
            del ref_list[local_idx : local_idx + del_len]
            injected_stats["Deletion"] += 1
            print(f"Injected Deletion at pos {pos}: Removed {db_ref} ({significance})")
            
        elif "insertion" in mut_type or db_ref == ".":
            # Insert the alternative sequence alleles directly into the list slot
            ins_seq = list(db_alt) if db_alt != "na" else ['A']
            ref_list[local_idx:local_idx] = ins_seq
            injected_stats["Insertion"] += 1
            print(f"Injected Insertion at pos {pos}: Added {db_alt} ({significance})")

    final_patient_string = "".join(ref_list)[:length]

    # 4. Save the target file to disk using standard gzip compression encoding
    with gzip.open(patient_file_path, "wt") as f_out:
        f_out.write(f">patient_chr{chrom_clean}_smart position={start_pos} length={length}\n")
        f_out.write(final_patient_string + "\n")

    print("\n" + "="*50)
    print(f"SUCCESS: Smart validation file generated from database.")
    print(f"Output File: {patient_file_path}")
    print(f"Total Injected: {injected_stats['SNP']} SNPs, {injected_stats['Deletion']} Del, {injected_stats['Insertion']} Ins")
    print("="*50)
    
    return patient_file_path

if __name__ == "__main__":
    # Provide the path to your compiled SQLite database file
    DB_FILE = lib.DATA_DIR / "disease_database.db" 
    generate_smart_test_dataset(
        db_path=DB_FILE, 
        chromosome_id="21", 
        start_pos=6100400, 
        length=20000, 
        max_mutations=15
    )