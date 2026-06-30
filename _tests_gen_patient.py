import gzip
import random
import sqlite3
from pathlib import Path
import modules.lib as lib
import modules.data_manager as dm

def generate_smart_test_dataset(db_path, chromosome_id, start_pos, reference_sequence, length=2000, max_mutations=5, output_dir="data"):
    chrom_clean = str(chromosome_id).lower().replace("chr", "").strip()
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    patient_file_path = output_path / f"_pat_chr{chrom_clean}_pos{start_pos}.fasta.gz"
    
    end_pos = start_pos + length
    window_reference = reference_sequence[start_pos-1 : end_pos-1]
    
    if len(window_reference) < length:
        print(f"Error: window shorter than requested length.")
        return None
        
    ref_list = list(window_reference)
    
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
        print("CRITICAL: No real clinical variants found. Test canceled.")
        return None

    sample_size = min(len(real_mutations), max_mutations)
    selected_mutations = random.sample(real_mutations, sample_size)
    selected_mutations.sort(key=lambda x: x[0], reverse=True)

    injected_stats = {"SNP": 0, "Length-changing": 0, "Skipped": 0}
    
    for mut in selected_mutations:
        pos, db_ref, db_alt, mut_type, significance, dis_name = mut
        
        # Skip incomplete entries without anchor
        if db_ref == "." or db_alt == "." or db_ref == "na" or db_alt == "na":
            injected_stats["Skipped"] += 1
            print(f"Skipped pos {pos}: incomplete ref/alt ({db_ref} -> {db_alt})")
            continue

        local_idx = pos - start_pos
        ref_len = len(db_ref)
        
        if local_idx < 0 or local_idx + ref_len > length:
            injected_stats["Skipped"] += 1
            continue

        # Verify reference matches before replacement (sanity check)
        actual_ref = "".join(ref_list[local_idx : local_idx + ref_len])
        if actual_ref != db_ref:
            print(f"WARNING pos {pos}: ref mismatch. Expected '{db_ref}', genome has '{actual_ref}'. Skipping.")
            injected_stats["Skipped"] += 1
            continue

        # Universal replacement — works for SNP, Deletion, Duplication, Insertion, Microsatellite, Indel
        ref_list[local_idx : local_idx + ref_len] = list(db_alt)
        
        if len(db_ref) == 1 and len(db_alt) == 1:
            injected_stats["SNP"] += 1
        else:
            injected_stats["Length-changing"] += 1
        
        print(f"Injected {mut_type} at pos {pos}: {db_ref} -> {db_alt} ({significance}-{dis_name})")

    final_patient_string = "".join(ref_list)

    with gzip.open(patient_file_path, "wt") as f_out:
        f_out.write(f">patient_chr{chrom_clean}_smart position={start_pos} length={length}\n")
        f_out.write(final_patient_string + "\n")

    print("\n" + "="*50)
    print(f"SUCCESS: Patient dataset generated.")
    print(f"Output File: {patient_file_path}")
    print(f"Total: {injected_stats['SNP']} SNPs, {injected_stats['Length-changing']} length-changing, {injected_stats['Skipped']} skipped")
    print("="*50)

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
        start_pos=32750830, 
        reference_sequence=FULL_CHROMOSOME_SEQUENCE,
        length=2000, 
        max_mutations=15
    )

    data_man.purge_temp()