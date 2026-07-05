### THIS TEST GENERATOR WAS FULLY CREATED BY AI ###

import gzip
import os
import random
from pathlib import Path
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_gene_test_dataset(gene_id, num_snps=2, num_dels=1, num_ins=1, output_dir="data"):
    """
    Downloads a reference sequence and its true matching GFF3 mapping layer directly 
    from NCBI Nucleotide database. Saves all artifacts in compressed .gz formats, 
    injects size-agnostic genomic modifications, and recalculates annotation boundaries.
    """
    output_path = Path(output_dir)
    genes_cache_path = output_path / "genes"
    genes_cache_path.mkdir(parents=True, exist_ok=True)
    
    patient_fasta_path = output_path / f"simulated_patient_{gene_id}.fasta.gz"
    patient_gff_path = output_path / f"simulated_patient_{gene_id}.gff3"
    
    ref_fasta_path = genes_cache_path / f"{gene_id}_reference.fasta.gz"
    ref_gff_path = genes_cache_path / f"{gene_id}_reference.gff3"
    Entrez.email = "your.email@example.com"  # Required by NCBI policy
    
    # ==========================================
    # 1. RESOLVE GENE AND FETCH ARTIFACTS
    # ==========================================
    # Search NCBI Nucleotide to get the correct internal structural ID
    try:
        search_handle = Entrez.esearch(db="nucleotide", term=f"{gene_id}[Gene Name] AND RefSeq[Keyword] AND Homo sapiens[Organism]")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if not search_results["IdList"]:
            print(f"CRITICAL: Target gene identifier {gene_id} not resolved on NCBI servers.")
            return None
        
        fetch_id = search_results["IdList"][0]
    except Exception as e:
        print(f"CRITICAL: NCBI search pipeline failed: {e}")
        return None

    # Download Reference FASTA directly into GZIP compressed storage
    if not ref_fasta_path.exists():
        print(f"Fetching and compressing reference FASTA for {gene_id}...")
        try:
            fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_id, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()
            
            with gzip.open(ref_fasta_path, "wt") as out_file:
                out_file.write(fasta_data)
        except Exception as e:
            print(f"CRITICAL: FASTA acquisition sequence interrupted: {e}")
            return None

    # Download genuine matching GFF3 structural mapping data from Nucleotide database
    if not ref_gff_path.exists():
        print(f"Fetching true structural GFF3 track layout for {gene_id}...")
        try:
            # Fetching from 'nucleotide' db guarantees an actual 9-column genomic mapping record
            fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_id, rettype="gff3", retmode="text")
            gff_data = fetch_handle.read()
            fetch_handle.close()
            
            with open(ref_gff_path, "w") as out_file:
                out_file.write(gff_data)
        except Exception as e:
            print(f"CRITICAL: GFF3 layout acquisition failure: {e}")
            return None

    # Read reference data transparently from our compressed cache archive
    with gzip.open(ref_fasta_path, "rt") as ref_fasta:
        ref_record = SeqIO.read(ref_fasta, "fasta")
    
    standard_sequence_id = ref_record.id
    ref_sequence_list = list(str(ref_record.seq))
    nucleotides = ["A", "T", "G", "C"]
    coordinate_shifts = []

    # ==========================================
    # 2. SIMULATE MUTATION INJECTIONS
    # ==========================================
    stats = {"SNP": 0, "Deletion": 0, "Insertion": 0}
    occupied_positions = set()

    # --- Phase A: Inject SNPs ---
    attempts = 0
    while stats["SNP"] < num_snps and attempts < 100:
        attempts += 1
        pos = random.randint(100, len(ref_sequence_list) - 100)
        if pos in occupied_positions:
            continue
        current_base = ref_sequence_list[pos]
        possible_alts = [b for b in nucleotides if b != current_base]
        alt_base = random.choice(possible_alts)
        ref_sequence_list[pos] = alt_base
        occupied_positions.add(pos)
        stats["SNP"] += 1

    # --- Phase B: Inject Deletions (Variable sizes) ---
    deletion_targets = []
    attempts = 0
    while len(deletion_targets) < num_dels and attempts < 100:
        attempts += 1
        del_length = random.randint(1, 5)
        pos = random.randint(100, len(ref_sequence_list) - 100)
        window = set(range(pos, pos + del_length))
        if window.intersection(occupied_positions):
            continue
        deletion_targets.append((pos, del_length))
        occupied_positions.update(window)

    deletion_targets.sort(key=lambda x: x[0], reverse=True)
    for pos, del_length in deletion_targets:
        del ref_sequence_list[pos : pos + del_length]
        coordinate_shifts.append((pos, -del_length))
        stats["Deletion"] += 1

    # --- Phase C: Inject Insertions (Variable sizes) ---
    insertion_targets = []
    attempts = 0
    while len(insertion_targets) < num_ins and attempts < 100:
        attempts += 1
        ins_length = random.randint(1, 5)
        pos = random.randint(100, len(ref_sequence_list) - 100)
        if pos in occupied_positions:
            continue
        insertion_targets.append((pos, ins_length))
        occupied_positions.add(pos)

    insertion_targets.sort(key=lambda x: x[0], reverse=True)
    for pos, ins_length in insertion_targets:
        inserted_chunk = [random.choice(nucleotides) for _ in range(ins_length)]
        ref_sequence_list[pos:pos] = inserted_chunk
        coordinate_shifts.append((pos, ins_length))
        stats["Insertion"] += 1

    # ==========================================
    # 3. SAVE MUTATED GENE AS COMPRESSED FASTA
    # ==========================================
    final_patient_string = "".join(ref_sequence_list)
    with gzip.open(patient_fasta_path, "wt") as f_out:
        f_out.write(f">{standard_sequence_id} simulated_patient_{gene_id} var_lengths\n")
        f_out.write(final_patient_string + "\n")

    # ==========================================
    # 4. RECALCULATE AND WRITE VALID MATCHING GFF
    # ==========================================
    coordinate_shifts.sort(key=lambda x: x[0])
    
    with open(ref_gff_path, "r") as src_gff, open(patient_gff_path, "w") as dst_gff:
        for line in src_gff:
            if line.startswith("#"):
                dst_gff.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                dst_gff.write(line)
                continue
            try:
                start_coord = int(parts[3])
                end_coord = int(parts[4])
                new_start = start_coord
                new_end = end_coord
                
                for mutation_pos, shift_val in coordinate_shifts:
                    if start_coord > mutation_pos:
                        new_start += shift_val
                    if end_coord > mutation_pos:
                        new_end += shift_val
                
                parts[3] = str(max(1, new_start))
                parts[4] = str(max(1, new_end))
                dst_gff.write("\t".join(parts) + "\n")
            except ValueError:
                dst_gff.write(line)

    print(f"SUCCESS: Generated compressed datasets matching clean RefSeq mapping tracking ID: {standard_sequence_id}")
    return patient_fasta_path, patient_gff_path

if __name__ == "__main__":
    generate_gene_test_dataset("TP53", num_snps=3, num_dels=2, num_ins=1, output_dir="data")