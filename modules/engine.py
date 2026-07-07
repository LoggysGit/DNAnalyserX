import os

import sqlite3
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from Bio import Entrez
from Bio.Align import *

import modules.lib as lib

class Core:
    def __init__(self, gui_cmd_buff, dman):
        self.data_manager = dman
        self.gui_command_buffer = gui_cmd_buff

    def run_comparing(self, patient_data_path, gene_id):
        self.gui_command_buffer.put(("PROGRESS", [0, "Opening patient file..."]))

        # 1. Open file
        in_seq_handle = gzip.open(patient_data_path, "rt")
        patient_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        if len(patient_seq_dict) != 1:
            lib.log(f"Expected exactly one sequence in patient file, got {len(patient_seq_dict)}")
            return []

        # 2. Check & Download reference
        Entrez.email = lib.USER_EMAIL

        search_handle = Entrez.esearch(
            db="nucleotide",
            term=f"{gene_id}[Gene Name] AND RefSeq[Keyword] AND Homo sapiens[Organism]"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            lib.log(f"Gene identifier {gene_id} not resolved on NCBI servers.")
            return []

        fetch_id = search_results["IdList"][0]

        ref_file_name = f"{gene_id}_reference.fasta.gz"
        ref_data_path = os.path.join(lib.GENES_CACHE_DIR, ref_file_name)
        
        self.gui_command_buffer.put(("PROGRESS", [1, "Loading a reference..."]))

        if not os.path.exists(ref_data_path):
            lib.log(f"Reference for {gene_id} not found locally. Fetching from NCBI...")

            try:
                search_handle = Entrez.esearch(fetch_id)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    fetch_id = search_results["IdList"][0]
                    fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_id, rettype="fasta", retmode="text")
                    fasta_data = fetch_handle.read()
                    fetch_handle.close()
                    
                    with gzip.open(ref_data_path, "wt") as out_file: out_file.write(fasta_data)

                else: raise FileNotFoundError(f"Gene identifier {gene_id} not resolved on NCBI servers.")

            except Exception as error:
                lib.log(f"Failed to download NCBI reference: {str(error)}")
                return []
            
        else: lib.log("Reference data file found. Skipping download.")

        # 3. Check & Download GFF3
        gff_file_name = f"{gene_id}_reference.gff3"
        ref_anno_path = os.path.join(lib.GENES_CACHE_DIR, gff_file_name)

        self.gui_command_buffer.put(("PROGRESS", [2, "Loading an annotation..."]))
        
        if not os.path.exists(ref_anno_path):
            lib.log(f"Annotation for {gene_id} not found locally. Fetching from NCBI...")
            
            try:
                search_handle = Entrez.esearch(fetch_id)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    fetch_id = search_results["IdList"][0]

                    fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_id, rettype="gff3", retmode="text")
                    gff_data = fetch_handle.read()
                    fetch_handle.close()
                    
                    with open(ref_anno_path, "w") as out_file: out_file.write(gff_data)
                    
                else: raise FileNotFoundError(f"Gene identifier {gene_id} not resolved on NCBI servers.")

            except Exception as error:
                lib.log(f"Failed to download NCBI annotation: {str(error)}")
                return []
            
        else: lib.log("Annotation file found. Skipping download.")

        in_seq_handle = gzip.open(ref_data_path, "rt")
        ref_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        lib.dbg(f"Patient seq dict: {patient_seq_dict}, Ref seq dict: {ref_seq_dict}")

        patient_raw_record = list(patient_seq_dict.values())[0]
        ref_raw_record = list(ref_seq_dict.values())[0]

        # 4. Align sequences
        self.gui_command_buffer.put(("PROGRESS", [3, "Aligning..."]))

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")

        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        dna_alignment = aligner.align(ref_raw_record.seq, patient_raw_record.seq)[0]

        ref_blocks = dna_alignment.aligned[0]
        pat_blocks = dna_alignment.aligned[1]

        lib.dbg(f"DNA alignment score: {dna_alignment.score}")
        lib.dbg(f"Aligned blocks: {len(ref_blocks)}")

        # 5. Apply annotation
        self.gui_command_buffer.put(("PROGRESS", [4, "Applying annotation..."]))

        with open(ref_anno_path, "r") as in_handle:
            gff_records = list(GFF.parse(in_handle))

        def collect_transcripts(gff_records, gene_id):
            candidates = []
            seen_feature_keys = set()
            for rec in gff_records:
                for feature in rec.features:
                    if feature.id == gene_id or gene_id in feature.qualifiers.get("Name", []):
                        key = (rec.id, feature.location.start, feature.location.end)
                        if key in seen_feature_keys: continue
                        seen_feature_keys.add(key)

                        transcripts = [sf for sf in feature.sub_features if sf.type in ("mRNA", "transcript")]
                        if transcripts:
                            for t in transcripts:
                                cds_feats = sorted(
                                    [sf for sf in t.sub_features if sf.type == "CDS"],
                                    key=lambda x: int(x.location.start),
                                    reverse=(feature.location.strand == -1),
                                )
                                if cds_feats: candidates.append((t, cds_feats))
                        else:
                            cds_feats = sorted(
                                [sf for sf in feature.sub_features if sf.type == "CDS"],
                                key=lambda x: int(x.location.start),
                                reverse=(feature.location.strand == -1),
                            )
                            if cds_feats: candidates.append((feature, cds_feats))

            return candidates

        def is_mane_select(transcript_feature):
            tags = transcript_feature.qualifiers.get("tag", [])
            return any("MANE Select" in t or "MANE_Select" in t for t in tags)

        candidates = collect_transcripts(gff_records, gene_id)
        if not candidates:
            lib.log(f"No CDS-bearing transcript found for {gene_id}")
            return []

        mane_candidates = [c for c in candidates if is_mane_select(c[0])]
        if mane_candidates:
            chosen_transcript, ref_cds_feats = mane_candidates[0]
            lib.log(f"Using MANE Select transcript {chosen_transcript.id}")
        else:
            chosen_transcript, ref_cds_feats = max(
                candidates, key=lambda c: sum(len(f.location) for f in c[1])
            )
            lib.log(f"No MANE Select tag, using longest-CDS transcript {chosen_transcript.id}")

        ref_cds_locs = [f.location for f in ref_cds_feats]

        def map_ref_pos_to_patient(pos):
            pos = int(pos)
            for (r_start, r_end), (p_start, p_end) in zip(ref_blocks, pat_blocks):
                r_start, r_end, p_start, p_end = int(r_start), int(r_end), int(p_start), int(p_end)
                if r_start <= pos < r_end: return p_start + (pos - r_start)

            for (r_start, r_end), (p_start, p_end) in zip(ref_blocks, pat_blocks):
                r_start = int(r_start)
                p_start = int(p_start)
                if pos < r_start:
                    lib.log(f"CDS boundary at ref pos {pos} falls inside an indel/gap; snapping to nearest block")
                    return p_start
                
            return int(pat_blocks[-1][1])

        def map_location_to_patient(loc):
            p_start = int(map_ref_pos_to_patient(int(loc.start)))
            p_end = int(map_ref_pos_to_patient(int(loc.end) - 1)) + 1

            return type(loc)(p_start, p_end, strand=loc.strand)

        patient_cds_locs = [map_location_to_patient(loc) for loc in ref_cds_locs]

        patient_cds_parts = [str(loc.extract(patient_raw_record.seq)) for loc in patient_cds_locs]
        ref_cds_parts = [str(loc.extract(ref_raw_record.seq)) for loc in ref_cds_locs]

        patient_dna_seq = Seq("".join(patient_cds_parts))
        ref_dna_seq = Seq("".join(ref_cds_parts))

        # 6. Turn sequence into aminoacids
        self.gui_command_buffer.put(("PROGRESS", [5, "Translating into aminoacids..."]))

        patient_amino_full = str(patient_dna_seq.translate())
        ref_amino_full = str(ref_dna_seq.translate())

        def truncate_at_stop(amino_str):
            stop_pos = amino_str.find("*")
            if stop_pos == -1: return amino_str, None
            return amino_str[:stop_pos], stop_pos

        patient_amino, patient_stop_pos = truncate_at_stop(patient_amino_full)
        ref_amino, ref_stop_pos = truncate_at_stop(ref_amino_full)

        patient_amino = Seq(patient_amino)
        ref_amino = Seq(ref_amino)

        lib.dbg(f"Patient CDS length: {len(patient_dna_seq)} bp, frame remainder: {len(patient_dna_seq) % 3}")
        lib.dbg(f"Ref stop at aa {ref_stop_pos}, patient stop at aa {patient_stop_pos}")

        protein_aligner = PairwiseAligner()
        protein_aligner.mode = "global"
        protein_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

        protein_aligner.open_gap_score = -10
        protein_aligner.extend_gap_score = -0.5

        protein_alignment = protein_aligner.align(ref_amino, patient_amino)[0]

        # 7. Extracting all mutations
        self.gui_command_buffer.put(("PROGRESS", [6, "Searching for mutations..."]))

        mutations = []

        # Nonsense / Stop-loss
        if patient_stop_pos is not None and ref_stop_pos is not None:
            if patient_stop_pos < ref_stop_pos:
                wt_aa = ref_amino_full[patient_stop_pos] if patient_stop_pos < len(ref_amino_full) else "?"
                mutations.append({
                    "pos": patient_stop_pos + 1,
                    "type": "Nonsense",
                    "hgvs": f"p.{wt_aa}{patient_stop_pos + 1}*",
                    "ref": wt_aa,
                    "alt": "*",
                })
            elif patient_stop_pos > ref_stop_pos:
                new_aa = patient_amino_full[ref_stop_pos] if ref_stop_pos < len(patient_amino_full) else "?"
                mutations.append({
                    "pos": ref_stop_pos + 1,
                    "type": "Stop-loss",
                    "hgvs": f"p.*{ref_stop_pos + 1}{new_aa}ext*?",
                    "ref": "*",
                    "alt": new_aa,
                })
        elif patient_stop_pos is None: lib.log("No stop codon found in patient CDS translation — possible read-through/large deletion")

        # Insertion / Deletion / Substuition
        for ref_range, pat_range in zip(protein_alignment.aligned[0], protein_alignment.aligned[1]):
            ref_start, ref_end = ref_range
            pat_start, pat_end = pat_range

            ref_slice = ref_amino[ref_start:ref_end]
            pat_slice = patient_amino[pat_start:pat_end]

            if ref_slice == pat_slice: continue

            # Substitution
            if len(ref_slice) == len(pat_slice):
                for index in range(len(ref_slice)):
                    if ref_slice[index] != pat_slice[index]:
                        pos = ref_start + index + 1
                        mutations.append({
                            "pos": int(pos),
                            "type": "Substitution",
                            "hgvs": f"p.{ref_slice[index]}{pos}{pat_slice[index]}",
                            "ref": str(ref_slice[index]),
                            "alt": str(pat_slice[index]),
                        })

            # Deletion
            elif len(ref_slice) > len(pat_slice):
                start_pos = ref_start + 1
                end_pos = ref_end

                if start_pos == end_pos:
                    hgvs = f"p.{ref_amino[ref_start]}{start_pos}del"
                else:
                    hgvs = (
                        f"p.{ref_amino[ref_start]}{start_pos}_"
                        f"{ref_amino[ref_end - 1]}{end_pos}del"
                    )

                mutations.append({
                    "pos": int(start_pos),
                    "type": "Deletion",
                    "hgvs": hgvs,
                    "ref": None,
                    "alt": None,
                })
            # Insertion
            elif len(ref_slice) < len(pat_slice):
                left_pos = ref_start
                right_pos = ref_start + 1

                left_aa = ref_amino[left_pos - 1] if left_pos - 1 >= 0 else None
                right_aa = ref_amino[right_pos - 1] if right_pos - 1 < len(ref_amino) else None

                if left_aa is not None and right_aa is not None: hgvs = f"p.{left_aa}{left_pos}_{right_aa}{right_pos}ins{str(pat_slice)}"
                elif left_aa is not None and right_aa is None: hgvs = f"p.{left_aa}{left_pos}_Terins{str(pat_slice)}"
                else: hgvs = f"p.0_{right_aa}{right_pos}ins{str(pat_slice)}"

                mutations.append({
                    "pos": int(left_pos),
                    "type": "Insertion",
                    "hgvs": hgvs,
                    "ref": None,
                    "alt": None,
                })

        self.gui_command_buffer.put(("PROGRESS", [7, "Analysis completed."]))

        return mutations
    
    def find_mutations(self, anomalies, gene):
        full_mutations_data = []
        lib.log(f"Parsing {len(anomalies)} mutations...")

        for a in range(len(anomalies)):
            try:
                anomaly = anomalies[a]
                local_pos, mutation_type, hgvs = anomaly["pos"], anomaly["type"], anomaly["hgvs"]
                ref_fetched, alt_fetched = anomaly["ref"], anomaly["alt"]

                #chromosome, position, ref_allele, alt_allele, clnvs, clinical_significance, disease_name 
                mutation = self.data_manager.disease_database.find_mutation(hgvs, gene)
                if mutation:
                    chrid, p, refer, alter, vs, sign, name = mutation

                    # Process values
                    pos : int = p if p != None else local_pos
                    clnvs : str = vs if vs != None else mutation_type

                    ref : str = refer if refer != None else ref_fetched
                    alt : str = alter if alter != None else alt_fetched

                    chrom : int = chrid if chrid != None else "?"

                    full_mutations_data.append([a+1, f"chr{chrom}", gene, hgvs, pos, ref, alt, clnvs, sign, name])

            except IndexError as e:
                lib.log(f"Data format structure error: {e}")
                continue
                
            except sqlite3.Error as e:
                lib.log(f"Database error: {e}")
                full_mutations_data.append([
                    anomaly[0] if len(anomaly) > 0 else 0,
                    anomaly[1] if len(anomaly) > 1 else "Unknown",
                    anomaly[2] if len(anomaly) > 2 else "-",
                    anomaly[3] if len(anomaly) > 3 else "-",
                    "Database Error",
                    "Database Error"
                ])
                
            except Exception as e:
                lib.log(f"Unexpected mutations fetcher error: {e}")
                continue

        return full_mutations_data