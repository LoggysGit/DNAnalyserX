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

        # Open file
        in_seq_handle = gzip.open(patient_data_path, "rt")
        patient_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        # Check & Download reference and GFF3
        ref_file_name = f"{gene_id}_reference.fasta.gz"
        ref_data_path = os.path.join(lib.GENES_CACHE_DIR, ref_file_name)
        
        self.gui_command_buffer.put(("PROGRESS", [1, "Loading a reference..."]))

        if not os.path.exists(ref_data_path):
            lib.log(f"Reference for {gene_id} not found locally. Fetching from NCBI...")

            Entrez.email = "blank.email@mail.com"
            try:
                search_handle = Entrez.esearch(db="nucleotide", term=f"{gene_id}[Gene Name] AND RefSeq[Keyword]")
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

        # Check & Download GFF3
        gff_file_name = f"{gene_id}_reference.gff3"
        ref_anno_path = os.path.join(lib.GENES_CACHE_DIR, gff_file_name)

        self.gui_command_buffer.put(("PROGRESS", [2, "Loading an annotation..."]))
        
        if not os.path.exists(ref_anno_path):
            lib.log(f"Annotation for {gene_id} not found locally. Fetching from NCBI...")
            Entrez.email = "blank.email@mail.com"
            try:
                search_handle = Entrez.esearch(db="nucleotide", term=f"{gene_id}[Gene Name] AND RefSeq[Keyword] AND Homo sapiens[Organism]")
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    fetch_id = search_results["IdList"][0]

                    fetch_handle = Entrez.efetch(db="nucleotide", id=fetch_id, rettype="gff3", retmode="text")
                    gff_data = fetch_handle.read()
                    fetch_handle.close()
                    
                    with open(ref_anno_path, "w") as out_file:
                        out_file.write(gff_data)
                else: raise FileNotFoundError(f"Gene identifier {gene_id} not resolved on NCBI servers.")

            except Exception as error:
                lib.log(f"Failed to download NCBI annotation: {str(error)}")
                return []
            
        else: lib.log("Annotation file found. Skipping download.")

        in_seq_handle = gzip.open(ref_data_path, "rt")
        ref_seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        lib.dbg(f"Patient seq dict: {patient_seq_dict}, Ref seq dict: {ref_seq_dict}")

        # Apply annotation
        self.gui_command_buffer.put(("PROGRESS", [3, "Applying annotation..."]))

        clear_patient_dna_str = ""
        clear_ref_dna_str = ""

        in_handle = open(ref_anno_path)
        for recP in GFF.parse(in_handle, base_dict=patient_seq_dict):
            for feature in recP.features:
                if feature.id == gene_id or gene_id in feature.qualifiers.get("Name", []):
                    for sub_feat in feature.sub_features:
                        if sub_feat.type == "exon":
                            clear_patient_dna_str += str(recP.seq[int(sub_feat.location.start):int(sub_feat.location.end)])
        
        in_handle.seek(0)
        for recR in GFF.parse(in_handle, base_dict=ref_seq_dict):
            for feature in recR.features:
                if feature.id == gene_id or gene_id in feature.qualifiers.get("Name", []):
                    for sub_feat in feature.sub_features:
                        if sub_feat.type == "exon":
                            clear_ref_dna_str += str(recR.seq[int(sub_feat.location.start):int(sub_feat.location.end)])

        in_handle.close()

        # Turn into aminoacids
        self.gui_command_buffer.put(("PROGRESS", [4, "Translating into aminoacids..."]))

        patient_dna_seq = Seq(clear_patient_dna_str)
        ref_dna_seq = Seq(clear_ref_dna_str)

        # Pad sequence lengths (MUTED)
        #if len(patient_proteins) % 3 != 0:  patient_proteins += "N" * (3 - (len(patient_proteins) % 3))
        #if len(ref_proteins) % 3 != 0: ref_proteins += "N" * (3 - (len(ref_proteins) % 3))

        patient_amino = patient_dna_seq.translate()
        ref_amino = ref_dna_seq.translate()

        # Align sequences
        self.gui_command_buffer.put(("PROGRESS", [5, "Aligning sequences..."]))

        aligner = PairwiseAligner()
        aligner.mode = "global"
        protein_alignment = aligner.align(ref_amino, patient_amino)[0]
        
        lib.dbg(protein_alignment)

        # Find all mutations
        self.gui_command_buffer.put(("PROGRESS", [6, "Searching for all objects..."]))

        mutations = []

        for ref_range, pat_range in zip(protein_alignment.aligned[0], protein_alignment.aligned[1]):
            ref_start, ref_end = ref_range
            pat_start, pat_end = pat_range
            
            ref_slice = ref_amino[ref_start:ref_end]
            pat_slice = patient_amino[pat_start:pat_end]
            
            if ref_slice != pat_slice:
                # Substitution
                if len(ref_slice) == len(pat_slice):
                    for index in range(len(ref_slice)):
                        if ref_slice[index] != pat_slice[index]:
                            pos = ref_start + index + 1
                            mutations.append({
                                pos, "Substitution",
                                f"p.{ref_slice[index]}{pos}{pat_slice[index]}"
                            })
                # Deletion
                elif len(ref_slice) > len(pat_slice):
                    pos = ref_start + 1
                    mutations.append({
                        pos, "Deletion",
                        f"p.{ref_slice}del"
                    })
                # Insertion
                elif len(ref_slice) < len(pat_slice):
                    pos = ref_start + 1
                    mutations.append({
                        pos, "Insertion",
                        f"p.{ref_amino[ref_start]}_{ref_amino[ref_start+1]}ins{pat_slice}"
                    })

        self.gui_command_buffer.put(("PROGRESS", [7, "Finishing..."]))

        return mutations
    
    def find_mutations(self, anomalies):
        full_mutations_data = []
        lib.log(f"Parsing {len(anomalies)} mutations...")

        for anomaly in anomalies:
            try:
                position, mutation_type, hgvs = anomaly

                lib.dbg(f"{anomaly}")

            except IndexError as e:
                lib.log(f"Data format structure error: {e}")
                continue
                
            except sqlite3.Error as e:
                lib.log(f"Database core engine error: {e}")
                full_mutations_data.append([
                    anomaly[0] if len(anomaly) > 0 else 0,
                    anomaly[1] if len(anomaly) > 1 else "Unknown",
                    anomaly[2] if len(anomaly) > 2 else "-",
                    anomaly[3] if len(anomaly) > 3 else "-",
                    "Database Error",
                    "Database Error"
                ])
                
            except Exception as e:
                lib.log(f"Unexpected application error: {e}")
                continue

        return full_mutations_data