import numpy as np
from numba import njit
import sqlite3

import modules.lib as lib

@njit(fastmath=True)
def generate_sw_matrix(ref_seq_len, pat_seq_len, ref_seq, patient_seq):
    sw_matrix = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)

    hor_gaps = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)
    ver_gaps = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)

    for i in range(1, ref_seq_len + 1):
        for j in range(1, pat_seq_len + 1):
            diag_score = lib.MATCH_SCORE if ref_seq[i-1] == patient_seq[j-1] else lib.MISMATCH_SCORE
            hor_gaps[i][j] = max(sw_matrix[i][j-1] + lib.GAP_OPEN_SCORE, hor_gaps[i][j-1] + lib.GAP_EXT_SCORE)
            ver_gaps[i][j] = max(sw_matrix[i-1][j] + lib.GAP_OPEN_SCORE, ver_gaps[i-1][j] + lib.GAP_EXT_SCORE)
            sw_matrix[i][j] = max(0, sw_matrix[i-1][j-1] + diag_score, hor_gaps[i][j], ver_gaps[i][j])

    return sw_matrix, hor_gaps, ver_gaps

class Core:
    def __init__(self, gui_cmd_buff, dman):
        self.data_manager = dman
        self.gui_command_buffer = gui_cmd_buff

    def run_comparing(self, patient_data, reference_data, position):
        # pos: 1-based genomic coordinate (FASTA standard / VCF format)
        if position < 1:
            lib.log("Position must be 1-based")
            return []

        # Parse Patient Data
        if isinstance(patient_data, list):
            clean_lines = [line.strip().upper() for line in patient_data if not line.startswith(">")]
            patient_seq = "".join(clean_lines)
        else: patient_seq = patient_data.strip().upper()

        patient_seq_len = len(patient_seq)

        # Hard Check for overflow
        if patient_seq_len > lib.MAX_NUCL_LENGTH:
            lib.log(f"Patient sequence too long ({patient_seq_len} > {lib.MAX_NUCL_LENGTH}). Skipping.")
            return []

        # Parse Reference Data
        if isinstance(reference_data, list):
            clean_ref_lines = [line.strip().upper() for line in reference_data if not line.startswith(">")]
            full_ref_str = "".join(clean_ref_lines)
        else: full_ref_str = reference_data.strip().upper()

        ref_seq = full_ref_str[position - 1:]
        ref_seq_len = len(ref_seq)
        # Cut reference data
        if ref_seq_len > lib.MAX_NUCL_LENGTH:
            ref_seq = ref_seq[:lib.MAX_NUCL_LENGTH]
            ref_seq_len = len(ref_seq)

        # Compare genome and extract VCF data
        results = self.compare_ref(position, patient_seq, ref_seq, patient_seq_len, ref_seq_len)

        return results
    
    def compare_ref(self, pos, patient_seq, ref_seq, pat_seq_len, ref_seq_len):
        # --- Validation --- #
        if not patient_seq or not ref_seq: return []

        lib.log(f"Ref: {ref_seq_len}, Pat: {pat_seq_len}. Start analyzing...")

        ref_seq = ref_seq.replace('N', '\x00')
        patient_seq = patient_seq.replace('N', '\x00')

        if len(ref_seq) != ref_seq_len or len(patient_seq) != pat_seq_len:
            lib.log(f"Sequence length mismatch")
            return []

        valid = set('ACGT\x00')
        if not set(ref_seq).issubset(valid) or not set(patient_seq).issubset(valid):
            lib.log(f"Invalid nucleotide symbols")
            return []

        # --- Algorithm --- #
        sw_matrix, hor_gaps, ver_gaps = generate_sw_matrix(ref_seq_len, pat_seq_len, ref_seq, patient_seq)
        lib.log(f"Smith-Waterman matrix build. Extracting...")

        # Backtrack
        results = []
                
        flat_idx = np.argmax(sw_matrix)
        max_row, max_col = np.unravel_index(flat_idx, sw_matrix.shape)

        bi, bj = max_row, max_col
        while bi > 0 and bj > 0 and sw_matrix[bi][bj] > 0:
            check_score = lib.MATCH_SCORE if ref_seq[bi-1] == patient_seq[bj-1] else lib.MISMATCH_SCORE

            # First priority - Gaps
            if sw_matrix[bi][bj] == ver_gaps[bi][bj]:
                while bi > 0 and ver_gaps[bi][bj] == ver_gaps[bi-1][bj] + lib.GAP_EXT_SCORE:
                    results.append([int(pos + bi - 1), "Deletion", ref_seq[bi-1], "."])
                    bi -= 1
                results.append([int(pos + bi - 1), "Deletion", ref_seq[bi-1], "."])
                bi -= 1
                
            elif sw_matrix[bi][bj] == hor_gaps[bi][bj]:
                while bj > 0 and hor_gaps[bi][bj] == hor_gaps[bi][bj-1] + lib.GAP_EXT_SCORE:
                    results.append([int(pos + bi - 1), "Insertion", ".", patient_seq[bj - 1]])
                    bj -= 1
                results.append([int(pos + bi - 1), "Insertion", ".", patient_seq[bj - 1]])
                bj-= 1
                
            # Second - mismatch
            elif sw_matrix[bi][bj] == sw_matrix[bi-1][bj-1] + check_score:
                if check_score != lib.MATCH_SCORE: results.append([int(pos + bi - 1), "SNP", ref_seq[bi - 1], patient_seq[bj - 1]])
                bi, bj = bi-1, bj-1

            # Endless loop protection
            else: break
        lib.log(f"Comparsion data extracted. Analyzing algorithm done.")

        results.reverse()
        lib.log(f"Raw results: {results}")
        return self.format_mutation_results(results)
    
    def format_mutation_results(self, raw_res):
        if not raw_res: return []
        lib.log(f"Result formatting...")

        res = []
        # Put first element
        temp_mut = [raw_res[0][0], raw_res[0][1], raw_res[0][2], raw_res[0][3]]

        for i in range(1, len(raw_res)):
            mutation = raw_res[i]

            stride = len(temp_mut[3]) if temp_mut[1] == "Insertion" else len(temp_mut[2])
            if mutation[0] == temp_mut[0] + stride and mutation[1] == temp_mut[1]:
                # Don't save gaps (dots)
                if mutation[2] != ".": temp_mut[2] += mutation[2]
                if mutation[3] != ".": temp_mut[3] += mutation[3]
            else:
                # Single -> Multiple
                if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"

                res.append(temp_mut)
                temp_mut = [mutation[0], mutation[1], mutation[2], mutation[3]]

        if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"
        res.append(temp_mut)

        lib.log("Results formatted. Run done.")
        return res
    
    def find_mutations(self, dna_anomalies, chromosome_id):
        full_mutations_data = []
        lib.log(f"Parsing {len(dna_anomalies)} mutations...")

        for anomaly in dna_anomalies:
            try:
                position = int(anomaly[0])
                ref_allele = anomaly[2]
                alt_allele = anomaly[3]
                
                db_result = self.data_manager.disease_database.find_mutation(
                    chromosome_id, position, ref_allele, alt_allele
                )
                
                if db_result: clinical_significance, disease_name = db_result
                else:
                    clinical_significance = "-"
                    disease_name = "Not found"

                full_mutations_data.append([
                    position,
                    anomaly[1],
                    ref_allele,
                    alt_allele,
                    clinical_significance,
                    disease_name
                ])

            except IndexError as e:
                lib.log(f"Data format structure indexing error inside anomaly tuple: {e}")
                continue
                
            except sqlite3.Error as e:
                lib.log(f"Database core engine execution timeout or processing error: {e}")
                full_mutations_data.append([
                    anomaly[0] if len(anomaly) > 0 else 0,
                    anomaly[1] if len(anomaly) > 1 else "Unknown",
                    anomaly[2] if len(anomaly) > 2 else "-",
                    anomaly[3] if len(anomaly) > 3 else "-",
                    "Database Error",
                    "Database Error"
                ])
                
            except Exception as e:
                lib.log(f"Unexpected application execution layer runtime error: {e}")
                continue

        return full_mutations_data