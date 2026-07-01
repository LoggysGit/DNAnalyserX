import numpy as np
from numba import njit
import sqlite3

import modules.lib as lib

@njit(fastmath=True)
def generate_sw_matrix(ref_seq_len, pat_seq_len, ref_seq, patient_seq):
    sw_matrix = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)

    hor_gaps = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)
    ver_gaps = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)

    traceback = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)
    # traceback codes: 0 = stop, 1 = diagonal, 2 = vertical (deletion), 3 = horizontal (insertion)

    for i in range(1, ref_seq_len + 1):
        for j in range(1, pat_seq_len + 1):
            diag_score = lib.MATCH_SCORE if ref_seq[i-1] == patient_seq[j-1] else lib.MISMATCH_SCORE
            diag_val = sw_matrix[i-1][j-1] + diag_score

            hor_gaps[i][j] = max(sw_matrix[i][j-1] + lib.GAP_OPEN_SCORE, hor_gaps[i][j-1] + lib.GAP_EXT_SCORE)
            ver_gaps[i][j] = max(sw_matrix[i-1][j] + lib.GAP_OPEN_SCORE, ver_gaps[i-1][j] + lib.GAP_EXT_SCORE)

            best = 0
            ptr = 0
            if diag_val > best:
                best = diag_val
                ptr = 1
            if ver_gaps[i][j] > best:
                best = ver_gaps[i][j]
                ptr = 2
            if hor_gaps[i][j] > best:
                best = hor_gaps[i][j]
                ptr = 3

            sw_matrix[i][j] = best
            traceback[i][j] = ptr

    return sw_matrix, hor_gaps, ver_gaps, traceback

def sw_backtrack(pos, ref_seq, patient_seq, sw_matrix, hor_gaps, ver_gaps, traceback):
    results = []

    flat_idx = np.argmax(sw_matrix)
    max_row, max_col = np.unravel_index(flat_idx, sw_matrix.shape)

    bi, bj = max_row, max_col
    while bi > 0 and bj > 0 and sw_matrix[bi][bj] > 0:
        ptr = traceback[bi][bj]

        if ptr == 1:  # diagonal
            check_score = lib.MATCH_SCORE if ref_seq[bi-1] == patient_seq[bj-1] else lib.MISMATCH_SCORE
            if check_score != lib.MATCH_SCORE:
                results.append([int(pos + bi - 1), "SNV", ref_seq[bi - 1], patient_seq[bj - 1]])
            bi, bj = bi-1, bj-1

        elif ptr == 2:  # vertical
            while bi > 0 and traceback[bi][bj] == 2 and ver_gaps[bi][bj] == ver_gaps[bi-1][bj] + lib.GAP_EXT_SCORE:
                results.append([int(pos + bi - 1), "Deletion", ref_seq[bi-1], "."])
                bi -= 1
            results.append([int(pos + bi - 1), "Deletion", ref_seq[bi-1], "."])
            bi -= 1

        elif ptr == 3:  # horizontal
            while bj > 0 and traceback[bi][bj] == 3 and hor_gaps[bi][bj] == hor_gaps[bi][bj-1] + lib.GAP_EXT_SCORE:
                results.append([int(pos + bi - 1), "Insertion", ".", patient_seq[bj - 1]])
                bj -= 1
            results.append([int(pos + bi - 1), "Insertion", ".", patient_seq[bj - 1]])
            bj -= 1

        else:
            lib.dbg(f"Backtrack stop at bi={bi}, bj={bj}, ptr={ptr}")
            break

    results.sort(key=lambda r: r[0])
    return results

class Core:
    def __init__(self, gui_cmd_buff, dman):
        self.data_manager = dman
        self.gui_command_buffer = gui_cmd_buff

    def run_comparing(self, patient_data, reference_data, position):
        # pos: 1-based genomic coordinate
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

        # Add padding
        padding_seq = full_ref_str[position - lib.START_POS_PADDING - 1 : position - 1]
        patient_seq = padding_seq + patient_seq
        patient_seq_len = len(patient_seq)

        ref_seq = full_ref_str[position - lib.START_POS_PADDING - 1:]
        ref_seq_len = len(ref_seq)
        max_len = patient_seq_len + lib.MAX_INDEL_SIZE
        if ref_seq_len > max_len:
            ref_seq = ref_seq[:max_len]
            ref_seq_len = len(ref_seq)

        results = self.compare_ref(
            position - lib.START_POS_PADDING,
            patient_seq,
            ref_seq,
            patient_seq_len,
            ref_seq_len
        )

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
        sw_matrix, hor_gaps, ver_gaps, traceback = generate_sw_matrix(ref_seq_len, pat_seq_len, ref_seq, patient_seq)
        lib.log(f"Smith-Waterman matrix build. Extracting...")

        results = sw_backtrack(pos, ref_seq, patient_seq, sw_matrix, hor_gaps, ver_gaps, traceback)
        lib.log(f"Comparsion data extracted. Analyzing algorithm done.")

        # --- Done --- #
        lib.dbg(f"Raw results: {results}")
        return self.format_mutation_results(results, ref_seq, pos)
    
    def format_mutation_results(self, raw_res, ref_seq, window_start_pos):
        if not raw_res: return []
        lib.log(f"Result formatting...")

        res = []
        temp_mut = [raw_res[0][0], raw_res[0][1], raw_res[0][2], raw_res[0][3]]

        for i in range(1, len(raw_res)):
            mutation = raw_res[i]
            stride = len(temp_mut[3]) if temp_mut[1] == "Insertion" else len(temp_mut[2])
            if mutation[0] == temp_mut[0] + stride and mutation[1] == temp_mut[1] and temp_mut[1] != "SNP":
                if mutation[2] != ".": temp_mut[2] += mutation[2]
                if mutation[3] != ".": temp_mut[3] += mutation[3]
            else:
                res.append(temp_mut)
                temp_mut = [mutation[0], mutation[1], mutation[2], mutation[3]]

        res.append(temp_mut)

        # Anchor convention for indels
        for mut in res:
            if mut[1] in ("Deletion", "Insertion"):
                anchor_idx = mut[0] - window_start_pos - 1
                if 0 <= anchor_idx < len(ref_seq):
                    anchor = ref_seq[anchor_idx]
                    mut[0] = mut[0] - 1
                    if mut[1] == "Deletion":
                        mut[2] = anchor + mut[2]
                        mut[3] = anchor
                    else:  # Insertion
                        mut[2] = anchor
                        mut[3] = anchor + mut[3]

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