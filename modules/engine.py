import numpy as np
import modules.lib as lib

class Core:
    def __init__(self, gui_cmd_buff, dman):
        self.data_manager = dman
        self.gui_command_buffer = gui_cmd_buff

    def run(self, patient_file, reference_file, position, chromosome):
        # --- Parse data --- #
        # Patient DNA sector
        raw_lines = patient_file.readlines()

        clean_lines = [line.strip().upper() for line in raw_lines if not line.startswith(">")]
        patient_seq = "".join(clean_lines)

        patient_seq_len = len(patient_seq)

        # Reference DNA sector
        reference_file.seek(0)
        
        letters_read = 0
        ref_buffer = []
        target_reached = False
        for line in reference_file:
            line = line.strip().upper()
            # Skip metadata and empty lines
            if not line or line.startswith(">"): continue
            line_len = len(line)
            
            if not target_reached and letters_read + line_len >= position:
                target_reached = True
                offset = position - letters_read - 1
                ref_buffer.append(line[offset:])
            elif target_reached: ref_buffer.append(line)
                
            letters_read += line_len
        ref_seq = "".join(ref_buffer)
        ref_seq_len = len(ref_seq)

        # --- Main check cycle --- #
        results = self.compare_ref(position, patient_seq, ref_seq, patient_seq_len, ref_seq_len)

        return results
    
    def compare_ref(self, pos, patient_seq, ref_seq, pat_seq_len, ref_seq_len):
        sw_matrix = np.zeros((ref_seq_len + 1, pat_seq_len + 1), dtype=np.int32)

        # DP
        for i in range(1, ref_seq_len + 1):
            for j in range(1, pat_seq_len + 1):
                gapped_v = sw_matrix[i-1][j] + lib.GAP_SCORE
                gapped_h = sw_matrix[i][j-1] + lib.GAP_SCORE

                diag_score = 0
                if ref_seq[i-1] == patient_seq[j-1]: diag_score = sw_matrix[i-1][j-1] + lib.MATCH_SCORE
                else: diag_score = sw_matrix[i-1][j-1] + lib.MISMATCH_SCORE

                sw_matrix[i][j] = max(0, diag_score, gapped_v, gapped_h)

        # Backtrack
        flat_idx = np.argmax(sw_matrix)
        max_row, max_col = np.unravel_index(flat_idx, sw_matrix.shape)

        results = []
        bi, bj = max_row, max_col
        while bi > 0 and bj > 0 and sw_matrix[bi][bj] > 0:
            check_score = lib.MATCH_SCORE if ref_seq[bi-1] == patient_seq[bj-1] else lib.MISMATCH_SCORE

            if sw_matrix[bi][bj] == sw_matrix[bi-1][bj-1] + check_score:
                bi, bj = bi-1, bj-1
                if check_score != lib.MATCH_SCORE:
                    results.append([pos + bi, "SNP", ref_seq[bi], patient_seq[bj]])

            elif sw_matrix[bi][bj] == sw_matrix[bi-1][bj] + lib.GAP_SCORE:
                bi, bj = bi-1, bj
                results.append([pos + bi, "Deletion", ref_seq[bi], "."])

            elif sw_matrix[bi][bj] == sw_matrix[bi][bj-1] + lib.GAP_SCORE:
                bi, bj = bi, bj-1
                results.append([pos + bi, "Insertion", ".", patient_seq[bj]])

        results.reverse()
        return self.format_mutation_results(results)
    
    def format_mutation_results(raw_res):
        if not raw_res: return []

        res = []
        # Put first element
        temp_mut = [raw_res[0][0], raw_res[0][1], raw_res[0][2], raw_res[0][3]]

        for i in range(1, len(raw_res)):
            mutation = raw_res[i]

            if mutation[0] == temp_mut[0] + len(temp_mut[2]) and mutation[1] == temp_mut[1]:
                # Don't copy gaps (dots)
                if mutation[2] != ".": temp_mut[2] += mutation[2]
                if mutation[3] != ".": temp_mut[3] += mutation[3]
            else:
                # Single -> Multiple
                if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"

                res.append(temp_mut)

                temp_mut = [mutation[0], mutation[1], mutation[2], mutation[3]]

        if temp_mut[1] == "SNP" and (len(temp_mut[2]) > 1 or len(temp_mut[3]) > 1): temp_mut[1] = "MNP"
        res.append(temp_mut)

        return res
    
    def find_mutations(self, dna_anomalies, chromosome_id):
        diseases = []
        
        return diseases