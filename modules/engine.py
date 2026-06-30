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

def sw_backtrack(pos, ref_seq, patient_seq, sw_matrix, hor_gaps, ver_gaps):
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
        else:
            lib.log(f"Backtrack break at bi={bi}, bj={bj}, sw={sw_matrix[bi][bj]}, ver={ver_gaps[bi][bj]}, hor={hor_gaps[bi][bj]}") 
            break

    return results

# DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG #
import customtkinter as ctk
def show_matrix_window(matrix_data):
    """
    Opens a separate CustomTkinter window displaying a massive 2D matrix (up to 2000x2000)
    using a Canvas with virtual scrolling to prevent UI freezing.
    """
    # Check if the NumPy array is empty or uninitialized
    if matrix_data is None or matrix_data.size == 0:
        return

    rows = len(matrix_data)
    cols = len(matrix_data[0])

    # Cell layout configuration
    cell_width = 100
    cell_height = 30
    header_width = 60
    header_height = 30

    # Total dimensions of the data grid
    total_width = header_width + (cols * cell_width)
    total_height = header_height + (rows * cell_height)

    # Top-level window setup
    top = ctk.CTkToplevel()
    top.title("Matrix Viewer")
    top.geometry("900x600")
    top.after(100, lambda: top.focus_set()) # Fix for window focus in CTkToplevel

    # Layout grid configurations for the container
    top.grid_rowconfigure(0, weight=1)
    top.grid_columnconfigure(0, weight=1)

    # Scrollable Canvas setup
    canvas = ctk.CTkCanvas(top, bg="#2b2b2b", highlightthickness=0)
    canvas.grid(row=0, column=0, sticky="nsew")

    def draw_visible_cells(event=None):
        """
        Clears the canvas and draws only the items currently visible in the viewport.
        """
        canvas.delete("all")

        # Get current viewport coordinates
        canvas_width = canvas.winfo_width()
        canvas_height = canvas.winfo_height()

        x_start = canvas.canvasx(0)
        y_start = canvas.canvasy(0)
        x_end = x_start + canvas_width
        y_end = y_start + canvas_height

        # Calculate visible column indices range
        col_start_idx = max(0, int((x_start - header_width) // cell_width))
        col_end_idx = min(cols, int((x_end - header_width) // cell_width) + 1)

        # Calculate visible row indices range
        row_start_idx = max(0, int((y_start - header_height) // cell_height))
        row_end_idx = min(rows, int((y_end - header_height) // cell_height) + 1)

        # Render visible data cells
        for r in range(row_start_idx, row_end_idx):
            for c in range(col_start_idx, col_end_idx):
                x1 = header_width + (c * cell_width)
                y1 = header_height + (r * cell_height)
                x2 = x1 + cell_width
                y2 = y1 + cell_height

                # Draw cell border and text value
                canvas.create_rectangle(x1, y1, x2, y2, outline="#4a4a4a", fill="#242424")
                canvas.create_text(x1 + 8, y1 + (cell_height // 2), text=str(matrix_data[r][c]), 
                                   anchor="w", fill="#ffffff", font=("Arial", 10))

        # Render floating horizontal header (Columns)
        for c in range(col_start_idx, col_end_idx):
            x1 = header_width + (c * cell_width)
            y1 = y_start if y_start > 0 else 0
            x2 = x1 + cell_width
            y2 = y1 + header_height

            canvas.create_rectangle(x1, y1, x2, y2, outline="#4a4a4a", fill="#1f1f1f")
            canvas.create_text(x1 + cell_width // 2, y1 + (header_height // 2), 
                               text=f"Col {c}", fill="#3b8ed0", font=("Arial", 10, "bold"))

        # Render floating vertical header (Rows)
        for r in range(row_start_idx, row_end_idx):
            x1 = x_start if x_start > 0 else 0
            y1 = header_height + (r * cell_height)
            x2 = x1 + header_width
            y2 = y1 + cell_height

            canvas.create_rectangle(x1, y1, x2, y2, outline="#4a4a4a", fill="#1f1f1f")
            canvas.create_text(x1 + header_width // 2, y1 + (cell_height // 2), 
                               text=f"Row {r}", fill="#3b8ed0", font=("Arial", 10, "bold"))

        # Render top-left corner intersection block
        corner_x = x_start if x_start > 0 else 0
        corner_y = y_start if y_start > 0 else 0
        canvas.create_rectangle(corner_x, corner_y, corner_x + header_width, corner_y + header_height, 
                                outline="#4a4a4a", fill="#1a1a1a")

    # Define custom scroll commands to trigger redraw on movement
    def on_vscroll(*args):
        canvas.yview(*args)
        draw_visible_cells()

    def on_hscroll(*args):
        canvas.xview(*args)
        draw_visible_cells()

    # Configure scrollbars to use custom scroll commands
    v_scroll = ctk.CTkScrollbar(top, orientation="vertical", command=on_vscroll)
    v_scroll.grid(row=0, column=1, sticky="ns")
    
    h_scroll = ctk.CTkScrollbar(top, orientation="horizontal", command=on_hscroll)
    h_scroll.grid(row=1, column=0, sticky="ew")

    # Link canvas updates back to the scrollbars
    canvas.configure(xscrollcommand=h_scroll.set, yscrollcommand=v_scroll.set)
    canvas.configure(scrollregion=(0, 0, total_width, total_height))

    # Bind window resizing/initial display to redraw
    canvas.bind("<Configure>", draw_visible_cells)

# DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG #

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
        sw_matrix, hor_gaps, ver_gaps = generate_sw_matrix(ref_seq_len, pat_seq_len, ref_seq, patient_seq)
        lib.log(f"Smith-Waterman matrix build. Extracting...")
        #show_matrix_window(sw_matrix)
        results = sw_backtrack(pos, ref_seq, patient_seq, sw_matrix, hor_gaps, ver_gaps)
        lib.log(f"Comparsion data extracted. Analyzing algorithm done.")

        # --- Done --- #
        results.reverse()
        lib.log(f"Raw results: {results}")
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