import os
import re

import time
import queue

import tkinter as tk
from tkinter import ttk, filedialog
import customtkinter as ctk

import modules.lib as lib

class MutationDetailWindow(ctk.CTkToplevel):
    def __init__(self, parent, ui_colors, mutation):
        super().__init__(parent)
        
        self.ui_colors = ui_colors
        
        self.title("Mutation Detailed Information")
        self.geometry("420x460")
        self.resizable(False, False)
        self.configure(fg_color=self.ui_colors.get("bg_main", "#1a1a1a"))
        
        self.transient(parent)
        self.grab_set()
        
        if len(mutation) < 10: return
        self.set_info(mutation)

        self.init_ui()

    def set_info(self, mut):
        self.id =           mut[0]
        self.chrn =         mut[1]
        self.gene =         mut[2]
        self.pos_start =    mut[3]
        self.pos_end =      mut[4]
        self.ref =          mut[5]
        self.alt =          mut[6]
        self.vs =           mut[7]
        self.sign =         mut[8]
        self.dname =        mut[9]

    def init_ui(self):
        main_container = ctk.CTkFrame(self, fg_color="transparent")
        main_container.pack(fill="both", expand=True, padx=20, pady=20)
        
        # --- 1. TOP METADATA BLOCK ---
        self.lbl_clnvs_idx = ctk.CTkLabel(
            main_container,
            text=f"{self.vs} № {self.id}",
            font=("Arial", 28, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        self.lbl_clnvs_idx.pack(anchor="w", pady=(0, 4))
        
        # Associated clinical disease description text line
        self.lbl_disease = ctk.CTkLabel(
            main_container,
            text=f"Disease: {self.dname}",
            font=("Arial", 16),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        self.lbl_disease.pack(anchor="w", pady=(0, 4))
        
        self.lbl_significance = ctk.CTkLabel(
            main_container,
            text=f"Significance: {self.sign}",
            font=("Arial", 16),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        self.lbl_significance.pack(anchor="w", pady=(0, 4))

        self.lbl_chromosome = ctk.CTkLabel(
            main_container,
            text=f"Chromosome : {self.chrn}",
            font=("Arial", 15),
            text_color=self.ui_colors.get("text_muted", "#aaaaaa")
        )
        self.lbl_chromosome.pack(anchor="w", pady=(0, 12))
        
        # --- 2. HORIZONTAL COORDINATES PANEL ---
        lbl_position_title = ctk.CTkLabel(
            main_container,
            text="Position",
            font=("Arial", 13, "bold"),
            text_color=self.ui_colors.get("text_muted", "#aaaaaa")
        )
        lbl_position_title.pack(anchor="w", pady=(0, 2))

        coords_frame = ctk.CTkFrame(
            main_container,
            fg_color=self.ui_colors.get("bg_panel", "#242424"),
            border_color=self.ui_colors.get("border", "#3a3a3a"),
            border_width=1,
            corner_radius=20
        )
        coords_frame.pack(fill="x", pady=(0, 20))
        
        coords_frame.grid_columnconfigure(0, weight=1)
        coords_frame.grid_columnconfigure(1, weight=1)
        
        self.lbl_start_pos = ctk.CTkLabel(
            coords_frame,
            text=f"Start: {self.pos_start}",
            font=("Arial", 15, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        self.lbl_start_pos.grid(row=0, column=0, sticky="w", padx=(20, 10), pady=10)
        
        self.lbl_end_pos = ctk.CTkLabel(
            coords_frame,
            text=f"End: {self.pos_end}",
            font=("Arial", 15, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        self.lbl_end_pos.grid(row=0, column=1, sticky="w", padx=(10, 20), pady=10)
        
        # --- 3. ALLELES SEQUENCES SPLIT VIEW ---
        alleles_container = ctk.CTkFrame(main_container, fg_color="transparent")
        alleles_container.pack(fill="both", expand=True, pady=(0, 15))
        
        alleles_container.grid_columnconfigure(0, weight=1, uniform="alleles")
        alleles_container.grid_columnconfigure(1, weight=1, uniform="alleles")
        alleles_container.grid_rowconfigure(1, weight=1)
        
        lbl_ref_header = ctk.CTkLabel(
            alleles_container,
            text="Ref",
            font=("Arial", 16, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        lbl_ref_header.grid(row=0, column=0, sticky="w", padx=(5, 5), pady=(0, 4))
        
        lbl_alt_header = ctk.CTkLabel(
            alleles_container,
            text="Alt",
            font=("Arial", 16, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff")
        )
        lbl_alt_header.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=(0, 4))
        
        self.txt_ref = ctk.CTkTextbox(
            alleles_container,
            font=("Consolas", 14),
            fg_color=self.ui_colors.get("bg_panel", "#242424"),
            border_color=self.ui_colors.get("border", "#3a3a3a"),
            border_width=1,
            corner_radius=8,
            activate_scrollbars=True
        )
        self.txt_ref.grid(row=1, column=0, sticky="nsew", padx=(0, 5))
        self.txt_ref.insert("1.0", self.ref)
        
        self.txt_alt = ctk.CTkTextbox(
            alleles_container,
            font=("Consolas", 14),
            fg_color=self.ui_colors.get("bg_panel", "#242424"),
            border_color=self.ui_colors.get("border", "#3a3a3a"),
            border_width=1,
            corner_radius=8,
            activate_scrollbars=True
        )
        self.txt_alt.grid(row=1, column=1, sticky="nsew", padx=(5, 0))
        self.txt_alt.insert("1.0", self.alt)
        
        # --- 4. BOTTOM DISMISSAL CONTROLS ---
        btn_close = ctk.CTkButton(
            main_container,
            text="Close",
            height=36,
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors.get("btn_file_fg", "#2d2d2d"),
            hover_color=self.ui_colors.get("btn_file_hover", "#404040"),
            text_color=self.ui_colors.get("text_main", "#ffffff"),
            command=self.destroy
        )
        btn_close.pack(fill="x")

class InfoWindow(ctk.CTkToplevel):
    def __init__(self, parent, ui_colors, md_path):
        super().__init__(parent)

        self.ui_colors = ui_colors
        self.md_file_path = md_path
                
        self.title("Information")
        self.geometry("700x500")
        self.resizable(True, True)
        
        self.transient(parent)
        self.grab_set()

        self.init_ui()
        self.render_md()

    def init_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.textbox = ctk.CTkTextbox(
            self,
            font=("Arial", 13),
            fg_color=self.ui_colors.get("bg_panel", "#1e1e1e"),
            text_color=self.ui_colors.get("text_main", "#ffffff"),
            border_color=self.ui_colors.get("border", "#3a3a3a"),
            border_width=1,
            wrap="word"
        )
        self.textbox.grid(row=0, column=0, sticky="nsew", padx=15, pady=15)

        self.textbox._textbox.tag_configure("H1", font=("Arial", 22, "bold"), foreground=self.ui_colors.get("btn_run_fg", "#ffffff"))
        self.textbox._textbox.tag_configure("H2", font=("Arial", 18, "bold"), foreground=self.ui_colors.get("text_main", "#ffffff"))
        self.textbox._textbox.tag_configure("H3", font=("Arial", 15, "bold"), foreground=self.ui_colors.get("text_main", "#ffffff"))
        self.textbox._textbox.tag_configure("BOLD", font=("Arial", 13, "bold"))
        self.textbox._textbox.tag_configure("BULLET", foreground=self.ui_colors.get("btn_run_fg", "#ffffff"))
        self.textbox._textbox.tag_configure("NORMAL", font=("Arial", 13))

    def render_md(self):
        if not os.path.exists(self.md_file_path):
            self.textbox.insert("end", f"Error: Markdown file not found at path:\n{self.md_file_path}")
            self.textbox.configure(state="disabled")
            return

        with open(self.md_file_path, "r", encoding="utf-8") as file: lines = file.readlines()

        self.textbox.configure(state="normal")
        
        for line in lines:
            stripped_line = line.strip()
            
            # # Header
            if stripped_line.startswith("# "):
                clean_text = stripped_line[2:] + "\n\n"
                self.textbox.insert("end", clean_text, "H1")
                continue
                
            # ## Header
            elif stripped_line.startswith("## "):
                clean_text = stripped_line[3:] + "\n\n"
                self.textbox.insert("end", clean_text, "H2")
                continue

            # ### Header
            elif stripped_line.startswith("### "):
                clean_text = stripped_line[4:] + "\n\n"
                self.textbox.insert("end", clean_text, "H3")
                continue

            # * Item
            elif stripped_line.startswith("* ") or stripped_line.startswith("- "):
                clean_text = stripped_line[2:]
                self.textbox.insert("end", "  • ", "BULLET")
                self.parse_inline_styles(clean_text)
                self.textbox.insert("end", "\n")
                continue
            
            # Regular text
            else:
                if line == "\n": self.textbox.insert("end", "\n")
                else: self.parse_inline_styles(line)

        self.textbox.configure(state="disabled")

    def parse_inline_styles(self, text_line):
        parts = re.split(r"(\*\*.*?\*\*)", text_line)
        for part in parts:
            if part.startswith("**") and part.endswith("**"):
                bold_content = part[2:-2]
                self.textbox.insert("end", bold_content, "BOLD")
            else: self.textbox.insert("end", part, "NORMAL")

class ProgressWindow(ctk.CTkToplevel):
    def __init__(self, parent, ui_colors, task_count, title=""):
        super().__init__(parent)

        self.ui_colors = ui_colors
        self.task_count = task_count
                
        self.title(title)
        self.geometry("380x140")
        self.resizable(False, False)

        self.attributes("-topmost", True)

        self.protocol("WM_DELETE_WINDOW", self.handle_close_attempt)
        
        self.transient(parent)
        self.grab_set()

        self.init_ui()

    def init_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        main_frame = ctk.CTkFrame(self, fg_color="transparent")
        main_frame.grid(row=0, column=0, sticky="nsew", padx=20, pady=20)
        main_frame.grid_columnconfigure(0, weight=1)

        self.lbl_status = ctk.CTkLabel(
            main_frame,
            text="Initializing...",
            font=("Arial", 12, "bold"),
            text_color=self.ui_colors.get("text_main", "#ffffff"),
            anchor="w"
        )
        self.lbl_status.grid(row=0, column=0, sticky="ew", pady=(0, 10))

        self.progress_bar = ctk.CTkProgressBar(
            main_frame,
            orientation="horizontal",
            mode="determinate",
            progress_color=self.ui_colors.get("btn_run_fg", "#1f538d"),
            fg_color=self.ui_colors.get("border", "#3a3a3a")
        )
        self.progress_bar.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        self.progress_bar.set(0.0)

        self.lbl_msg = ctk.CTkLabel(
            main_frame,
            text="",
            font=("Arial", 11),
            text_color=self.ui_colors.get("text_muted", "#8a8a8a"),
            anchor="w"
        )
        self.lbl_msg.grid(row=2, column=0, sticky="ew")

    def update_progress(self, tnum, msg):
        progress_value = min(max(tnum / self.task_count, 0.0), 1.0)
        self.progress_bar.set(progress_value)

        self.lbl_status.configure(text=f"Processing ({tnum}/{self.task_count})")
        self.lbl_msg.configure(text=str(msg))
        
        self.update_idletasks()

    def handle_close_attempt(self): pass

class App(ctk.CTk):
    def __init__(self, gui_cmd_buff, sys_cmd_buff, dman):
        super().__init__()

        self.data_manager = dman
        
        self.command_queue = gui_cmd_buff
        self.system_command_buffer = sys_cmd_buff

        self.title("Gene Analyzer X")
        self.geometry("1000x650")

        ctk.set_appearance_mode("dark")

        self.current_data_file_path = None

        self.analysis_progress_window = None
        self.db_update_progress_window = None

        self.init_colors()
        self.configure(fg_color=self.ui_colors["bg_main"])

        self.init_ui()
        
        self.read_buffer()

    def init_colors(self):
        # Tk/CTk color structure
        self.ui_colors = {
            "bg_main": "#141416",       # Main window
            "bg_panel": "#1d1d22",      # Control panel
            "text_main": "#e4e4e7",     # Labels
            "text_muted": "#71717a",    # Secondary info
            "border": "#2e2e33",        # Subtle border
            
            # Button states
            "btn_file_fg": "#2b2b36",
            "btn_file_hover": "#3f3f4e",
            "btn_run_fg": "#16a34a",
            "btn_run_hover": "#15803d",
            
            # Table specific colors
            "table_bg": "#18181b",
            "table_header": "#27272a"
        }
        # Apply style overrides
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("Treeview", 
                        background=self.ui_colors["table_bg"], 
                        foreground=self.ui_colors["text_main"], 
                        fieldbackground=self.ui_colors["table_bg"], 
                        rowheight=26,
                        bd=0)
        style.configure("Treeview.Heading", 
                        background=self.ui_colors["table_header"], 
                        foreground=self.ui_colors["text_main"], 
                        font=("Arial", 10, "bold"),
                        bd=0)

    def init_ui(self):
        self.grid_columnconfigure(0, weight=0, minsize=230)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # --- LEFT PANEL ---
        left_panel = ctk.CTkFrame(
            self,
            fg_color=self.ui_colors["bg_panel"],
            border_color=self.ui_colors["border"],
            border_width=1,
            corner_radius=0
        )
        left_panel.grid(row=0, column=0, sticky="nsew", padx=(0, 5), pady=0)
        
        header_frame = ctk.CTkFrame(left_panel, fg_color="transparent")
        header_frame.pack(fill="x", padx=15, pady=(20, 15))

        lbl_panel_title = ctk.CTkLabel(
            header_frame,
            text="GeneAnalyzerX",
            font=("Arial", 20, "bold"),
            text_color=self.ui_colors["text_main"]
        )
        lbl_panel_title.pack(side="left")

        # Circular Info Button
        self.btn_info = ctk.CTkButton(
            header_frame,
            text="i",
            width=20,
            height=20,
            corner_radius=12,
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors.get("btn_file_fg", "#242424"),
            hover_color=self.ui_colors.get("btn_file_hover", "#3a3a3a"),
            text_color=self.ui_colors["text_main"],
            command=self.show_info_popup
        )
        self.btn_info.pack(side="right", padx=(5, 0))

        # --- PATIENT DATA FILE SECTION ---
        lbl_data_section = ctk.CTkLabel(
            left_panel,
            text="PATIENT DATA FILE (.fasta.gz)",
            font=("Arial", 10, "bold"),
            text_color=self.ui_colors["text_muted"]
        )
        lbl_data_section.pack(anchor="w", padx=15, pady=(10, 2))

        self.btn_select_data = ctk.CTkButton(
            left_panel,
            text="Choose Data File",
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors["btn_file_fg"],
            hover_color=self.ui_colors["btn_file_hover"],
            text_color=self.ui_colors["text_main"],
            height=32,
            command=self.select_data_file
        )
        self.btn_select_data.pack(fill="x", padx=15, pady=5)

        self.lbl_data_file = ctk.CTkLabel(
            left_panel,
            text="No file selected",
            font=("Arial", 11),
            text_color=self.ui_colors["text_muted"],
            wraplength=230,
            justify="left"
        )
        self.lbl_data_file.pack(anchor="w", padx=15, pady=(0, 10))

        # --- TARGET GENE NAME SECTION ---
        lbl_gene_section = ctk.CTkLabel(
            left_panel,
            text="TARGET GENE NAME",
            font=("Arial", 10, "bold"),
            text_color=self.ui_colors["text_muted"]
        )
        lbl_gene_section.pack(anchor="w", padx=15, pady=(10, 2))

        self.gene_var = tk.StringVar()
        self.gene_var.trace_add("write", self.gene_input_validate)

        self.entry_gene_id = ctk.CTkEntry(
            left_panel,
            placeholder_text="Name",
            font=("Arial", 12),
            textvariable=self.gene_var,
            fg_color=self.ui_colors["bg_main"],
            border_color=self.ui_colors["border"],
            text_color=self.ui_colors["text_main"],
            height=32
        )
        self.entry_gene_id.pack(fill="x", padx=15, pady=5)

        # --- CLEAR ALL ACTION BUTTON ---
        self.btn_clear_all = ctk.CTkButton(
            left_panel,
            text="Clear All Data",
            font=("Arial", 11, "bold"),
            fg_color="transparent",
            hover_color=self.ui_colors.get("btn_file_hover", "#3a3a3a"),
            text_color=self.ui_colors["text_muted"],
            border_width=1,
            border_color=self.ui_colors["border"],
            height=28,
            command=self.clear_all_inputs
        )
        self.btn_clear_all.pack(fill="x", padx=15, pady=(5, 10))

        spacer = ctk.CTkLabel(left_panel, text=" ")
        spacer.pack(fill="both", expand=True)

        self.btn_analyse = ctk.CTkButton(
            left_panel,
            text="Run Gene Analysis",
            font=("Arial", 13, "bold"),
            fg_color=self.ui_colors["btn_run_fg"],
            hover_color=self.ui_colors["btn_run_hover"],
            text_color=self.ui_colors["text_main"],
            height=40,
            command=self.run_analysis
        )
        self.btn_analyse.pack(fill="x", padx=15, pady=(0, 10))

        self.btn_export_vcf = ctk.CTkButton(
            left_panel,
            text="Export Results to VCF",
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors.get("btn_file_fg", "#242424"),
            hover_color=self.ui_colors.get("btn_file_hover", "#3a3a3a"),
            text_color=self.ui_colors["text_muted"],
            height=35,
            state="disabled",
            command=self.export_to_vcf
        )
        self.btn_export_vcf.pack(fill="x", padx=15, pady=(0, 20))

        # --- RIGHT PANEL ---
        right_panel = ctk.CTkFrame(self, fg_color="transparent")
        right_panel.grid(row=0, column=1, sticky="nsew", padx=15, pady=15)

        table_frame = ctk.CTkFrame(right_panel, fg_color="transparent")
        table_frame.pack(fill="both", expand=True)

        style = ttk.Style()
        style.theme_use("clam") 

        style.configure(
            "Treeview",
            font=("Arial", 8),
            rowheight=24
        )
        style.configure( "Treeview.Heading",
            font=("Arial", 10, "bold")
        )

        columns = ("id", "chr", "gene", "position", "ref", "alt", "clnvs", "clnsign", "name")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings")
        
        self.tree.tag_configure(
            'in_database',
            foreground=self.ui_colors["text_main"],
            background=self.ui_colors["border"]
        )
        self.tree.tag_configure(
            'missing',
            foreground=self.ui_colors["text_muted"]
        )

        self.tree.heading("id", text="№")
        self.tree.heading("chr", text="CHR")
        self.tree.heading("gene", text="Gene")
        self.tree.heading("position", text="Position")
        self.tree.heading("ref", text="Ref")
        self.tree.heading("alt", text="Alt")
        self.tree.heading("clnvs", text="CLNVS")
        self.tree.heading("clnsign", text="Significance")
        self.tree.heading("name", text="Disease Name")

        self.tree.column("id", width=20, minwidth=20, anchor="center")
        self.tree.column("chr", width=35, minwidth=35, anchor="center")
        self.tree.column("gene", width=55, minwidth=45, anchor="center")
        self.tree.column("position", width=110, minwidth=100, anchor="w")
        self.tree.column("ref", width=75, minwidth=65, anchor="center")
        self.tree.column("alt", width=75, minwidth=65, anchor="center")
        self.tree.column("clnvs", width=105, minwidth=95, anchor="center")
        self.tree.column("clnsign", width=120, minwidth=100, anchor="center")
        self.tree.column("name", width=240, minwidth=180, anchor="w")

        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.bind("<Double-1>", self.on_tree_double_click)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
    
    # === Event functions === #
    def show_info_popup(self): InfoWindow(self, self.ui_colors, lib.APP_INFO_DIR)

    def gene_input_validate(self, *args):
        current_text = self.gene_var.get()

        if any(char.islower() for char in current_text):
            cursor_pos = self.entry_gene_id._entry.index(tk.INSERT)
            self.gene_var.set(current_text.upper())
            self.entry_gene_id._entry.icursor(cursor_pos)

    def select_data_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta.gz *.fa.gz")])
        if file_path:
            filename = file_path.split("/")[-1]
            # Update UI
            self.lbl_data_file.configure(text=filename, text_color=self.ui_colors["text_main"])
            # Save path
            self.current_data_file_path = file_path

    def clear_all_inputs(self):
        self.lbl_data_file.configure(text="No file selected")
        self.entry_gene_id.delete(0, "end")

        for item in self.tree.get_children(): self.tree.delete(item)
            
        self.btn_export_vcf.configure(state="disabled")

    def on_tree_double_click(self, e):
        selected_item = self.tree.selection()[0]
        if not selected_item: return

        row_values = self.tree.item(selected_item, "values")

        mock_data = [
            row_values[0],
            row_values[1],
            row_values[2],
            row_values[3],
            str(int(row_values[3]) + len(row_values[5]) - 1),
            row_values[4],
            row_values[5],
            row_values[6],
            row_values[7],
            row_values[8]
        ]

        MutationDetailWindow(self, self.ui_colors, mock_data)

    def run_analysis(self):
        gene_id = self.entry_gene_id.get()
        if self.current_data_file_path and gene_id.strip() != "":
            # Interface
            self.btn_analyse.configure(state="disabled")
            self.analysis_progress_window = ProgressWindow(self, self.ui_colors, 7, "Analysis Progress")

            # Send a command
            self.after(50, lambda: self.system_command_buffer.put(("RUN", [self.current_data_file_path, gene_id])))
            
        else: lib.dbg("Cannot start the analysis.")

        self.clear_all_inputs()
            
    def export_to_vcf(self):
        tree_items = self.tree.get_children()
        if not tree_items: return

        default_filename = f"gene_{self.entry_gene_id.get()}_export.vcf"
        file_path = filedialog.asksaveasfilename(
            defaultextension=".vcf",
            filetypes=[("VCF Files", "*.vcf"), ("All Files", "*.*")],
            initialfile=default_filename,
            title="Export Analysis"
        )
        
        if not file_path: return

        parsed_mutations = []
        for item_id in tree_items:
            row_values = self.tree.item(item_id, "values")

            parsed_mutations.append({
                "id": row_values[0],
                "chr": row_values[1],
                "position": row_values[2],
                "ref": row_values[3],
                "alt": row_values[4],
                "clnvs": row_values[5],
                "clnsign": row_values[6],
                "name": row_values[7]
            })

        target_chrom_name = f"chr{self.entry_chr.get()}"
        self.system_command_buffer.put(("EXPORT", [parsed_mutations, target_chrom_name, file_path]))

    # === Command buffer handler === #
    def read_buffer(self):
        try:
            while True:
                command, payload = self.command_queue.get_nowait()

                match command:
                    case "MUTATION":
                        try:
                            i, chr, gene, pos, clnvs, ref, alt, sign, name = payload
                            self.tree.insert("", "end", 
                                values=(i, chr, gene, pos, ref, alt, clnvs, sign, name), 
                                tags=("missing" if name in ["Not found", "Unknown", "None"] else "in_database")
                            )
                        except Exception as e: lib.log(f"Error parsing mutation {payload}: {e}.")

                    case "DB_UPDATE":
                        self.btn_analyse.configure(state="disabled")
                        self.db_update_progress_window = ProgressWindow(self, self.ui_colors, 1, "DB Update Progress")
                        self.db_update_progress_window.update_progress(0, "Updating disease data...")

                    case "PROGRESS":
                        task_num, msg = payload
                        if self.analysis_progress_window: self.analysis_progress_window.update_progress(task_num, msg)

                    case "DONE":
                        self.btn_analyse.configure(state="normal")

                        # Close analysis progressbar
                        if self.analysis_progress_window and self.analysis_progress_window.winfo_exists():
                            self.analysis_progress_window.grab_release()
                            self.analysis_progress_window.destroy()
                            self.analysis_progress_window = None
                        # Close DB update process progressbar
                        if self.db_update_progress_window and self.db_update_progress_window.winfo_exists():
                            self.db_update_progress_window.update_progress(1, "DB update completed.")
                            self.db_update_progress_window.grab_release()
                            self.db_update_progress_window.destroy()
                            self.db_update_progress_window = None

                    case _: pass
                        
                self.command_queue.task_done()
        except queue.Empty: pass

        self.after(100, self.read_buffer)