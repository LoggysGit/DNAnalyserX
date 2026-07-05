import queue

import customtkinter as ctk
from tkinter import ttk, filedialog

import modules.lib as lib

class ValidatedNumberEntry(ctk.CTkEntry):
    def __init__(self, master, min_val=0, max_val=100, allow_float=False, default_val=None, **kwargs):
        if default_val is not None and "placeholder_text" not in kwargs: 
            kwargs["placeholder_text"] = str(default_val)
            
        super().__init__(master, **kwargs)
        self.min_val = min_val
        self.max_val = max_val
        self.allow_float = allow_float
        self.default_val = default_val
        
        vcmd = (self.register(self._validate_input), '%P')
        self.configure(validate="key", validatecommand=vcmd)
        self.bind("<FocusOut>", self._on_focus_out)

    def resize_min_val(self, val): self.min_val = val
    def resize_max_val(self, val): self.max_val = val
        
    def _validate_input(self, new_value):
        if new_value == "": return True
        if new_value == "-" and self.min_val < 0: return True
            
        if self.allow_float:
            if new_value.endswith(".") and new_value.count(".") == 1:
                try:
                    base = new_value[:-1]
                    if base == "" or base == "-": 
                        return True
                    return self.min_val <= float(base) <= self.max_val
                except ValueError: 
                    return False
            try: 
                return self.min_val <= float(new_value) <= self.max_val
            except ValueError: 
                return False
        else:
            if not new_value.lstrip('-').isdigit(): 
                return False
            return self.min_val <= int(new_value) <= self.max_val
        
    def _on_focus_out(self, event):
        raw_value = self.get()
        if raw_value == "" or raw_value == "-":
            if self.default_val is not None:
                self.delete(0, "end")
                self.insert(0, str(self.default_val))

    def get_val(self):
        raw_value = self.get()
        if raw_value == "" or raw_value == "-":
            if self.default_val is not None: 
                return self.default_val
            return self.min_val
        return float(raw_value) if self.allow_float else int(raw_value)
    
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
        
        if len(mutation) < 8: return
        self.set_info(mutation)

        self.init_ui()

    def set_info(self, mut):
        self.id = mut[0]
        self.chrn = mut[1]
        self.pos_start = mut[2]
        self.pos_end = mut[3]
        self.ref = mut[4]
        self.alt = mut[5]
        self.vs = mut[6]
        self.sign = mut[7]
        self.dname = mut[8]

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
        self.geometry("500x350")
        self.resizable(False, False)
        
        self.transient(parent)
        self.grab_set()

        self.init_ui()

    def init_ui(self):
        pass

class App(ctk.CTk):
    def __init__(self, gui_cmd_buff, sys_cmd_buff, dman):
        super().__init__()

        self.data_manager = dman
        
        self.command_queue = gui_cmd_buff
        self.system_command_buffer = sys_cmd_buff

        self.title("Gene Analyzer X")
        self.geometry("800x650")

        ctk.set_appearance_mode("dark")

        self.current_data_file_path = None
        self.current_annotation_file_path = None

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
        self.grid_columnconfigure(0, weight=0, minsize=260)
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

        # Circular Info Button moved to the left panel next to title
        self.btn_info = ctk.CTkButton(
            header_frame,
            text="?",
            width=24,
            height=24,
            corner_radius=12,
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors.get("btn_file_fg", "#242424"),
            hover_color=self.ui_colors.get("btn_file_hover", "#3a3a3a"),
            text_color=self.ui_colors["text_main"],
            command=self.show_info_popup
        )
        self.btn_info.pack(side="right", padx=(5, 0))

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

        lbl_anno_section = ctk.CTkLabel(
            left_panel,
            text="ANNOTATION MAP FILE (.gff3 / .gtf)",
            font=("Arial", 10, "bold"),
            text_color=self.ui_colors["text_muted"]
        )
        lbl_anno_section.pack(anchor="w", padx=15, pady=(10, 2))

        self.btn_select_anno = ctk.CTkButton(
            left_panel,
            text="Choose Map File",
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors["btn_file_fg"],
            hover_color=self.ui_colors["btn_file_hover"],
            text_color=self.ui_colors["text_main"],
            height=32,
            command=self.select_annotation_file
        )
        self.btn_select_anno.pack(fill="x", padx=15, pady=5)

        self.lbl_anno_file = ctk.CTkLabel(
            left_panel,
            text="No file selected",
            font=("Arial", 11),
            text_color=self.ui_colors["text_muted"],
            wraplength=230,
            justify="left"
        )
        self.lbl_anno_file.pack(anchor="w", padx=15, pady=(0, 10))

        lbl_gene_section = ctk.CTkLabel(
            left_panel,
            text="TARGET GENE NAME / ID",
            font=("Arial", 10, "bold"),
            text_color=self.ui_colors["text_muted"]
        )
        lbl_gene_section.pack(anchor="w", padx=15, pady=(10, 2))

        self.entry_gene_id = ctk.CTkEntry(
            left_panel,
            placeholder_text="Name/ID",
            font=("Arial", 12),
            fg_color=self.ui_colors["bg_main"],
            border_color=self.ui_colors["border"],
            text_color=self.ui_colors["text_main"],
            height=32
        )
        self.entry_gene_id.pack(fill="x", padx=15, pady=5)

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

        self.progress_bar = ctk.CTkProgressBar(
            right_panel,
            height=8,
            corner_radius=4,
            fg_color=self.ui_colors["bg_panel"],
            progress_color=self.ui_colors["btn_run_fg"]
        )
        self.progress_bar.set(0.0)
        self.progress_bar.pack(fill="x", pady=(0, 12))

        # --- 2. DATA VIEW AREA (TREEVIEW) ---
        table_frame = ctk.CTkFrame(right_panel, fg_color="transparent")
        table_frame.pack(fill="both", expand=True)

        # Treeview setting
        style = ttk.Style()
        style.theme_use("clam") 

        style.configure(
            "Treeview",
            font=("Arial", 6),
            rowheight=24
        )
        style.configure(
            "Treeview.Heading",
            font=("Arial", 8, "bold")
        )

        columns = ("id", "chr", "position", "ref", "alt", "clnvs", "clnsign", "name")
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
        self.tree.heading("position", text="Position")
        self.tree.heading("ref", text="Ref")
        self.tree.heading("alt", text="Alt")
        self.tree.heading("clnvs", text="CLNVS")
        self.tree.heading("clnsign", text="Significance")
        self.tree.heading("name", text="Disease Name")

        self.tree.column("id", width=20, minwidth=20, anchor="center")
        self.tree.column("chr", width=45, minwidth=40, anchor="center")
        self.tree.column("position", width=100, minwidth=100, anchor="w")
        self.tree.column("ref", width=65, minwidth=65, anchor="center")
        self.tree.column("alt", width=65, minwidth=65, anchor="center")
        self.tree.column("clnvs", width=95, minwidth=95, anchor="center")
        self.tree.column("clnsign", width=100, minwidth=100, anchor="center")
        self.tree.column("name", width=220, minwidth=180, anchor="w")

        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.bind("<Double-1>", self.on_tree_double_click)

        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def show_info_popup(self):
        InfoWindow(self, self.ui_colors, lib.APP_INFO_DIR)

    def select_data_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta.gz *.fa.gz")])
        if file_path:
            filename = file_path.split("/")[-1]
            # Update UI
            self.lbl_data_file.configure(text=filename, text_color=self.ui_colors["text_main"])
            # Save path
            self.current_data_file_path = file_path
    def select_annotation_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("GTF files", "*.gff3 *.gtf")])
        if file_path:
            filename = file_path.split("/")[-1]
            # Update UI
            self.lbl_anno_file.configure(text=filename, text_color=self.ui_colors["text_main"])
            # Save path
            self.current_annotation_file_path = file_path

    def on_tree_double_click(self, e):
        selected_item = self.tree.selection()[0]
        if not selected_item: return

        row_values = self.tree.item(selected_item, "values")

        mock_data = [
            row_values[0],
            row_values[1],
            row_values[2],
            str(int(row_values[2]) + len(row_values[3]) - 1),
            row_values[3],
            row_values[4],
            row_values[5],
            row_values[6],
            row_values[7],
        ]

        MutationDetailWindow(self, self.ui_colors, mock_data)

    def run_analysis(self):
        self.tree.delete(*self.tree.get_children())
        if self.current_data_file_path and self.current_annotation_file_path:
            self.btn_analyse.configure(state="disabled")
            self.system_command_buffer.put(("RUN", [self.current_data_file_path, self.current_annotation_file_path, self.entry_gene_id.get()]))

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

    def read_buffer(self):
        try:
            while True:
                command, payload = self.command_queue.get_nowait()

                match command:
                    case "MUTATION":
                        try:
                            i, pos, clnvs, ref, alt, sign, name = payload
                            tag = "missing" if name == "Not found" else "in_database"
                            self.tree.insert("", "end", 
                                values=(i, f"chr{self.entry_chr.get()}", pos, ref, alt, clnvs, sign, name), 
                                tags=(tag)
                            )
                        except Exception as e: lib.log(f"Error parsing mutation: {e}.")

                    case "DB_UPDATE":
                        self.btn_analyse.configure(state="disabled")
                        # Debug

                    case "DONE":
                        self.btn_analyse.configure(state="normal")

                    case _: pass
                        
                self.command_queue.task_done()
        except queue.Empty: pass

        self.after(100, self.read_buffer)