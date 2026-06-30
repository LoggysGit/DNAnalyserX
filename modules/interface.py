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

class App(ctk.CTk):
    def __init__(self, gui_cmd_buff, sys_cmd_buff, dman):
        super().__init__()

        self.data_manager = dman
        
        self.command_queue = gui_cmd_buff
        self.system_command_buffer = sys_cmd_buff

        self.title("DNA Analyzer")
        self.geometry("800x650")

        ctk.set_appearance_mode("dark")

        self.current_file_path = None

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
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=0)

        # Workspace Container (File Info + Treeview)
        workspace_frame = ctk.CTkFrame(self, fg_color="transparent")
        workspace_frame.grid(row=0, column=0, sticky="nsew", padx=15, pady=(15, 5))

        # 1. FILE METADATA SUB-PANEL
        file_info_frame = ctk.CTkFrame(workspace_frame, fg_color="transparent")
        file_info_frame.pack(fill="x", pady=(0, 10))

        metadata_center_container = ctk.CTkFrame(file_info_frame, fg_color="transparent")
        metadata_center_container.pack(anchor="center")

        self.lbl_status_filename = ctk.CTkLabel(
            metadata_center_container, 
            text="Not selected", 
            font=("Arial", 32, "bold"),
            text_color=self.ui_colors["text_main"]
        )
        self.lbl_status_filename.pack(side="left")

        # 3. ANALYSIS OUTPUT CONSOLE
        table_frame = ctk.CTkFrame(workspace_frame, fg_color="transparent")
        table_frame.pack(fill="both", expand=True)
        
        columns = ("id", "chr", "position", "ref", "alt", "clnvs", "clnsign", "name")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings")

        # Set tree styles
        self.tree.tag_configure(
            'in_database', 
            foreground=self.ui_colors["text_main"], 
            background=self.ui_colors["border"]
        )
        self.tree.tag_configure(
            'missing', 
            foreground=self.ui_colors["text_muted"]
        )
        
        # Tree headers
        self.tree.heading("id", text="№")
        self.tree.heading("chr", text="CHR")
        self.tree.heading("position", text="Position")
        self.tree.heading("ref", text="Ref")
        self.tree.heading("alt", text="Alt")
        self.tree.heading("clnvs", text="CLNVS")
        self.tree.heading("clnsign", text="Significance")
        self.tree.heading("name", text="Disease name")
        
        self.tree.column("id", width=5, anchor="center")
        self.tree.column("chr", width=10, anchor="center")
        self.tree.column("position", width=80, anchor="w")
        self.tree.column("ref", width=50, anchor="center")
        self.tree.column("alt", width=50, anchor="center")
        self.tree.column("clnvs", width=70, anchor="center")
        self.tree.column("clnsign", width=80, anchor="center")
        self.tree.column("name", width=130, anchor="center")
        
        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.bind("<Double-1>", self.on_tree_double_click)
        
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Export to VCF Button
        self.btn_export_vcf = ctk.CTkButton(
            table_frame,
            text="Export to VCF",
            width=110,
            height=28,
            font=("Arial", 11, "bold"),
            fg_color=self.ui_colors.get("btn_file_fg", "#242424"),
            hover_color=self.ui_colors.get("btn_file_hover", "#3a3a3a"),
            text_color=self.ui_colors["text_main"],
            corner_radius=6,
            command=self.export_to_vcf
        )
        self.btn_export_vcf.place(relx=1.0, rely=1.0, anchor="se", x=-35, y=-15)

        # BOTTOM AREA: Control Panel
        control_panel = ctk.CTkFrame(
            self, 
            height=70, 
            fg_color=self.ui_colors["bg_panel"],
            border_color=self.ui_colors["border"],
            border_width=1
        )
        control_panel.grid(row=1, column=0, sticky="ew", padx=15, pady=15)
        
        # File selection button
        self.btn_select = ctk.CTkButton(
            control_panel, 
            text="Select File", 
            width=120, 
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors["btn_file_fg"],
            hover_color=self.ui_colors["btn_file_hover"],
            text_color=self.ui_colors["text_main"],
            command=self.select_file
        )
        self.btn_select.pack(side="left", padx=15, pady=15)

        # File path display label
        self.lbl_file = ctk.CTkLabel(
            control_panel, 
            text="No file selected", 
            font=("Arial", 12), 
            text_color=self.ui_colors["text_muted"]
        )
        self.lbl_file.pack(side="left", padx=5, pady=15)

        # Coordinates input frame
        inputs_frame = ctk.CTkFrame(control_panel, fg_color="transparent")
        inputs_frame.pack(side="left", expand=True, pady=15)

        # Chr Input
        lbl_chr = ctk.CTkLabel(inputs_frame, text="Chr:", font=("Arial", 12), text_color=self.ui_colors["text_main"])
        lbl_chr.pack(side="left", padx=2)
        self.entry_chr = ValidatedNumberEntry(
            inputs_frame, 
            min_val=1,
            max_val=23,
            default_val=1,
            width=60, 
            fg_color=self.ui_colors["bg_main"],
            border_color=self.ui_colors["border"],
            text_color=self.ui_colors["text_main"]
        )
        self.entry_chr.pack(side="left", padx=10)

        # Pos Input
        lbl_pos = ctk.CTkLabel(inputs_frame, text="Pos:", font=("Arial", 12), text_color=self.ui_colors["text_main"])
        lbl_pos.pack(side="left", padx=2)
        self.entry_pos = ValidatedNumberEntry(
            inputs_frame,
            min_val=0,
            max_val=1000000000,
            default_val=0,
            width=120,
            fg_color=self.ui_colors["bg_main"],
            border_color=self.ui_colors["border"],
            text_color=self.ui_colors["text_main"]
        )
        self.entry_pos.pack(side="left", padx=10)

        # Analysis execution button
        self.btn_analyse = ctk.CTkButton(
            control_panel, 
            text="Analyse", 
            width=120, 
            font=("Arial", 12, "bold"),
            fg_color=self.ui_colors["btn_run_fg"], 
            hover_color=self.ui_colors["btn_run_hover"],
            text_color=self.ui_colors["text_main"],
            command=self.run_analysis
        )
        self.btn_analyse.pack(side="right", padx=15, pady=15)

    def select_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.fasta.gz *.fa.gz")])
        if file_path:
            filename = file_path.split("/")[-1]
            
            # Update UI
            self.lbl_file.configure(text=filename, text_color=self.ui_colors["text_main"])
            self.lbl_status_filename.configure(text=filename)
            
            # Save path
            self.current_file_path = file_path

    def on_tree_double_click(self, event):
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
        if self.current_file_path:
            self.btn_analyse.configure(state="disabled")
            self.system_command_buffer.put(("RUN", [self.current_file_path, self.entry_chr.get(), self.entry_pos.get()]))

    def export_to_vcf(self):
        tree_items = self.tree.get_children()
        if not tree_items: return

        default_filename = f"mutations_{self.entry_chr.get()}_{self.entry_pos.get()}_export.vcf"
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