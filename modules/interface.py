import os
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
        if new_value == "": 
            return True
        if new_value == "-" and self.min_val < 0: 
            return True
            
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
        # Grid configuration for the main window
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)  # Top working area (expands)
        self.grid_rowconfigure(1, weight=0)  # Bottom control panel (fixed)

                # =====================================================================
        # TOP AREA: Workspace Container (File Info + Progress + Treeview)
        # =====================================================================
        workspace_frame = ctk.CTkFrame(self, fg_color="transparent")
        workspace_frame.grid(row=0, column=0, sticky="nsew", padx=15, pady=(15, 5))

        # 1. FILE METADATA SUB-PANEL (Centered layout with scaled fonts)
        file_info_frame = ctk.CTkFrame(workspace_frame, fg_color="transparent")
        file_info_frame.pack(fill="x", pady=(0, 10))

        # Center container to group filename and size horizontally in the middle
        metadata_center_container = ctk.CTkFrame(file_info_frame, fg_color="transparent")
        metadata_center_container.pack(anchor="center")

        self.lbl_status_filename = ctk.CTkLabel(
            metadata_center_container, 
            text="Not selected", 
            font=("Arial", 32, "bold"), # Scaled up for prominent display
            text_color=self.ui_colors["text_main"]
        )
        self.lbl_status_filename.pack(side="left")

        self.lbl_status_filesize = ctk.CTkLabel(
            metadata_center_container, 
            text="0.0 MB", 
            font=("Arial", 16), # Proportionally adjusted secondary text
            text_color=self.ui_colors["text_muted"]
        )
        self.lbl_status_filesize.pack(side="left", padx=12, pady=(4, 0))

        # 3. ANALYSIS OUTPUT CONSOLE
        table_frame = ctk.CTkFrame(workspace_frame, fg_color="transparent")
        table_frame.pack(fill="both", expand=True)
        
        columns = ("chr", "position", "ref", "alt")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings")
        
        self.tree.heading("chr", text="CHR")
        self.tree.heading("position", text="Position")
        self.tree.heading("ref", text="Ref")
        self.tree.heading("alt", text="Alt")
        
        self.tree.column("chr", width=80, anchor="center")
        self.tree.column("position", width=200, anchor="w")
        self.tree.column("ref", width=100, anchor="center")
        self.tree.column("alt", width=100, anchor="center")
        
        scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

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

        # Central container for coordinates input
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
            
            file_size_bytes = os.path.getsize(file_path)
            file_size_mb = file_size_bytes / (1024 * 1024)
            
            self.lbl_status_filename.configure(text=filename)
            self.lbl_status_filesize.configure(text=f"{file_size_mb:.1f} MB")
            
            # Save path
            self.current_file_path = file_path

    def run_analysis(self):
        if self.current_file_path:
            self.btn_analyse.configure(state="disabled")
            self.system_command_buffer.put(("RUN", [self.current_file_path, self.entry_chr.get(), self.entry_pos.get()]))

    def read_buffer(self):
        try:
            while True:
                command, payload = self.command_queue.get_nowait()

                match command:
                    case "DNA_ANOMALY":
                        try:
                            pos, ref, alt = payload
                            self.tree.insert("", "end", values=(f"chr{self.entry_chr.get()}", pos, ref, alt))
                        except Exception as e: lib.log(f"Error parsing DNA anomaly: {e}.")

                    case "DISEASE":
                        try:
                            p = payload
                        except Exception as e: lib.log(f"Error parsing disease: {e}.")

                    case "DONE":
                        self.btn_analyse.configure(state="normal")

                    case _: pass
                        
                self.command_queue.task_done()
        except queue.Empty: pass

        self.after(100, self.read_buffer)