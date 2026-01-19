import os
import sys
import subprocess
import threading
import webbrowser
import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from io import StringIO

# --- NOUVEL IMPORT IMPORTANT POUR L'AFFICHAGE GRAPHIQUE ---
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# --- BIOPYTHON IMPORTS ---
from Bio import Entrez, SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment

# =============================================================================
# PART 1: SCIENTIFIC ENGINE (BACKEND)
# =============================================================================

class BioEngine:
    def __init__(self):
        self.full_alignment = None
        self.records_list = []
        self.divergence_dict = {}
        self.final_ref = ""
        self.clean_count_msg = ""

    def load_input_data(self, input_path):
        input_path = input_path.strip('"').strip("'")
        if not os.path.exists(input_path): return None, None
        try:
            if input_path.lower().endswith('.csv'):
                df = pd.read_csv(input_path)
                df.columns = [c.strip() for c in df.columns]
                return df, None
            elif input_path.lower().endswith(('.fa', '.fasta', '.fas')):
                records = list(SeqIO.parse(input_path, "fasta"))
                SeqIO.write(records, "sequences_to_align.fasta", "fasta")
                df = pd.DataFrame({'Common name': [rec.id for rec in records]})
                return df, "sequences_to_align.fasta"
        except Exception as e:
            print(f"Read Error: {e}")
            return None, None
        return None, None

    def fetch_ncbi_sequences(self, df, email, progress_callback=None):
        Entrez.email = email
        sequences_records = []
        
        # Data Cleaning
        initial_count = len(df)
        if 'Common name' in df.columns:
            df['Common name'] = df['Common name'].astype(str).str.strip().str.title()
            df = df.drop_duplicates(subset=['Common name'], keep='first')
        if 'RefSeq Protein accessions' in df.columns:
            df = df.drop_duplicates(subset=['RefSeq Protein accessions'])
            
        self.clean_count_msg = f"Data Cleaned: {initial_count} rows -> {len(df)} unique species kept."
        
        total = len(df)
        name_counts = {}
        count = 0
        
        for _, row in df.iterrows():
            common_name = str(row['Common name'])
            if not common_name or common_name.lower() == 'nan': continue
            
            if 'RefSeq Protein accessions' in row:
                acc = str(row['RefSeq Protein accessions']).split(',')[0].strip()
                if not acc or acc == 'nan': continue
                
                try:
                    handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                    record = SeqIO.read(handle, "fasta")
                    handle.close()
                    
                    if common_name in name_counts:
                        name_counts[common_name] += 1
                        unique_name = f"{common_name}_{name_counts[common_name]}"
                    else:
                        name_counts[common_name] = 1
                        unique_name = common_name
                    
                    record.id = unique_name.replace(" ", "_").replace(":", "")
                    record.description = "" 
                    sequences_records.append(record)
                except:
                    pass
            count += 1
            if progress_callback: progress_callback(count, total, f"Downloading: {common_name}")
                
        outfile = "sequences_to_align.fasta"
        if sequences_records: SeqIO.write(sequences_records, outfile, "fasta")
        return outfile

    def run_muscle_alignment(self, fasta_file, muscle_exe_path):
        output_afa = "alignment_results.afa"
        muscle_exe = muscle_exe_path.strip('"').strip("'")
        if not os.path.exists(muscle_exe): return None
        try:
            subprocess.run([muscle_exe, '-in', fasta_file, '-out', output_afa], check=True, capture_output=True)
            self.full_alignment = AlignIO.read(output_afa, "fasta")
            self.records_list = list(self.full_alignment)
            return self.full_alignment
        except Exception as e:
            print(f"Muscle Error: {e}")
            return None

    def calculate_divergences(self, alignment, ref_species):
        calculator = DistanceCalculator('identity') 
        dm = calculator.get_distance(alignment)
        names_list = [rec.id for rec in alignment]
        
        if ref_species not in names_list:
            matches = [n for n in names_list if ref_species.lower() in n.lower()]
            ref_species = matches[0] if matches else names_list[0]
            
        divergences_dict = {name: dm[ref_species, name] * 100 for name in names_list}
        return names_list, divergences_dict, ref_species

# =============================================================================
# PART 2: GRAPHICAL USER INTERFACE
# =============================================================================

class App:
    def __init__(self, root):
        self.root = root
        self.logic = BioEngine()
        self.root.title("Phylogenetic Analysis Tool (Fixed)")
        self.root.geometry("850x950")
        self.root.configure(bg="#f0f0f0")
        
        self.muscle_path_var = tk.StringVar()
        self.df_data = None
        self.font_header = ("Segoe UI", 12, "bold")

        self.check_muscle_startup()
        self.setup_ui()

    def check_muscle_startup(self):
        if getattr(sys, 'frozen', False):
            app_path = os.path.dirname(sys.executable)
        else:
            app_path = os.path.dirname(os.path.abspath(__file__))
        muscle_loc = os.path.join(app_path, "muscle.exe")
        
        if not os.path.exists(muscle_loc):
            if messagebox.askyesno("Muscle Missing", "Muscle.exe not found.\nDownload it now?"):
                webbrowser.open("https://drive5.com/muscle/downloads_v3.htm")
            self.muscle_path_var.set("")
        else:
            self.muscle_path_var.set(muscle_loc)

    def setup_ui(self):
        # Header
        tk.Label(self.root, text="🧬 Automated Phylogeny & Alignment", 
                 font=("Segoe UI", 16, "bold"), bg="#1565c0", fg="white", pady=15).pack(fill="x")

        # --- 1. CONFIGURATION ---
        self.f_conf = tk.LabelFrame(self.root, text=" 1. Configuration ", font=self.font_header, bg="white", padx=10, pady=10)
        self.f_conf.pack(fill="x", padx=15, pady=10)

        tk.Label(self.f_conf, text="NCBI Email:", bg="white").grid(row=0, column=0, sticky="w")
        self.email_ent = tk.Entry(self.f_conf, width=40)
        self.email_ent.grid(row=0, column=1, padx=5, pady=5)

        tk.Label(self.f_conf, text="Data File:", bg="white").grid(row=1, column=0, sticky="w")
        self.path_ent = tk.Entry(self.f_conf, width=40)
        self.path_ent.grid(row=1, column=1, padx=5, pady=5)
        tk.Button(self.f_conf, text="📂", command=self.browse_data).grid(row=1, column=2)

        tk.Label(self.f_conf, text="Muscle.exe:", bg="white").grid(row=2, column=0, sticky="w")
        self.muscle_ent = tk.Entry(self.f_conf, textvariable=self.muscle_path_var, width=40)
        self.muscle_ent.grid(row=2, column=1, padx=5, pady=5)
        tk.Button(self.f_conf, text="📂", command=self.browse_muscle).grid(row=2, column=2)

        tk.Button(self.f_conf, text="LOAD CONFIGURATION", bg="#bbdefb", command=self.load_config).grid(row=3, columnspan=3, sticky="we", pady=10)

        # --- 2. EXECUTION ---
        self.f_run = tk.LabelFrame(self.root, text=" 2. Execution ", font=self.font_header, bg="white", padx=10, pady=10)
        
        tk.Label(self.f_run, text="Reference Species:", bg="white").pack(anchor="w")
        self.ref_combo = ttk.Combobox(self.f_run, state="readonly", width=50)
        self.ref_combo.pack(fill="x", pady=5)

        self.btn_run = tk.Button(self.f_run, text="RUN ANALYSIS", bg="#e0e0e0", state="disabled",
                                 font=("Segoe UI", 10, "bold"), command=self.start_analysis_thread)
        self.btn_run.pack(fill="x", pady=10)

        self.lbl_status = tk.Label(self.f_run, text="Waiting...", bg="white", fg="gray")
        self.lbl_status.pack(anchor="w")
        self.pg_bar = ttk.Progressbar(self.f_run, orient="horizontal", mode="determinate")
        self.pg_bar.pack(fill="x", pady=5)

        # --- 3. RESULTS ---
        self.f_vis = tk.LabelFrame(self.root, text=" 3. Results & Tree ", font=self.font_header, bg="white", padx=10, pady=10)

        tk.Label(self.f_vis, text="Select Species for the Tree (Ctrl + Click):", bg="white").pack(anchor="w")
        
        # TABLE (Treeview)
        columns = ("species", "diff")
        self.tree = ttk.Treeview(self.f_vis, columns=columns, show="headings", height=12, selectmode="extended")
        self.tree.heading("species", text="Species Name")
        self.tree.heading("diff", text="% Difference (vs Ref)")
        self.tree.column("species", width=350)
        self.tree.column("diff", width=150, anchor="center")
        self.tree.pack(fill="both", expand=True, pady=5)

        sb = tk.Scrollbar(self.tree, command=self.tree.yview)
        sb.pack(side="right", fill="y")
        self.tree.config(yscrollcommand=sb.set)

        btn_frame = tk.Frame(self.f_vis, bg="white")
        btn_frame.pack(fill="x", pady=5)

        # Bouton Arbre
        tk.Button(btn_frame, text="🌳 DRAW TREE", bg="#ff9800", fg="white", font=("Segoe UI", 9, "bold"),
                  command=self.draw_tree_window).pack(side="left", fill="x", expand=True, padx=2)
        
        # Bouton Texte
        tk.Button(btn_frame, text="📄 VIEW ALIGNMENT", bg="#9c27b0", fg="white", font=("Segoe UI", 9, "bold"),
                  command=self.show_text).pack(side="left", fill="x", expand=True, padx=2)

    # --- LOGIC ---

    def browse_data(self):
        f = filedialog.askopenfilename(filetypes=[("Data", "*.csv *.fasta")])
        if f: 
            self.path_ent.delete(0, tk.END)
            self.path_ent.insert(0, f)

    def browse_muscle(self):
        f = filedialog.askopenfilename(filetypes=[("EXE", "*.exe")])
        if f: self.muscle_path_var.set(f)

    def load_config(self):
        path = self.path_ent.get().strip()
        user_email = self.email_ent.get().strip()
        Entrez.email = user_email if user_email else None
        
        df, _ = self.logic.load_input_data(path)
        if df is not None:
            if 'Common name' in df.columns:
                clean_names = sorted(df['Common name'].astype(str).str.strip().str.title().drop_duplicates().tolist())
            else:
                clean_names = []
            
            self.df_data = df
            self.ref_combo['values'] = clean_names
            if clean_names: self.ref_combo.current(0)
            
            self.f_run.pack(fill="x", padx=15, pady=5)
            self.btn_run.config(state="normal", bg="#43a047", fg="white")
            messagebox.showinfo("Success", f"File Loaded!\n{len(df)} rows detected.")
        else:
            messagebox.showerror("Error", "Invalid File.")

    def update_pg(self, cur, tot, msg):
        self.pg_bar['maximum'] = tot
        self.pg_bar['value'] = cur
        self.lbl_status.config(text=f"{msg} ({cur}/{tot})")
        self.root.update_idletasks()

    def start_analysis_thread(self):
        if not self.muscle_path_var.get():
            messagebox.showerror("Error", "Select Muscle.exe first")
            return
        self.btn_run.config(state="disabled", text="PROCESSING...")
        threading.Thread(target=self.run_analysis, daemon=True).start()

    def run_analysis(self):
        try:
            fasta = self.logic.fetch_ncbi_sequences(self.df_data, Entrez.email, self.update_pg)
            if self.logic.clean_count_msg:
                 self.root.after(0, lambda: messagebox.showinfo("Data Cleaning", self.logic.clean_count_msg))

            self.root.after(0, lambda: self.lbl_status.config(text="Aligning with MUSCLE..."))
            self.root.after(0, lambda: self.pg_bar.config(mode='indeterminate'))
            self.root.after(0, self.pg_bar.start)
            
            align = self.logic.run_muscle_alignment(fasta, self.muscle_path_var.get())
            
            self.root.after(0, self.pg_bar.stop)
            self.root.after(0, lambda: self.pg_bar.config(mode='determinate', value=self.pg_bar['maximum']))
            
            if align:
                ref_choisie = self.ref_combo.get().replace(" ", "_").strip()
                _, div_dict, final_ref = self.logic.calculate_divergences(align, ref_choisie)
                self.logic.divergence_dict = div_dict
                self.logic.final_ref = final_ref
                self.root.after(0, self.on_success)
            else:
                self.root.after(0, lambda: messagebox.showerror("Error", "Muscle Failed"))
                self.root.after(0, lambda: self.btn_run.config(state="normal", text="RETRY"))
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Critical Error", str(e)))
            self.root.after(0, lambda: self.btn_run.config(state="normal", text="RETRY"))

    def on_success(self):
        self.lbl_status.config(text="Done! Select species below.")
        self.f_vis.pack(fill="both", expand=True, padx=15, pady=10)
        
        # Clear & Fill Table
        for i in self.tree.get_children():
            self.tree.delete(i)
            
        for r in self.logic.records_list:
            div = self.logic.divergence_dict.get(r.id, 0.0)
            self.tree.insert("", tk.END, iid=r.id, values=(r.id, f"{div:.2f} %"))
        
        messagebox.showinfo("Success", "Alignment Complete!\nSelect species from the table and click 'DRAW TREE' or 'SHOW TEXT'.")

    def get_selected(self):
        selected_ids = self.tree.selection()
        return [r for r in self.logic.records_list if r.id in selected_ids]

    def draw_tree_window(self):
        """Draws the tree INSIDE a new application window (Crash-proof method)"""
        recs = self.get_selected()
        if len(recs) < 3:
            messagebox.showwarning("Info", "Select at least 3 species in the table.")
            return
        
        try:
            # 1. Calculation
            aln_sub = MultipleSeqAlignment(recs)
            calc = DistanceCalculator('identity')
            tree = DistanceTreeConstructor(calc, 'nj').build_tree(aln_sub)
            
            for clade in tree.find_clades():
                if clade.is_terminal():
                    div = self.logic.divergence_dict.get(clade.name, 0)
                    clade.name = f"{clade.name} ({div:.1f}%)"
                    clade.color = 'red' if clade.name.startswith(str(self.logic.final_ref)) else 'blue'
            
            # 2. Window Setup
            top = tk.Toplevel(self.root)
            top.title(f"Phylogenetic Tree (Ref: {self.logic.final_ref})")
            top.geometry("900x700")
            
            # 3. Embedding Matplotlib
            fig = plt.Figure(figsize=(10, 8), dpi=100)
            ax = fig.add_subplot(111)
            Phylo.draw(tree, axes=ax, do_show=False)
            
            canvas = FigureCanvasTkAgg(fig, master=top)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

        except Exception as e:
            messagebox.showerror("Tree Error", f"Could not draw tree:\n{str(e)}")

    def show_text(self):
        recs = self.get_selected()
        if not recs: return
        top = tk.Toplevel(self.root)
        top.title("Alignment Text")
        top.geometry("900x600")
        txt = tk.Text(top, font=("Courier New", 10), wrap="none")
        txt.pack(fill="both", expand=True, side="left")
        sb = tk.Scrollbar(top, command=txt.yview)
        sb.pack(side="right", fill="y")
        txt.config(yscrollcommand=sb.set)
        with StringIO() as h:
            AlignIO.write(MultipleSeqAlignment(recs), h, "clustal")
            txt.insert("1.0", h.getvalue())

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()