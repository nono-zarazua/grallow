#!/usr/bin/env python3

import tkinter as tk
from tkinter import messagebox
from tkinterdnd2 import DND_FILES, TkinterDnD
import subprocess
import os
import re
import glob

class LabApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Grallow Local QC Tool")
        self.root.geometry("500x400")

        self.label = tk.Label(root, text="Drag & Drop Batch Folder for QC", 
                              pady=20, font=("Helvetica", 12, "bold"))
        self.label.pack()

        # Fixed relief to "groove" to avoid TclError
        self.drop_zone = tk.Label(root, text="\n\n[ DROP FOLDER HERE ]\n\n", 
                                  bg="#f0f0f0", width=40, height=10,
                                  relief="groove", bd=2)
        self.drop_zone.pack(padx=20, pady=10)

        self.drop_zone.drop_target_register(DND_FILES)
        self.drop_zone.dnd_bind('<<Drop>>', self.process_files)

        self.status = tk.Label(root, text="Ready to run QC", fg="green")
        self.status.pack(pady=10)

    def get_paths(self, event_data):
        """Standardizes paths from TkinterDnD (handles spaces and braces)."""
        # This regex splits paths correctly even if they contain spaces and are wrapped in {}
        return re.findall(r'{(.*?)}|(\S+)', event_data)

    def process_files(self, event):
        dropped_items = []
        for match in self.get_paths(event.data):
            path = match[0] if match[0] else match[1]
            dropped_items.append(path.strip('{}'))

        if not dropped_items: return

        # 1. Identify BAM files
        bam_files = []
        if os.path.isdir(dropped_items[0]):
            import glob
            bam_files = glob.glob(os.path.join(dropped_items[0], "*.bam"))
            batch_name = os.path.basename(dropped_items[0])
        else:
            bam_files = [f for f in dropped_items if f.lower().endswith('.bam')]
            batch_name = "manual_batch"

        if not bam_files:
            messagebox.showerror("Error", "No BAM files found.")
            return

        # 2. GENERATE THE CSV (Matching your required format)
        csv_path = os.path.join(os.path.dirname(__file__), "samples.csv")
        with open(csv_path, 'w') as f:
            f.write("sample_id,bam_path\n")
            for bam in bam_files:
                sample_id = os.path.basename(bam).replace(".sorted.aligned.bam", "").replace(".bam", "")
                f.write(f"{sample_id},{os.path.abspath(bam)}\n")

    def process_files2(self, event):
        dropped_items = []
        for match in self.get_paths(event.data):
            path = match[0] if match[0] else match[1]
            dropped_items.append(path.strip('{}'))

        if not dropped_items: return

        # 1. Identify BAM files
        bam_files = []
        if os.path.isdir(dropped_items[0]):
            import glob
            bam_files = glob.glob(os.path.join(dropped_items[0], "*.bam"))
            batch_name = os.path.basename(dropped_items[0])
        else:
            bam_files = [f for f in dropped_items if f.lower().endswith('.bam')]
            batch_name = "manual_batch"

        if not bam_files:
            messagebox.showerror("Error", "No BAM files found.")
            return

        # 2. GENERATE THE CSV (Matching your required format)
        csv_path = os.path.join(os.path.dirname(__file__), "samples.csv")
        with open(csv_path, 'w') as f:
            f.write("sample_id,bam_path\n")
            for bam in bam_files:
                sample_id = os.path.basename(bam).replace(".sorted.aligned.bam", "").replace(".bam", "")
                f.write(f"{sample_id},{os.path.abspath(bam)}\n")

        # 3. Trigger the Bash script
        self.status.config(text=f"Running QC for {batch_name}...", fg="blue")
        try:
            script_path = os.path.join(os.path.dirname(__file__), 'run_qc.sh')
            # We only need to pass the batch name now; the CSV is already saved
            subprocess.run(['bash', script_path, batch_name], check=True)
            messagebox.showinfo("Success", f"QC Complete for {batch_name}")
        except subprocess.CalledProcessError:
            messagebox.showerror("Error", "Nextflow failed.")
            

if __name__ == "__main__":
    root = TkinterDnD.Tk()
    app = LabApp(root)
    root.mainloop()