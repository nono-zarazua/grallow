#!/usr/bin/env python3

import tkinter as tk
from tkinter import messagebox
import customtkinter as ctk
from tkinterdnd2 import DND_FILES, TkinterDnD
import subprocess
import os
import re
import pandas as pd
import glob

# Enable CustomTkinter to work with the Drag-and-Drop wrapper
class Tk(TkinterDnD.Tk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class LabApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Grallow Local QC Tool")
        self.root.geometry("550x450")

        # UI Styling
        self.label = ctk.CTkLabel(root, text="Grallow Lab QC Launcher", 
                                  font=("Helvetica", 20, "bold"))
        self.label.pack(pady=20)

        self.drop_zone = ctk.CTkLabel(root, text="\n\n[ DROP FOLDER HERE ]\n\n", 
                                      width=400, height=200, fg_color="#333333",
                                      corner_radius=10)
        self.drop_zone.pack(padx=20, pady=20)

        # FIX: Register and Bind using the correct methods
        self.drop_zone.drop_target_register(DND_FILES)
        self.drop_zone.dnd_bind('<<Drop>>', self.process_files)

        self.status = ctk.CTkLabel(root, text="System Ready", text_color="green")
        self.status.pack(pady=10)

    def clean_path(self, path_str):
        """Standardizes paths from TkinterDnD (handles spaces and braces)."""
        # Logic to remove {} often added by Linux/TkinterDnD for paths with spaces
        paths = re.findall(r'{(.*?)}|(\S+)', path_str)
        if not paths:
            return path_str.strip('{}')
        return (paths[0][0] if paths[0][0] else paths[0][1]).strip()

    def process_files(self, event):
        input_path = self.clean_path(event.data)
        
        if not os.path.isdir(input_path):
            messagebox.showerror("Error", "Please drop the Batch FOLDER, not individual files.")
            return

        # 1. Find the Lab CSV (Sample Sheet)
        csv_files = glob.glob(os.path.join(input_path, "*.csv"))
        if not csv_files:
            messagebox.showerror("Missing CSV", "No .csv sample sheet found in the folder!")
            return

        lab_csv = csv_files[0]
        batch_name = os.path.basename(lab_csv).replace(".csv", "")

        # 2. Translate Lab CSV to Nextflow CSV
        try:
            df_lab = pd.read_csv(lab_csv)
            nextflow_data = []

            # Check for the 'alias' column we saw in your sample file
            if 'alias' not in df_lab.columns:
                messagebox.showerror("CSV Error", "CSV must contain an 'alias' column.")
                return

            for _, row in df_lab.iterrows():
                alias = str(row['alias'])
                # Search for any BAM file starting with the alias
                search_pattern = os.path.join(input_path, f"{alias}*.bam")
                found_files = glob.glob(search_pattern)

                if found_files:
                    nextflow_data.append({
                        'sample_id': alias,
                        'bam_path': os.path.abspath(found_files[0])
                    })

            if not nextflow_data:
                messagebox.showerror("Matching Error", "No BAM files in the folder match the 'alias' column in your CSV.")
                return

            # 3. Save the 'trial_samples.csv' Nextflow expects
            df_nextflow = pd.DataFrame(nextflow_data)
            target_csv = os.path.join(os.getcwd(), "trial_samples.csv")
            df_nextflow.to_csv(target_csv, index=False)

            # 4. Trigger the Bash/Nextflow Execution
            self.status.configure(text=f"🚀 Running QC: {batch_name}...", text_color="orange")
            self.root.update()

            script_path = os.path.join(os.getcwd(), 'run_qc.sh')
            # This opens a NEW GNOME terminal (common on Ubuntu/ThinkPads) to run the script
            subprocess.run(['gnome-terminal', '--wait', '--', 'bash', script_path, batch_name], check=True)

            messagebox.showinfo("Success", f"QC Complete for {batch_name}!\nResults in ./results/qc/{batch_name}")
            self.status.configure(text="✅ QC Finished", text_color="green")

        except Exception as e:
            messagebox.showerror("Pipeline Error", f"Failed to process batch:\n{str(e)}")
            self.status.configure(text="❌ Failed", text_color="red")

if __name__ == "__main__":
    # Initialize the special DnD-enabled window
    root = TkinterDnD.Tk() # This is the most direct way in the new version
    app = LabApp(root)
    root.mainloop()