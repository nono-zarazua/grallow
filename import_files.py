#!/usr/bin/env python3

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as fd
from tkinter import messagebox
import customtkinter as ctk
from tkinterdnd2 import DND_FILES, TkinterDnD
import subprocess
import os
import sys
import re
import pandas as pd
import glob
import time
from selector import SelectionWindow
from guibaseclass import GuiBaseClass


class LabApp(GuiBaseClass):
    def __init__(self, root):
        super().__init__(root)

        self.setAppTitle("Grallow Lab QC Launcher")

        root.minsize(550,450)

        menu_directory = self.getMenu('Directory')
        menu_directory.add_command(label="Open", command=self.open_directory, underline=0)

        menu_help = self.getMenu('Help')

        self.addStatusBar()
        self.message("Select or drag a directory or use Directory > Open")        


        # Place GUI's main frame
        frame = self.getFrame()

        # Add drop zone to main frame
        self.drop_zone = ctk.CTkLabel(frame, text="\n\n[ DROP FOLDER HERE ]\n\n", 
                                      width=400, height=200, fg_color="white",
                                      corner_radius=10)
        self.drop_zone.pack(padx=20, pady=20, side = "top", fill = "both", expand = True)

        self.select_button = ttk.Button(self.frame, text="Select Folder", command=self.open_directory)
        self.select_button.pack(padx=10, pady=10)

        self.submit_button = ttk.Button(self.frame, text="Submit", command=self.process_files)
        self.submit_button.pack(padx=10, pady=10)

        # Register and Bind using the correct methods
        self.drop_zone.drop_target_register(DND_FILES)
        self.drop_zone.dnd_bind('<<Drop>>', self.handle_drop)

        # Add to LabApp.__init__
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Variables to store selected directory
        self.selected_dir = None

    def handle_drop(self, event):
        self.selected_dir = self.clean_path(event.data)
        self.drop_zone.configure(
            text=f"\n\n[ Folder Loaded ]\n{self.selected_dir}\n\n",
            text_color="green"
        )
        self.message(f"Ready to process: {self.selected_dir}")

    def on_closing(self):
        # Stop any subprocesses or loops here if needed
        self.root.destroy()
        os._exit(0)

    def clean_path(self, path_str):
        """Standardizes paths from TkinterDnD (handles spaces and braces)."""
        # Logic to remove {} often added by Linux/TkinterDnD for paths with spaces
        paths = re.findall(r'{(.*?)}|(\S+)', path_str)
        if not paths:
            return path_str.strip('{}')
        return (paths[0][0] if paths[0][0] else paths[0][1]).strip()

    def update_nextflow_config(self, batch_name):
        """Edits nextflow.config to set params.batch to the current run name."""
        config_path = os.path.join(os.getcwd(), "nextflow.config")
        if not os.path.exists(config_path):
            return

        with open(config_path, 'r') as f:
            content = f.read()

        # Regex replaces the value inside quotes for the 'batch' parameter
        new_content = re.sub(r'batch\s*=\s*["\'].*?["\']', f'batch = "{batch_name}"', content)

        with open(config_path, 'w') as f:
            f.write(new_content)
    
    def open_directory(self):
        # Opens a dialog to select a directory
        dir_path = fd.askdirectory(title="Select a Directory")
        if dir_path:
            self.selected_dir = dir_path
            # self.label.config(text=f"Selected Directory: {dir_path}")
            self.message(f"Directory set to: {dir_path}")

            self.drop_zone.configure(
                text=f"\n\n[ FOLDER LOADED ]\n{self.selected_dir}\n\n", 
                text_color="lightgreen"
            )
        else:
            self.message("No directory selected")

    def process_files(self):
        if not self.selected_dir:
            messagebox.showwarning("Warning", "Please drop a folder or use the menu to select one first.")
            return

        self.root.focus_force()
        self.root.after(100, lambda: self.execute_pipeline(self.selected_dir))
    
    def execute_pipeline(self, data):
        input_path = self.clean_path(data)
        
        if not os.path.isdir(input_path):
            messagebox.showerror("Error", "Please drop the Batch FOLDER, not individual files.")
            return

        # 1. Identify Batch Name from CSV (corrida_YYYYMMDD_#)
        csv_files = glob.glob(os.path.join(input_path, "*.csv"))
        corrida_pattern = re.compile(r'corrida_\d{8}_\d+')

        batch_name = "trial" # Fallback
        lab_csv = None

        if "sample_sheet.csv" in csv_files:
            df_sample_sheet = pd.read_csv("sample_sheet.csv")
            batch_name = df_sample_sheet["sample_id"][1] 
            lab_csv = batch_name + ".csv"
            os.system(f"mv sample_sheet.csv {lab_csv}")
        else:
            for f in csv_files:
                name = os.path.basename(f)
                if corrida_pattern.search(name):
                    batch_name = name.replace(".csv", "")
                    lab_csv = f
                    break

        if not lab_csv:
            messagebox.showerror("Missing CSV", "No 'corrida_...' sample sheet found in the folder!")
            return

        # 2. Update config and prepare samples.csv
        try:
            self.update_nextflow_config(batch_name)

            df_lab = pd.read_csv(lab_csv)
            nextflow_data = []

            # Check for the 'alias' column we saw in your sample file
            if 'alias' not in df_lab.columns:
                messagebox.showerror("CSV Error", "CSV must contain an 'alias' column.")
                return

            missing_samples = []
            for _, row in df_lab.iterrows():
                alias = str(row['alias'])
                bam_files = glob.glob(os.path.join(input_path, f"{alias}*.bam"))
                bai_files = glob.glob(os.path.join(input_path, f"{alias}*.bai"))

                if bam_files and bai_files:
                    nextflow_data.append({'sample_id': alias, 'bam_path': os.path.abspath(bam_files[0]), 'bai_path': os.path.abspath(bai_files[0])})
                else:
                    missing_samples.append(alias)

            if missing_samples:
                messagebox.showerror("Missing Files", f"Missing .bam or .bai for aliases:\n{', '.join(missing_samples)}")
                return

            if not nextflow_data:
                messagebox.showerror("Matching Error", "No matching BAM files found in the folder.")
                return

            # 3. Save the 'samples.csv' Nextflow expects
            df_nextflow = pd.DataFrame(nextflow_data)
            target_csv = os.path.join(os.getcwd(), "samples.csv")
            df_nextflow.to_csv(target_csv, index=False)

            # 4. Trigger the Bash/Nextflow Execution
            self.status.configure(text=f"🚀 Running QC: {batch_name}...", text_color="orange")
            self.root.update()

            script_path = os.path.join(os.getcwd(), 'run_pipeline.sh')
            # Pass BOTH the batch name and the original folder path to the bash script
            subprocess.Popen(['gnome-terminal', '--', 'bash', script_path, batch_name, input_path])

            # 2. Wait for the Nextflow evaluation file to be created
            eval_csv = os.path.join(os.getcwd(), "results", batch_name, "evaluation", "qc_summary.csv")

            self.status.configure(text="⏳ Processing... Window will pop when ready", text_color="orange")
            self.root.update()

            # Polling loop: Wait up to 10 minutes for the file
            found = False
            for _ in range(900):
                if os.path.exists(eval_csv):
                    found = True
                    break
                time.sleep(1) 
                self.root.update()

            if found:
                # 3. Load data and launch the Selection Window
                df_eval = pd.read_csv(eval_csv)
                selector = SelectionWindow(self.root, df_eval, batch_name)
                self.root.wait_window(selector)
                
                approved_samples = selector.selected_samples
                
                if approved_samples is None:
                    self.status.configure(text="❌ Operation Cancelled", text_color="red")
                else:
                    # ---> CREATE EXPLICIT INCLUDE LIST FOR AWS <---
                    include_file = os.path.join(os.getcwd(), "aws_includes.txt")
                    with open(include_file, "w") as f:
                        # Ensure the sample CSV is also uploaded
                        f.write(f"*{batch_name}*.csv\n")
                        # Add approved BAMs and BAIs
                        for s in approved_samples:
                            f.write(f"*{s}*.bam\n")
                            f.write(f"*{s}*.bai\n")
                            
                    self.status.configure(text=f"✅ {len(approved_samples)} accepted. Ready for AWS Sync.", text_color="green")
            else:
                messagebox.showerror("Timeout", "Nextflow took too long or failed to create the summary.")

        except Exception as e:
            messagebox.showerror("Pipeline Error", f"Failed to process batch:\n{str(e)}")
            self.status.configure(text="❌ Failed", text_color="red")

def main(argv):
    # Initialize the special DnD-enabled window
    root = TkinterDnD.Tk()
    app = LabApp(root)

    if len(argv) == 1:
        app.message("Drop a folder in the box or select via Directory > Open menu.")
    elif len(argv) == 2:
        # Handle argument in ""
        match = re.match(r'^(["\'])(.*)\1$', argv[-1])
        if match:
            the_path = match.group(2)
        else:
            the_path = argv[-1]
        abs_the_path = os.path.abspath(the_path)
        app.selected_dir = abs_the_path
        app.message(f"Dir: {abs_the_path}")
    else:
        print("Error: wrong number of arguments.")
        sys.exit(1)


    app.mainLoop()

if __name__ == "__main__":
    main(sys.argv)
