#!/usr/bin/env python3
import tkinter as tk
import customtkinter as ctk
import pandas as pd

class SelectionWindow(ctk.CTkToplevel):
    def __init__(self, parent, qc_data, batch_name):
        super().__init__(parent)
        self.title("Approval")
        self.geometry("900x700")
        self.attributes("-topmost", True)
        
        self.configure(fg_color="#FFFFFF") 
        self.selected_samples = None
        self.vars = {}

        # 1. TOP: Title
        ctk.CTkLabel(self, text=f"Batch: {batch_name}", 
                     font=("Arial", 22, "bold"), text_color="#000000").pack(pady=20)

        # 2. BOTTOM: Buttons
        btn_f = ctk.CTkFrame(self, fg_color="#F0F0F0", height=100)
        btn_f.pack(side="bottom", fill="x")

        ctk.CTkButton(btn_f, text="CANCEL & EXIT", fg_color="#666666", width=200, height=45,
                      text_color="white", command=self.destroy).pack(side="left", padx=100, pady=25)
        
        ctk.CTkButton(btn_f, text="CONFIRM UPLOAD", fg_color="#2ecc71", width=200, height=45,
                      text_color="white", command=self.on_confirm).pack(side="right", padx=100, pady=25)

        # 3. MIDDLE: Scroll Frame
        self.scroll = ctk.CTkScrollableFrame(self, fg_color="#FFFFFF", border_width=1)
        self.scroll.pack(padx=20, pady=10, fill="both", expand=True)

        # --- ADDING THE TABLE HEADER ---
        header_f = ctk.CTkFrame(self.scroll, fg_color="transparent")
        header_f.pack(fill="x", padx=10, pady=5)
        
        # Grid columns for alignment
        ctk.CTkLabel(header_f, text="UP?", width=50, text_color="#000000", font=("Arial", 12, "bold")).grid(row=0, column=0)
        ctk.CTkLabel(header_f, text="SAMPLE ID", width=150, text_color="#000000", font=("Arial", 12, "bold"), anchor="w").grid(row=0, column=1)
        ctk.CTkLabel(header_f, text="STATUS", width=100, text_color="#000000", font=("Arial", 12, "bold"), anchor="w").grid(row=0, column=2)
        ctk.CTkLabel(header_f, text="REASONS", width=300, text_color="#000000", font=("Arial", 12, "bold"), anchor="w").grid(row=0, column=3)

        # --- THE MISSING LOOP: Populating the samples ---
        for _, row in qc_data.iterrows():
            s_id = row['sample_id']
            status = row['status']
            reason = row['fail_reasons'] if pd.notna(row['fail_reasons']) else ""
            
            row_f = ctk.CTkFrame(self.scroll, fg_color="transparent")
            row_f.pack(fill="x", padx=10, pady=2)

            var = tk.BooleanVar(value=(status == "PASS"))
            self.vars[s_id] = var
            
            # Checkbox
            ctk.CTkCheckBox(row_f, text="", variable=var, width=50).grid(row=0, column=0)
            
            # Sample ID (Forced Black)
            ctk.CTkLabel(row_f, text=s_id, width=150, anchor="w", 
                         text_color="#000000", font=("Arial", 12)).grid(row=0, column=1)

            # Status (Green for PASS, Red for FAIL)
            st_color = "#008000" if status == "PASS" else "#FF0000"
            ctk.CTkLabel(row_f, text=status, width=100, anchor="w", 
                         text_color=st_color, font=("Arial", 12, "bold")).grid(row=0, column=2)

            # Reasons (Forced Dark Gray)
            ctk.CTkLabel(row_f, text=reason, text_color="#444444", 
                         font=("Arial", 11), anchor="w", wraplength=350).grid(row=0, column=3)

    def on_confirm(self):
        self.selected_samples = [s for s, v in self.vars.items() if v.get()]
        self.destroy()
