# Grallow Lab QC Pipeline

The Grallow Lab QC Pipeline is a semi-automated bioinformatics workflow designed to process raw sequencing data, perform quality control (QC) checks, and facilitate a human-in-the-loop review process before uploading data to AWS S3.

## Overview

This tool chain coordinates three main stages:
1. **Automated QC (Nextflow):** A Nextflow pipeline runs `mosdepth`, `flagstat`, and `NanoPlot` to assess sequencing quality[cite: 6].
2. **QC Evaluation (Python):** `evaluate_qc.py` parses MultiQC and Mosdepth summaries to flag samples based on defined thresholds including depth, mapping rate, N50, read quality, and sex-discordance checks.
3. **Human Review (GUI):** A Tkinter-based GUI allows users to visually inspect results and selectively approve samples for final cloud synchronization[cite: 2, 3].

---

## Project Structure

* `import_files.py`: The main GUI launcher and entry point for starting the QC process.
* `main.nf`: Nextflow pipeline orchestration script[cite: 6].
* `nextflow.config`: Configuration profile managing local Conda environments and cloud Docker containers.
* `evaluate_qc.py`: Logic parser for evaluating sample metrics against quality control thresholds[cite: 9].
* `run_pipeline.sh`: Bash wrapper script that executes the Nextflow pipeline, prompts for AWS SSO login, and handles the S3 upload configuration[cite: 8].
* `selector.py`: UI component for the sample approval and verification table.
* `statusbar.py`: Status and progress tracking widget helper[cite: 4].
* `guibaseclass.py`: Base architectural template class for the application interface[cite: 5].

---

## Requirements

* **Nextflow**: Required to run the primary QC processes[cite: 6].
* **Conda / Docker**: Used to manage tool environments seamlessly across local or cloud setups[cite: 7].
* **AWS CLI**: Required for syncing finalized sample batches to your target cloud infrastructure[cite: 8].

---

## Usage

### 1. Launching the GUI
Start the pipeline interface by executing the master script from your terminal:
```bash
python3 import_files.py
