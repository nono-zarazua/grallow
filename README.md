# Grallow Lab QC Pipeline

The Grallow Lab QC Pipeline is a semi-automated bioinformatics workflow designed to process raw sequencing data, perform quality control (QC) checks, and facilitate a human-in-the-loop review process before uploading data to AWS S3.

## Overview

This tool chain coordinates three main stages:
1. **Automated QC (Nextflow):** A Nextflow pipeline runs `mosdepth`, `flagstat`, and `NanoPlot` to assess sequencing quality.
2. **QC Evaluation (Python):** `evaluate_qc.py` parses MultiQC and Mosdepth summaries to flag samples based on defined thresholds including depth, mapping rate, N50, read quality, and sex-discordance checks.
3. **Human Review (GUI):** A Tkinter-based GUI allows users to visually inspect results and selectively approve samples for final cloud synchronization.

---

## Project Structure

* `import_files.py`: The main GUI launcher and entry point for starting the QC process.
* `main.nf`: Nextflow pipeline orchestration script.
* `nextflow.config`: Configuration profile managing local Conda environments and cloud Docker containers.
* `evaluate_qc.py`: Logic parser for evaluating sample metrics against quality control thresholds.
* `run_pipeline.sh`: Bash wrapper script that executes the Nextflow pipeline, prompts for AWS SSO login, and handles the S3 upload configuration.
* `selector.py`: UI component for the sample approval and verification table.
* `statusbar.py`: Status and progress tracking widget helper.
* `guibaseclass.py`: Base architectural template class for the application interface.

---

## Requirements

* **Nextflow**: Required to run the primary QC processes.
* **Conda / Docker**: Used to manage tool environments seamlessly across local or cloud setups.
* **AWS CLI**: Required for syncing finalized sample batches to your target cloud infrastructure.

---

## Usage

### 1. Launching the GUI
Start the pipeline interface by executing the master script from your terminal:
```bash
python3 import_files.py
