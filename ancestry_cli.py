#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# This script is solely for the purpose of establishing the console app and checking arguments
"""
$ python3 ancestry-cli.py
Usage: ancestry-cli.py VCF | OUTPUT.txt | COMPILED.tsv [--help |--manifest MANIFEST.csv | --reference REF.csv]
GSA-1 parser for imputation data by Roberto Zarazuanav, 2026.
Output an Ancestry-like .tsv using data derived from QUILT2 for lcWGS and using Illumina GSA-1 manifest as reference.
-----------------------------------------------

Optional arguments are:
--help - Display this help page.
--manifest GSA-1-MANIFEST-PATH - Path to GSA-1 manifest csv. Must be accompanied of a reference.
--reference GSA-1-RSIDS-REFERENCE-PATH - Path to GSA-1 reference txt. Must be accompanied of a manifest.
--compile PATH-COMPILED-REFERENCE.TSV - Path to a previously created compiled GSA-1 reference using the AncestryWriter class. If file in input path does not
                                         exist, one will be created.
Mandatory arguments are:
VCF - Path to VCF file.
OUTPUT.txt - Path to output text file. Default name: grallow-ancestry-output.txt
COMPILED-TSV-PATH - Use this argument to declare where the compiled, reference tsv is or should be stored after building it. If tsv already exists,
                        the manifest and reference will be ignored.
"""

import os
import sys
import re
from pathlib import Path


def usage(argv):
    """Show the usage line."""
    print(f"Usage: {argv[0]} VCF | OUTPUT.txt | COMPILED.tsv [--help | --manifest MANIFEST.csv | --reference REF.csv]")


def mhelp(argv):
    """Display help menu."""
    print(__doc__.format(argv[0]))


def main(argv):
    """Main function to handle argument validation."""
    if len(argv) == 1:
        usage(argv)
        sys.exit(1)

    elif "--help" in argv:
        mhelp(argv)
        sys.exit(1)

    
    valid = ["--manifest","--reference"]
    # Initialize dictionary to store arguments
    valid_arguments = {
        "VCF": None,
        "TXT": None,
        "TSV": None,
        "--compile": None,
        "--manifest": None,
        "--reference": None
    }

    if len(argv) < 3:
        print("Error: Minimum arguments required are input VCF and output TXT.")
        usage(argv)
        sys.exit(1)
    
    input_vcf = Path(argv[1])
    if not re.search(r"\.vcf(?:\.gz)?$", input_vcf.name):
        print(f"Invalid input VCF: {input_vcf.name} (supported formats: '.vcf', '.vcf.gz')")
        usage(argv)
        sys.exit(1)
    if not input_vcf.is_file():
        print(f"Error: The file '{input_vcf}' does not exist!")
        sys.exit(1)
    valid_arguments["VCF"] = input_vcf

    output_arg = Path(argv[2])
    if re.search(r"\.txt(?:\.gz)?$", output_arg.name):
        output_filename = output_arg
        target_dir = output_filename.parent
    else:
        target_dir = output_arg
        output_filename = target_dir / "grallow-ancestry-output.txt"
    
    valid_arguments["TXT"] = output_filename
    
    if str(target_dir) != ".":
        target_dir.mkdir(parents=True, exist_ok=True)
    
    i = 3
    while i < len(argv):
        if i + 1 >= len(argv):
            print(f"Error: Missing file path after {argv[i]}")
            usage(argv)
            sys.exit(1)
        
        flag = argv[i]
        val = Path(argv[i+1])

        if not flag in ["--manifest", "--reference", "--compile"]:
            print(f"Invalid optional argument {argv[i]}, '--manifest' and '--reference' are the valid options.")
            usage(argv)
            sys.exit(1)
        elif flag == "--manifest":
            if not re.search(r"\.csv(?:\.gz)?$", val.name):
                print(f"Invalid input GSA-1 manifest: {val.name} (supported formats: '.csv', '.csv.gz')")
                usage(argv)
                sys.exit(1)
            elif not val.is_file():
                print(f"{val.name} does not exist.")
                usage(argv)
                sys.exit(1)
            else:
                valid_arguments["--manifest"] = val

        elif flag == "--reference":
            if not re.search(r"\.txt(?:\.gz)?$", val.name):
                print(f"Invalid input GSA-1 rsid reference: {val.name} (supported formats: '.txt', '.txt.gz')")
                usage(argv)
                sys.exit(1)
            elif not val.is_file():
                print(f"{val.name} does not exist.")
                usage(argv)
                sys.exit(1)
            else:
                valid_arguments["--reference"] = val

        elif flag == "--compile":
            if not re.search(r"\.tsv(?:\.gz)?$", val.name):
                print(f"Error: Invalid compiled TSV path: {val.name} (supported formats: '.tsv', '.tsv.gz')")
                sys.exit(1)
            valid_arguments["--compile"] = val
            
        i+=2
    
    tsv_path = valid_arguments["--compile"]
    if tsv_path and tsv_path.is_file():
        valid_arguments["TSV"] = tsv_path
    else:
        if not valid_arguments["--manifest"] or not valid_arguments["--reference"]:
            print("Error: You must provide an existing '--compile' TSV, OR both '--manifest' and '--reference' files to generate one.")
            usage(argv)
            sys.exit(1)

        if not tsv_path:
            valid_arguments["--compile"] = target_dir / "compiled_gsa_map.tsv"
            print(f"Notice: No '--compile' path provided. Defaulting creation to '{valid_arguments['--compile']}'")

        valid_arguments["TSV"] = None

    tsv_target_dir = valid_arguments["--compile"].parent
    if str(tsv_target_dir) != ".":
        tsv_target_dir.mkdir(parents=True, exist_ok=True)

    #print(valid_arguments)
    return valid_arguments


if __name__ == "__main__":
    argv = sys.argv
    arguments = main(argv)
