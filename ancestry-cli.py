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

Mandatory arguments are:
VCF - Path to VCF file.
OUTPUT.txt - Path to output text file.
COMPILED-TSV-PATH - Use this argument to declare where the compiled, reference tsv is or should be stored after building it. If tsv already exists,
                        the manifest and reference will be ignored.
"""

import os
import sys
import re


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
        "VCF": "",
        "TXT": "",
        "TSV": "",
        "--manifest": "",
        "--reference": ""
    }

    if len(argv) < 4 or len(argv) > 8:
        print(f"Inocrrect number of arguments {len(argv)}.")
        usage(argv)
        sys.exit(1)
    elif not re.search(r"\.vcf(?:\.gz)?$", argv[1]):
        print(f"Invalid input VCF: {argv[1]} (supported formats: '.vcf', '.vcf.gz'")
        usage(argv)
        sys.exit(1)
    elif not os.path.isfile(argv[1]):
        print(f"Error: The file '{argv[1]}' does not exist!")
        sys.exit(1)
    elif True:#TODO:Aqui vas
    else:
        i = 3
        while i < len(argv):
            val = argv[i+1]
            if argv[i] = "--manifest":

        # Validate file extension
        for filename in input_files:
            if not re.search(r"\.dat(?:\.gz)?$", filename):
                print(f"Invalid filename: {filename} (supported formats: '.dat', '.dat.gz'")
                usage(argv)
                sys.exit(1)
            # Validate file existence
            if not os.path.isfile(filename):
                print(f"Error: The file '{filename}' does not exist!")
                sys.exit(1)
            valid_arguments["FILE"].append(filename)

    #print(valid_arguments)
    return valid_arguments


if __name__ == "__main__":
    argv = sys.argv
    arguments = main(argv)

    if argv[1] == '--doi':
        print(f"GO ID: {arguments['--doi']}, File: {arguments['FILE']}")
    elif argv[1] == '--pubmed':
        print(f"GO ID: {arguments['--pubmed']}, File: {arguments['FILE']}")
