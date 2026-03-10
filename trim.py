"""Trim Raw Reads Script
Fletcher Falk
Updated March 6, 2026
Script takes input path to folder (of fastq.gz files) and runs FastQC, trims, 
and FastQC with standard settings."""

"""Note: Adapted for specific adapters. Run one FastQC test to see what is needed
and modify the name of the adapters file on line 89"""

"""Importing libraries"""
import os, subprocess, argparse, sys
from pathlib import Path
from collections import defaultdict

"""Function to execute bash commands"""
def execute(command):
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    print(result.stdout, "\n", result.stderr)

"""Arguments for running the program"""
def argument_parser():
    """Arguments"""
    parser = argparse.ArgumentParser(description="Trim Argument List")
    parser.add_argument("--input", "-i",
                        help="Path to folder", required=True)
    parser.add_argument("--threads", "-t",
                        help="Specify total number of threads to use. Default is 1.", default=1)
    """Return arguments"""
    return parser.parse_args()

"""Check valid fasta
Validates specific input format specified in GitHub."""
def valid_fastq(file):
    if os.path.splitext(file)[1] not in (".fastq.gz"):
        print("Invalid fastq.gz files in directory...")
        sys.exit()

"""Check valid directory"""
def valid_dir(path):
    """Check if directory exists before starting"""
    if Path(path).is_dir() == False:
        print("Error... input is not a directory...")
        sys.exit()
    files = os.listdir(path)
    """Check valid fasta files in directory"""
    for file in files:
        valid_fastq(file)

"""Main"""
def trim(args):
    """Assign arguments"""
    path = args.input
    threads = args.threads
    """Check inputs"""
    valid_dir(path)

    """Start"""
    print("--- Running Trimmer --- \n",
    "Using ", threads, " threads... \n", "File(s) path = ", path, "\n")

    """Running FastQC on input path"""
    print("Initial Quality Check")
    os.mkdir("initial-fastqc")
    qualitycheck = "fastqc -t " + threads + " " + path + "*.fastq.gz -o initial-fastqc"
    execute(qualitycheck)

    """Running Trimmomatic"""
    os.mkdir("trimmed-paired")
    os.mkdir("trimmed-unpaired")
    print("Running trimming: ensure to set adapter sequences. \n")
    print("Currently using: ")

    filepairs = defaultdict(dict)
    print("Pairs of sequences found: \n")
    for file in os.listdir(path):
        if file.endswith("-MarmotData-R1.fastq.gz"):
            key = file.replace("-MarmotData-R1.fastq.gz", "")
            filepairs[key]["R1"] = file
        elif file.endswith("-MarmotData-R2.fastq.gz"):
            key = file.replace("-MarmotData-R2.fastq.gz", "")
            filepairs[key]["R2"] = file

    for sample, reads in filepairs.items():
        forward = reads["R1"]
        reverse = reads["R2"]
        trim = "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads " + threads + " " + path + forward + " " + path + reverse + " " \
            "./trimmed-paired/" + sample + "_F_paired_R1.fastq.gz " + "./trimmed-unpaired/" + sample + "_F_unpaired_R1.fastq.gz " + \
            "./trimmed-paired/" + sample + "_R_paired_R2.fastq.gz " + "./trimmed-unpaired/" + sample + "_R_unpaired_R2.fastq.gz " + \
            "ILLUMINACLIP:./adapters/Illumina-Adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
        execute(trim)

    """Running FastQC on trimmed reads"""
    os.mkdir("final-fastqc")
    qualitycheck = "fastqc -t " + threads + " ./trimmed-paired/*.fastq.gz -o final-fastqc"
    execute(qualitycheck)

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    trim(args)

"""Start Program"""
if __name__ == "__main__":
    main()
