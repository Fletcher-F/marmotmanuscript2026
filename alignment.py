"""Alignment Script
Fletcher Falk
Updated March 6, 2026
Script takes input path to folder and runs bwa-mem2 onto a passed reference genome"""

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
    parser = argparse.ArgumentParser(description="Alignment Argument List")
    parser.add_argument("--input", "-i",
                        help="Path to folder", required=True)
    parser.add_argument("--reference", "-r",
                        help="Specify reference genome", required=True)
    parser.add_argument("--threads", "-t",
                        help="Specify total number of threads to use. Default is 1.", default=1)
    """Return arguments"""
    return parser.parse_args()

"""Main"""
def alignment(args):
    """Assign arguments"""
    path = args.input
    reference = args.reference
    threads = args.threads

    """Start"""
    print("--- Running Alignment --- \n",
    "Using ", threads, " threads... \n", "File(s) path = ", path, "\n",
    "Reference genome: ", reference, "\n")

    """Running BWA"""
    filepairs = defaultdict(dict)
    print("Pairs of sequences found: \n")
    for file in os.listdir(path):
        if file.endswith("_F_paired_R1.fastq.gz"):
            key = file.replace("_F_paired_R1.fastq.gz", "")
            filepairs[key]["R1"] = file
        elif file.endswith("_R_paired_R2.fastq.gz"):
            key = file.replace("_R_paired_R2.fastq.gz", "")
            filepairs[key]["R2"] = file

    print("Indexing reference genome for BWA")
    """Using bwa-mem2"""
    bwaindex = "bwa-mem2 index " + reference
    execute(bwaindex)

    """Making directories for output bam"""
    os.mkdir("./trimmed-paired/bam-align")
    os.mkdir("./trimmed-paired/bcfcalls")
    os.mkdir("./trimmed-paired/bcfnorm")
    os.mkdir("./trimmed-paired/bcfinal")
    os.mkdir("./trimmed-paired/consensus-mitogenomes")

    """Run each program on the given sample"""
    for sample, reads in filepairs.items():
        forward = reads["R1"]
        reverse = reads["R2"]
        bwa = "bwa-mem2 mem -t " + threads + " " + reference + " " + path + forward + " " + path + reverse + " | samtools sort --write-index -@ " + threads + " -o " + path + "bam-align/" + sample + "align.bam"
        print("Running bwa-mem2 on: " + sample)
        execute(bwa)
        """Write consensus for each mitogenome using bcftools
        Call variants and index
        These commands could be piped together but are separate for error testing"""
        bcftools = "bcftools mpileup -Ou --min-MQ 30 --min-BQ 30 -f " + reference + " " + path + "bam-align/" + sample + "align.bam | bcftools call -mv -Oz --threads 96 --ploidy 1 -o " + path + "bcfcalls/" + sample + "calls.vcf.gz"
        execute(bcftools)
        bcfindex = "bcftools index --threads " + threads + " " + path + "bcfcalls/" + sample + "calls.vcf.gz"
        execute(bcfindex)
        bcfnorm = "bcftools norm -f " + reference + " --threads " + threads + " " + path + "bcfcalls/" + sample + "calls.vcf.gz -Ob -o " + path + "bcfnorm/" + sample + "calls.norm.vcf.gz"
        execute(bcfnorm)
        bcfilter = "bcftools filter -e 'QUAL<30' --threads " + threads + " " + path + "bcfnorm/" + sample + "calls.norm.vcf.gz -Ob -o " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz"
        execute(bcfilter)
        bcfindex2 = "bcftools index --threads " + threads + " " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz"
        execute(bcfindex2)
        consensus = "cat " + reference + " | bcftools consensus " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz > " + path + "consensus-mitogenomes/" + sample + "consensus.fa"
        execute(consensus)

    """SAMtools stats for review"""
    os.mkdir("bwa-samtools-results")
    for sample, reads in filepairs.items():
        bwasamtools = "samtools stats -@ " + threads + " " + path + "bam-align/" + sample + "align.bam > " + path + "/bwa-samtools-results/" + sample + "_stats.txt"
        execute(bwasamtools)

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    alignment(args)

"""Start Program"""
if __name__ == "__main__":
    main()