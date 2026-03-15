"""Alignment Script
Fletcher Falk
Updated March 13, 2026
Script takes input path to folder and runs bbmap onto a passed reference genome
Raw data is located in trimmed-paired directory. Modify to directory of use.
Also reads are required to have specific naming format see line 43."""

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

    """Building sequence list
    Raw data must be organized with same file format or edit below"""
    filepairs = defaultdict(dict)
    print("Pairs of sequences found: \n")
    for file in os.listdir(path):
        if file.endswith("_F_paired_R1.fastq.gz"):
            key = file.replace("_F_paired_R1.fastq.gz", "")
            filepairs[key]["R1"] = file
        elif file.endswith("_R_paired_R2.fastq.gz"):
            key = file.replace("_R_paired_R2.fastq.gz", "")
            filepairs[key]["R2"] = file

    """Making directories for output
    In this case only writing mapped reads prevents overusing disk space"""
    os.mkdir("./trimmed-paired/sam-align")
    os.mkdir("./trimmed-paired/sorted-align")
    os.mkdir("./trimmed-paired/bcfcalls")
    os.mkdir("./trimmed-paired/bcfnorm")
    os.mkdir("./trimmed-paired/bcfinal")
    os.mkdir("./trimmed-paired/consensus-mitogenomes")

    """Run each program on the given sample"""
    for sample, reads in filepairs.items():
        forward = reads["R1"]
        reverse = reads["R2"]
        """Map on vslow algorithm to ensure precision"""
        bbmap = "bbmap.sh in1=" + path + forward + " in2=" + path + reverse + " mappedonly=t ref=" + reference + " vslow k=8 maxindel=200 out=" + path + "sam-align/" + sample + "align.sam"
        sort = "samtools view -bS " + path + "sam-align/" + sample + "align.sam | samtools sort -o " + path + "sorted-align/" + sample + "salign.bam"
        """Convert to bam for BCFtools steps"""
        sindex = "samtools index " + path + "sorted-align/" + sample + "salign.bam"
        execute(bbmap)
        execute(sort)
        execute(sindex)
        """Write consensus for each mitogenome using bcftools
        First call variants and index
        These commands could be piped together but were separated during error testing"""
        bcftools = "bcftools mpileup -Ou --min-MQ 30 --min-BQ 30 -f " + reference + " " + path + "sorted-align/" + sample + "salign.bam | bcftools call -mv -Oz --threads 96 --ploidy 1 -o " + path + "bcfcalls/" + sample + "calls.vcf.gz"
        execute(bcftools)
        bcfindex = "bcftools index --threads " + threads + " " + path + "bcfcalls/" + sample + "calls.vcf.gz"
        execute(bcfindex)
        bcfnorm = "bcftools norm -f " + reference + " --threads " + threads + " " + path + "bcfcalls/" + sample + "calls.vcf.gz -Ob -o " + path + "bcfnorm/" + sample + "calls.norm.vcf.gz"
        execute(bcfnorm)
        bcfilter = "bcftools filter -e 'QUAL<40' --threads " + threads + " " + path + "bcfnorm/" + sample + "calls.norm.vcf.gz -Ob -o " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz"
        execute(bcfilter)
        bcfindex2 = "bcftools index --threads " + threads + " " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz"
        execute(bcfindex2)
        consensus = "cat " + reference + " | bcftools consensus " + path + "bcfinal/" + sample + "calls.norm.flt.vcf.gz > " + path + "consensus-mitogenomes/" + sample + "consensus.fa"
        execute(consensus)
        """Transfer consensus to Geneious: annotate, translate, and align for comparison"""

    """SAMtools stats for review"""
    os.mkdir("samtools-stats")
    for sample, reads in filepairs.items():
        samstats = "samtools stats -@ " + threads + " " + path + "bam-align/" + sample + "align.bam > " + path + "/bwa-samtools-results/" + sample + "_stats.txt"
        execute(samstats)

"""Main: sets up arguments and runs program"""
def main():
    args = argument_parser()
    alignment(args)

"""Start Program"""
if __name__ == "__main__":
    main()
