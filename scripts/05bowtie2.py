#for parsing arguments
import argparse
import sys

#for executing Unix commands in python script
import os

#to work with a gzipped file
import gzip

#do download FASTA file from NCBI database
from Bio import Entrez

def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="compare strains") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #input argument
    parser.add_argument('-a', '--accession', help='output file name', required='True') #accession number arguement
    parser.add_argument('-p', '--path', help='path to existing fastq file') #path to external fastq files argument
    parser.add_argument("--no-wget", action="store_true", help="Use alternate file path") #if external fastq files are used argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
dir = args.dir #saves each arguement
infile = args.input
path = args.path
accession = args.accession

#save each SRR ID to a list
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create lists of sampleids, timepoints/conditions, and donor ids
samples = []
conditions = []
patients = []
for read in reads:
    patient, condition, sample = read.strip().split(",")
    samples.append(sample.strip())
    conditions.append(condition.strip())
    patients.append(patient.strip())

#create new directory for bowtie outputs
os.system(f"mkdir {dir}/bowtie")

#retrive HCMV fasta file
Entrez.email ="ckrueger2@luc.edu" #needed to use esearch 
handle = Entrez.efetch(db = "nucleotide", id = accession, rettype = 'fasta', retmode = "text") #search nucleotide database
with open(f"{dir}/refseq.fasta", "w") as f: #save downloaded fasta file
    f.write(handle.read())

#create genome index
os.system(f"bowtie2-build {dir}/refseq.fasta {dir}/bowtie/bowtie_build")

#write sam files and find only reads that map to reference
for sample in samples: #for each sample
    if args.no_wget: #if external fastq files are used
        #flags; supress non-essential output, path to bowtie2 index, path to first fastq, path to second fastq, output file
        os.system(f"bowtie2 --quiet -x {dir}/bowtie/bowtie_build -1 {path}/{sample}_1.fastq.gz -2 {path}/{sample}_2.fastq.gz -S {dir}/bowtie/{sample}.sam")
        os.system(f"bowtie2 -x {dir}/bowtie/bowtie_build -1 {path}/{sample}_1.fastq.gz -2 {path}/{sample}_2.fastq.gz -S {dir}/bowtie/{sample}.sam --al-conc-gz {dir}/bowtie/{sample}_mapped%.fastq.gz")
    else: #if files downloaded from NCBI in pipeline are used
        os.system(f"bowtie2 --quiet -x {dir}/bowtie/bowtie_build -1 {dir}/reads/{sample}_1.fastq.gz -2 {dir}/reads/{sample}_2.fastq.gz -S {dir}/bowtie/{sample}.sam")
        os.system(f"bowtie2 -x {dir}/bowtie/bowtie_build -1 {dir}/reads/{sample}_1.fastq.gz -2 {dir}/reads/{sample}_2.fastq.gz -S {dir}/bowtie/{sample}.sam --al-conc-gz {dir}/bowtie/{sample}_mapped%.fastq.gz")

#add spacer in log file
with open((f'{dir}/PipelineProject.log'), 'a') as f: 
    f.write("\n")

#find number of reads before and after Bowtie2 filtering
for patient, condition, sample in zip(patients, conditions, samples):  #combine the three lists
    if args.no_wget: #if external fastq files are used
        with gzip.open((f"{path}/{sample}_1.fastq.gz"), "r") as f: #pre-bowtie file; first is equivalent to second
            pre_reads = f.readlines()  #read each line
            pre_num_reads = len(pre_reads) // 4  #fastq format is 4 lines per read
        with gzip.open((f"{dir}/bowtie/{sample}_mapped1.fastq.gz"), "r") as g: #post-bowtie file
            post_reads = g.readlines()  
            post_num_reads = len(post_reads) // 4 
        #write pre and post read pairs to log file
        with open((f'{dir}/PipelineProject.log'), 'a') as h: #"with" closes file after it is written to, "open()" opens the file path, "w" opens the file and writes to it, overwritting if the file already exists, "as f" assigns file to f object for f.write
            h.write(f"{patient} ({condition}) had {pre_num_reads} read pairs before Bowtie2 filtering and {post_num_reads} read pairs after.\n")
    else: #if files downloaded from NCBI in pipeline are used
        with gzip.open((f"{dir}/reads/{sample}_1.fastq.gz"), "r") as f:
            pre_reads = f.readlines() 
            pre_num_reads = len(pre_reads) // 4 
        with gzip.open((f"{dir}/bowtie/{sample}_mapped1.fastq.gz"), "r") as g:
            post_reads = g.readlines()
            post_num_reads = len(post_reads) // 4 
        with open((f'{dir}/PipelineProject.log'), 'a') as h:
            h.write(f"{patient} ({condition}) had {pre_num_reads} read pairs before Bowtie2 filtering and {post_num_reads} read pairs after.\n")