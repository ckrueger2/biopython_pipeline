#for parsing arguments
import argparse
import sys

#for executing Unix commands in python script
import os

#for downloading files from NCBI databse
from Bio import SeqIO

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="quantify the TPM of each CDS") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    parser.add_argument('-o', '--organism', help='path to input and output file directory', required='True') #organism group argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
infile = args.input #saves each arguement
dir = args.dir
organism = args.organism
organism = organism.strip()

#save each SRR ID to a list
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create donor/patient list
patients = []
for read in reads:
    patient, condition, sample = read.strip().split(",")
    patient_clean = patient.strip().replace(' ', '')
    if patient_clean not in patients: #only unique values added
        patients.append(patient_clean)

#make directory for blast output
os.system(f"mkdir {dir}/blast")

#download dataset for organism and make a database from the files
os.system(f"datasets download virus genome taxon {organism.lower()} --include genome --filename {dir}/blast/ncbi_dataset.zip")
os.system(f"unzip -o {dir}/blast/ncbi_dataset.zip -d {dir}/blast") #unzip downloaded file
os.system(f"makeblastdb -in {dir}/blast/ncbi_dataset/data/genomic.fna -out {dir}/blast/{organism} -title {organism} -dbtype nucl") #make local database of genomic data called the organisms name

#find longest contig
for patient in patients: ##parse through each donor
    with open(f"{dir}/{patient}_assembly/contigs.fasta") as f: 
        contigs =  SeqIO.parse(f, "fasta") #save each contig
        longest_contig = None #start with no longest contig
        max_length = 0 #and a longest length of zero
        for contig in contigs: #parse through each contig
            contig_length = len(contig.seq)
            if contig_length > max_length: #save if new longest contig
                longest_contig = contig
                max_length = contig_length
        if longest_contig: #once longest contig is found
            with open(f"{dir}/{patient}_assembly/{patient}_longest_contig.fasta", "w") as g: #save to a fasta file
                SeqIO.write(longest_contig, g, "fasta")
    query_seqfile = f"{dir}/{patient}_assembly/{patient}_longest_contig.fasta" #save path to longest contig
    #blasting from Python, using local dictonary; --max_hsps 1 saves only best alignment from each query-subject pair of sequences
    blast_command = f'blastn -query {query_seqfile} -db {dir}/blast/{organism} -out {dir}/blast/{patient}_{organism}.tsv -max_hsps 1 -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    os.system(blast_command)

    #define the headers
    headers = "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n"

    #output file path
    output_file = f"{dir}/{patient}_blast_assembly"

    #read the original content of the tsv file
    with open(f"{dir}/blast/{patient}_{organism}.tsv", 'r') as g:
        file = g.readlines()

    #write to log
    with open((f'{dir}/PipelineProject.log'), 'a') as f: 
        f.write(f"\n{patient}:\n{headers}{''.join(file[:10])}")