#for parsing arguments
import argparse
import sys

#for executing Unix commands in python script
import os

from Bio import Entrez #for downloading data from NCBI database
from Bio import SeqIO #for parsing through genbank file

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="build a transctiptome index for HCMV") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    parser.add_argument('-a', '--accession', help='output file name', required='True') #accession number argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
dir = args.dir #saves each arguement
accession = args.accession

#use entrez to find reference sequence fasta file
Entrez.email ="ckrueger2@luc.edu" #needed to use esearch 
handle = Entrez.efetch(db = "nucleotide", id = accession, rettype = 'gb', retmode = "text") #search nucleotide database

#isolate coding regions
cds = [] #initialize empty list
for region in SeqIO.parse(handle, "genbank"): #for each section in the genbank file
    for feature in region.features: #for each feature in the section
        if feature.type == "CDS": #if the feature is a CDS
            id = feature.qualifiers['protein_id'][0] #save the RefSeq protein_id
            sequence = feature.extract(region.seq) #save the amino acid sequence
            cds.append(f'>{id}\n{sequence}') #add id and sequence to the cds list

#this is how many cds there are
num_cds = (len(cds)) 

#where to save the cds file to
output_file = f"{dir}/cds.fasta" 

#save fasta file
with open(output_file, "w") as f:
    f.write("\n".join(cds) + "\n")

#gzip fasta file
os.system(f"gzip {dir}/cds.fasta")

#save how many cds are present to log file
with open((f'{dir}/PipelineProject.log'), 'w') as f:
    f.write(f"The HCMV genome ({accession}) has {num_cds} CDS.\n\n")