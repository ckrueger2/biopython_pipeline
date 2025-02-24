#for parsing arguments
import argparse 
import sys

#for executing Unix commands in python script
import os

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="retrieve fastq files from NCBI") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input file directory arguement, 'True' means the argument must be provided
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
infile =args.input #saves each arguement
dir = args.dir

#open input file and save by line
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create lists of sampleids and timepoints/conditions
samples = [] #intitialize empty list
conditions = []
for read in reads: #for each line in input file
    patient, condition, sample = read.strip().split(",") #input file is comma-delimited
    conditions.append(condition.strip()) #strip any whitespace
    samples.append(sample.strip())

#parse through each sample id
for sample in samples:
    #download SRA files from the NCBI database
    os.system(f"wget -P {dir}/reads https://sra-pub-run-odp.s3.amazonaws.com/sra/{sample}/{sample}")
    #split downloaded SRA file to retrieve paired-end read fasta files
    os.system(f"fasterq-dump --split-3 --threads 2 -O {dir}/reads/ {dir}/reads/{sample}")
    #gzip fastq files
    os.system(f"gzip {dir}/reads/{sample}_1.fastq")
    os.system(f"gzip {dir}/reads/{sample}_2.fastq")