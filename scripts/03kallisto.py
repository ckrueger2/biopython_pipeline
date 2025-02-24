#for parsing arguments
import argparse
import sys

#for executing Unix commands in python script
import os

#for calculating median and mean
import statistics

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="run kallisto") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    parser.add_argument('-p', '--path', help='path to existing fastq file') #path to external fastq files argument
    parser.add_argument("--no-wget", action="store_true", help="Use alternate file path") #if external fastq files are used argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
dir = args.dir #saves each arguement
infile = args.input
path = args.path

#create an index of the transcriptome reference needed for kallisto quantification
os.system(f'kallisto index -i {dir}/index.idx {dir}/cds.fasta.gz')

#save each SRR ID to a list
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create lists of sampleids and timepoints/conditions
samples = [] #intitialize empty list
conditions = []
for read in reads: #for each line in input file
    patient, condition, sample = read.strip().split(",") #input file is comma-delimited
    conditions.append(condition.strip()) #strip any whitespace
    samples.append(sample.strip())

#parse through SRRs to generate abundance.tsvs for each SRR using 30 bootstraps and 2 cores
for sample in samples:
    if args.no_wget: #if external fastq files are used
        os.system(f"time kallisto quant -i {dir}/index.idx -o {dir}/quantified_{sample} -b 30 -t 2 {path}/{sample}_1.fastq.gz {path}/{sample}_2.fastq.gz")
    else: #if files downloaded from NCBI in pipeline are used
        os.system(f"time kallisto quant -i {dir}/index.idx -o {dir}/quantified_{sample} -b 30 -t 2 {dir}/reads/{sample}_1.fastq.gz {dir}/reads/{sample}_2.fastq.gz")

#save header line to log file
with open((f'{dir}/PipelineProject.log'), 'a') as f:
    f.write(f"sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

#parse through sample-condition pairs to retrieve statistics
for sample, condition in zip(samples, conditions):
    with open(f"{dir}/quantified_{sample}/abundance.tsv", "r") as g: #using kallisto abundance file
        header = g.readline() #save each line
        tpm_nums = [] #initialize empty list
        for line in g: #parse through each line
            column = line.strip().split("\t") #to index
            tpm = float(column[4]) #save with decimal place
            tpm_nums.append(tpm) #save all tpm values to list
        min_tpm = min(tpm_nums) #minimum tpm value
        med_tpm = statistics.median(tpm_nums) #median tpm value
        mean_tpm = statistics.mean(tpm_nums)  #mean tpm value
        max_tpm = max(tpm_nums) #max tpm value
    #save each statistic to log file, tab delimited
    with open((f'{dir}/PipelineProject.log'), 'a') as f:
        f.write(f"{sample}\t{condition}\t{min_tpm}\t{med_tpm}\t{mean_tpm}\t{max_tpm}\n")