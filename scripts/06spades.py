#for parsing arguments
import argparse
import sys

#for executing Unix commands in python script
import os

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="run spades") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    parser.add_argument('-k', '--kmer', help='path to input and output file directory', required='True') #kmer length argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
infile = args.input #saves each arguement
dir = args.dir
kmer = args.kmer

#save each SRR ID to a list
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create dictionary of patient:list of samples and another list of samples
samples = []
patient_samp_dict = {}
for read in reads:
    patient, condition, sample = read.strip().split(",")
    patient = patient.strip()
    sample = sample.strip()
    samples.append(sample) #add sample to list
    #if patient is already in dictionary, append the sample
    if patient in patient_samp_dict:
        patient_samp_dict[patient].append(sample) #add to existing list
    else:
        patient_samp_dict[patient] = [sample] #start new list

#add spacer in log file
with open((f'{dir}/PipelineProject.log'), 'a') as f:
    f.write("\n")

#run spades for each donor/patient
for patient, patient_samples in patient_samp_dict.items(): #parse through dictionary
    command = [] #initialize empty list
    for i, sample in enumerate(patient_samples): #parse through each list within a key of the dictionary
        command.append(f"--pe{i + 1}-1 {dir}/bowtie/{sample}_mapped1.fastq.gz") #command for first paired-end fastq
        command.append(f"--pe{i + 1}-2 {dir}/bowtie/{sample}_mapped2.fastq.gz") #command for second paired-end fastq
        patient_assembly_dir = f"{dir}/{patient.replace(' ', '')}_assembly" #directory to output spades to
    #perform spades assembly with every sample each donor/patient has
    os.system(f"spades.py -k {kmer} -t 2 --only-assembler {' '.join(command)} -o {patient_assembly_dir}/") #join each paired-end fastq command
    #save unix command to the log file
    with open((f'{dir}/PipelineProject.log'), 'a') as f:
        f.write(f"spades.py -k {kmer} -t 2 --only-assembler {' '.join(command)} -o {patient_assembly_dir}/\n")