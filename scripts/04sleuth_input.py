#for parsing arguments
import argparse
import sys

#funtion to parse command line arguments
def check_arg(args=None): #defines function named 'check_arg' with 'args' argument that is used when no file is provided, meaning providing an argument is optional and script will still run
    parser = argparse.ArgumentParser(description="write sleuth input file") #creates argument parser for command line arguments, description shows purpose of the parsing argument
    parser.add_argument('-i', '--input', help='path to input file', required='True') #input argument
    parser.add_argument('-d', '--dir', help='path to input and output file directory', required='True') #directory argument
    return parser.parse_args(args) #parses the provided arguments and returns it to be used in the script

#retrieve command line arguments and assign them to variables
args = check_arg(sys.argv[1:]) #takes in all provided arguments
dir = args.dir #saves each arguement
infile = args.input

#save each SRR ID to a list
with open(infile, "r") as f:
    reads = f.read().strip().split("\n") #read in file, strip any whitespace, split by line

#create lists of sampleids, timepoints/conditions, and make path to kallisto output list
samples = []
conditions = []
paths = []
for read in reads:
    patient, condition, sample = read.strip().split(",")
    samples.append(sample.strip())
    conditions.append(condition.strip())
    path = (f"{dir}/quantified_{sample.strip()}") #kallisto output
    paths.append(path)

#make each sample/condition/path combo its own string
sleuth_input = [] #initialize empty list
for sample, condition, path in zip(samples, conditions, paths):  #parse through list parallelly
    final_string = (f"{sample}\t{condition}\t{path}")  #format the string
    sleuth_input.append(final_string)  #add the formatted string to the list

#format header
header = ("sample"+"\t"+"condition"+"\t"+"path")

#write header and all cample/condition/path combos to log file
with open((f'{dir}/sleuth_input.txt'), 'w') as f:
    f.write(header + "\n")
    for group in sleuth_input:
        f.write(group + "\n")