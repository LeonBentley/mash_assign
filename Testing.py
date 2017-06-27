#Imports relevant stuff
import os,sys
import argparse
import subprocess
import itertools
import random

#Imports functions from a local script file
import assign_functions

separator = "\t"

# Read in options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

io.add_argument("--inputs_num", dest="inputs_num", help="Numbers of samples to be assigned", default="10")
io.add_argument("--cluster_file", help="Cluster label file")
io.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
io.add_argument("-o","--output", dest="output_file", help="Output file", default="accuracies.txt")
options.add_argument("--repeats", dest="repeats", help="Numbers of times whole test should be run",default='10')
options.add_argument("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
options.add_argument("--kmer_size", dest="kmer_size", help="K-mer size for mash sketches", default="21")
options.add_argument("--sketch_size", dest="sketch_size", help="Size of mash sketch", default="1000")
options.add_argument("--max_kmer_size", dest="max_kmer_size", help="Max K-mer size for mash sketches", default="1000")
options.add_argument("--min_kmer_size", dest="min_kmer_size", help="Minimum K-mer size for mash sketches", default="3")
options.add_argument("--max_sketch_size", dest="max_sketch_size", help="Max sketch size for mash sketches", default="100000")
options.add_argument("--min_sketch_size", dest="min_sketch_size", help="Minimum sketch size for mash sketches", default="500")
args = parser.parse_args()

#Checks arguments are intergers in sensible range
if not args.kmer_size.isdigit:
    sys.stderr.write(str(args.kmer_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.kmer_size) < int(args.min_kmer_size) or int(args.kmer_size) > int(args.max_kmer_size):
    sys.stderr.write(str(args.kmer_size) + " is an invalid value\n")
    sys.exit(1)

if not args.sketch_size.isdigit:
    sys.stderr.write(str(args.sketch_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.sketch_size) < int(args.min_sketch_size) or int(args.sketch_size) > int(args.max_sketch_size):
    sys.stderr.write(str(args.sketch_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.inputs_num) < 0 or int(args.inputs_num) > 2000:
    sys.stderr.write(str(args.inputs_num) + " is an invalid value\n")
    sys.exit(1)

if not args.inputs_num.isdigit:
    sys.stderr.write(str(args.inputs_num) + " is an invalid value\n")
    sys.exit(1)

# Check input files exist
if not os.path.isfile(str(args.cluster_file)):
    if args.cluster_file == None:
        sys.stderr.write("No cluster file entered\n")
        sys.exit(1)
    else:
        sys.stderr.write(str(args.cluster_file) + " does not exist\n")
        sys.exit(1)

if not os.path.isfile(str(args.assembly_file)):
    if args.assembly_file == None:
        sys.stderr.write("No assembly file entered\n")
        sys.exit(1)
    else:
        sys.stderr.write(str(args.assembly_file) + " does not exist\n")
        sys.exit(1)

inputs = assign_functions.input_reader(args.cluster_file, args.assembly_file)#Runs the function with provided inputs and makes a list of its outputs

assign_functions.reference_creator(inputs[1], args.kmer_size, args.sketch_size, args.mash_exec)#Runs the function and provides inputs

accuracy_list = []#Creates empty list

for count in range(0, int(args.repeats)):
    #Takes 'k' random files from provided set
    items = inputs[1]
    k = int(args.inputs_num)
    random_input = items[0:k]

    for i in range(k, len(items)):
        j = random.randrange(0, i)
        if j < k:
            random_input[j] = items[i]

    accuracycount = 0#Sets variable as 0

    for random_sample in random_input:#Takes first value in list and names in randon_sample
        assign_functions.sample_sketching (args.mash_exec, args.kmer_size, args.sketch_size, random_sample)#Runs the function and provides inputs
        cluster_output = assign_functions.cluster_identification(args.mash_exec, inputs[0] , inputs[2], random_sample, random_input)#Runs the function with provided inputs and makes a list of its outputs

        #Prints outputs with coloured text to make manual testing easier
        print("Assigned Cluster is " + str(cluster_output[2][0]))
        print("Actual Cluster is " + str(inputs[2][inputs[0][random_sample]]))
        if int(cluster_output[2][0]) == int(inputs[2][inputs[0][random_sample]]):
            print('\x1b[6;30;42m' + 'Correct' + '\x1b[0m')
            accuracycount = accuracycount + 1

        else:
            print('\033[6;30;91m' + 'Incorrect' + '\033[00m')

        print(count + '/' + (inputs_num))

    #Calculates and prints accuracy, appends accuracy to list
    accuracy = accuracycount/float(k)
    print(accuracy)
    accuracy_list.append(str(accuracy))

#Writes list to file and prints
with open(str(args.accuracies_file), "w+") as file:
        file.write(','.join(accuracy_list) + '\n')#Joins different values in list with seperator ','
print(accuracy_list)