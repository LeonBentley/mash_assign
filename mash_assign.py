#Imports relevant stuff
import os,sys
import argparse
import subprocess
import itertools

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
args = parser.parse_args()

#Checks arguments are intergers in sensible range
if not args.kmer_size.isdigit:
    sys.stderr.write(str(args.kmer_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.kmer_size) < 3 or int(args.kmer_size) > 1000:
    sys.stderr.write(str(args.kmer_size) + " is an invalid value\n")
    sys.exit(1)

if not args.sketch_size.isdigit:
    sys.stderr.write(str(args.sketch_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.sketch_size) < 500 or int(args.sketch_size) > 100000:
    sys.stderr.write(str(args.sketch_size) + " is an invalid value\n")
    sys.exit(1)

if int(args.inputs_num) < 0 or int(args.inputs_num) > 2000:
    sys.stderr.write(str(args.inputs_num) + " is an invalid value\n")
    sys.exit(1)

if not args.inputs_num.isdigit:
    sys.stderr.write(str(args.inputs_num) + " is an invalid value\n")
    sys.exit(1)

inputs = assign_functions.input_reader(args.cluster_file, args.assembly_file)

assign_functions.reference_creator(inputs[1], args.kmer_size, args.sketch_size, args.mash_exec)

assign_functions.sample_sketching (args.mash_exec, args.kmer_size, args.sketch_size, args.input_file)

cluster_output = assign_functions.cluster_identification(args.mash_exec, inputs[0] , inputs[2], args.input_file, None, args.output_count)

#Formats the output (purely aesthetic) and then prints
for y in range(0, len(cluster_output[0])):
    output_list = [str(cluster_output[1][y]), str(cluster_output[0][y]), str(cluster_output[2][y])]
    output_list = '\t'.join(output_list)
    print(output_list)