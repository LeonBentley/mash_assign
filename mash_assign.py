import os,sys
import argparse
import subprocess
import itertools

import assign_functions

separator = "\t"

# Read in options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

io.add_argument("--cluster_file", help="Cluster label file")
io.add_argument("--input_file", help="Fasta file to assign")
io.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("--output_count", dest="output_count", help="Numbers of samples outputted", default="10")
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

# Check input files exist
if not os.path.isfile(str(args.input_file)):
    if args.input_file == None:
        sys.stderr.write("No input file entered\n")
        sys.exit(1)
    else:
        sys.stderr.write(str(args.input_file) + " does not exist\n")
        sys.exit(1)

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

inputs = assign_functions.input_reader(args.cluster_file, args.assembly_file)

assign_functions.reference_creator(inputs[1], args.kmer_size, args.sketch_size, args.mash_exec)

assign_functions.sample_sketching (args.mash_exec, args.kmer_size, args.sketch_size, args.input_file)

cluster_output = assign_functions.cluster_identification(args.mash_exec, inputs[0] , inputs[2], args.input_file, None, args.output_count)
for y in range(0, len(cluster_output[0])):
    output_list = [str(cluster_output[1][y]), str(cluster_output[0][y]), str(cluster_output[2][y])]
    output_list = '\t'.join(output_list)
    print(output_list)