import os,sys
import argparse
import subprocess
import itertools

separator = "\t"

# Read in options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

io.add_argument("--cluster_file", help="Cluster label file")
io.add_argument("--input_file", help="Fasta file to assign")
io.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
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

# Read input files
# Read figInfo.txt
clusters = {}#creates dictonary
with open(str(args.cluster_file)) as f:
    line_num = 0
    for line in f:
        line_num = line_num + 1
        if line_num > 2:
            line = line.rstrip()#removes \n
            (row, cluster, sampleID) = line.split(separator)#splits into 3 variables
            if row.isdigit:#checks line is not header
                clusters.update({sampleID:cluster})#adds key and item, unless they already exists

file_names = []#creates empty list
sampleIDs = []#creates empty list
with open(str(args.assembly_file)) as f:
    for line in f:
        line = line.rstrip()#removes \n
        (sampleID, file_name) = line.split(separator)#splits into 2 variables
        sampleIDs.append(sampleID)#adds sampleID
        file_names.append(file_name)#adds address/location

location_dict = {}#creates empty dictionary
with open(str(args.assembly_file)) as f:
    for line in f:
        line = line.rstrip()#removes \n
        (sampleID, file_name) = line.split(separator)#splits into 2 variables
        location_dict.update({file_name:sampleID})#adds new key and item


#If there is no reference.msh, one is created
#Splits into managable chunks(500)
mash_chunk_size = 500
num_chunks = ((len(file_names) - 1) // mash_chunk_size) + 1
if not os.path.isfile("reference.msh"):#checks if file already exists
    mash_sep = " "#seperator
    for chunk in range(num_chunks):
        chunk_start = chunk * mash_chunk_size#sets chunk start
        chunk_end = (chunk+1) * mash_chunk_size#and chunk end
        if chunk_end > len(file_names):#if chunk end is beyond the limits of the data provided...
            chunk_end = len(file_names)#chunk end is set according to the length of the data

    mash_command = str(args.mash_exec) + " sketch -k " + args.kmer_size + " -s " + args.sketch_size + " -o reference" + str(chunk) + " " + mash_sep.join(file_names[chunk_start:chunk_end])
    retcode = subprocess.call(mash_command, shell=True)#makes sketch for assembly file
    if retcode != 0:#checks for code failure
        sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
        sys.exit(1)

    if ((len(file_names) - 1) // mash_chunk_size) > 0:#checks if more than one chunk
        paste_join = ".msh reference"
        mash_paste_command = str(args.mash_exec) + " paste reference reference" + paste_join.join([str(chunk) for chunk in range(num_chunks)]) + ".msh"
        retcode = subprocess.call(mash_paste_command, shell=True)#runs command to put all chunks together
        if retcode != 0:#checks for code failure
            sys.stderr.write("Mash paste failed with signal " + str(retcode) + "\n")
            sys.exit(1)
        else:
            for chunk in range(num_chunks):
                os.remove("reference" + str(chunk) + ".msh")#deletes the temporary files used for chunks
    else:
        os.rename("reference0.msh", "reference.msh")#if there are less clusters than the size of one chunk, only one file is made, named without the number

mash_command = str(args.mash_exec) + " sketch -k " + args.kmer_size + " -s " + args.sketch_size + " -o input " + (args.input_file)
retcode = subprocess.call(mash_command, shell=True)#makes sketch for new sample
if retcode != 0:#checks for code failure
    sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
    sys.exit(1)

min_dist = 1
p = subprocess.Popen([str(args.mash_exec) + ' dist reference.msh input.msh'], stdout=subprocess.PIPE, shell=True)
for line in iter(p.stdout.readline, ''):
    line = line.decode('utf-8')
    line = line.rstrip()
    if line != '':
        (name1, name2, dist, p, matches) = line.split(separator)
        if float(dist) < min_dist:
            min_dist = float(dist)
            location_file = name1
    else:
        break

print("Minimum distance is " + str(min_dist))
print("SampleID is " + str(location_dict[location_file]))
print("Cluster is " + str(clusters[(location_dict[location_file])]))