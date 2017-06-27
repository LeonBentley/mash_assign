import os,sys
import argparse
import subprocess
import itertools
import operator

separator = "\t"

def input_reader(cluster_file, assembly_file):
    # Read input files
    # Read figInfo.txt
    clusters = {}#creates dictonary
    with open(str(cluster_file)) as f:
        line_num = 0
        for line in f:
            line_num = line_num + 1
            if line_num > 2:
                line = line.rstrip()#removes \n
                (row, cluster, sampleID) = line.split(separator)#splits into 3 variables
                if row.isdigit:#checks line is not header
                    clusters.update({sampleID:cluster})#adds key and item, unless they already exists

    location_dict = {}#creates empty dictionary
    file_names = []#creates empty list
    with open(str(assembly_file)) as f:
        for line in f:
            line = line.rstrip()#removes \n
            (sampleID, file_name) = line.split(separator)#splits into 2 variables
            file_names.append(file_name)#adds address/location
            location_dict.update({file_name:sampleID})#adds new key and item

    return location_dict, file_names, clusters
    
def reference_creator(file_names, kmer_size, sketch_size, mash_exec):
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

            mash_command = str(mash_exec) + " sketch -k " + kmer_size + " -s " + sketch_size + " -o reference" + str(chunk) + " " + mash_sep.join(file_names[chunk_start:chunk_end])
            retcode = subprocess.call(mash_command, shell=True)#makes sketch for assembly file
            if retcode != 0:#checks for code failure
                sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
                sys.exit(1)

        if ((len(file_names) - 1) // mash_chunk_size) > 0:#checks if more than one chunk
            paste_join = ".msh reference"
            mash_paste_command = str(mash_exec) + " paste reference reference" + paste_join.join([str(chunk) for chunk in range(num_chunks)]) + ".msh"
            retcode = subprocess.call(mash_paste_command, shell=True)#runs command to put all chunks together
            if retcode != 0:#checks for code failure
                sys.stderr.write("Mash paste failed with signal " + str(retcode) + "\n")
                sys.exit(1)
            else:
                for chunk in range(num_chunks):
                    os.remove("reference" + str(chunk) + ".msh")#deletes the temporary files used for chunks
        else:
            os.rename("reference0.msh", "reference.msh")#if there are less clusters than the size of one chunk, only one file is made, named without the number

def cluster_identification (mash_exec, location_dict, clusters, msh_file_name, random_testing = None, output_count = 1):

    #Establishes variables, lists and dictonary
    distances = {}
    min_clusters = []
    min_IDs = []
    min_dists = []
    random_IDs = []

    #Uses list of random files to make a random list of IDs
    if random_testing != None:
        for random_sample in random_testing:
            random_IDs.append(location_dict[random_sample])

    p = subprocess.Popen([str(mash_exec) + ' dist reference.msh ' + str(msh_file_name)], stdout=subprocess.PIPE, shell=True)#Creates command to calculate distance between two files
    for line in iter(p.stdout.readline, ''):#Reads all lines in output of command above
        line = line.decode('utf-8')#Decodes line to specific encoding
        line = line.rstrip()#Removes \n from end of line

        if line != '':
            (name1, name2, dist, p, matches) = line.split(separator)#Splits line at tabs, into 5 variables
            distances.update({location_dict[name1]:dist})#Adds two newly created variables to dictionary

        else:
            break

    sorted_distances = sorted(distances.items(), key = operator.itemgetter(1))#Makes list of dictionary contents (tuples) ordered by value

    #Adds contents derived from list of tuples to lists to be used in final output
    count = 0
    for x in range(0, len(sorted_distances)):
        if random_testing == None or sorted_distances[x][0] not in random_IDs:
            if count < int(output_count):#Determines when loop ends
                min_dists.append(sorted_distances[x][1])
                min_IDs.append(sorted_distances[x][0])
                min_clusters.append(clusters[sorted_distances[x][0]])
                count += 1
            else:
                break

    return(min_dists, min_IDs, min_clusters)#outputs relevant lists

def sample_sketching(mash_exec, kmer_size, sketch_size, input_file):
    (split, apart) = input_file.split(".")
    stitch = split + '.msh'
    if not os.path.isfile(stitch):
        mash_command = str(mash_exec) + " sketch -k " + kmer_size + " -s " + sketch_size + " -o input " + str(input_file)
        retcode = subprocess.call(mash_command, shell=True)#makes sketch for new sample
        if retcode != 0:#checks for code failure
            sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
            sys.exit(1)
    return(stitch)#Outputs mash file of sample to be assigne