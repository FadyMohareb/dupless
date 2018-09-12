#!/bin/python

from subprocess import call, Popen, PIPE
from Bio import SeqIO
import sys
import os

def usage():
    print("python filter_blast_results.py blast_results.tab min_identity min_length assembly.fasta output_folder")


def get_assembly_coordinates_from_blast_results(region_name, region_start, region_stop, assembly_dict):
    """
    From a blast result "scaffold_start_stop start stop"
    Returns:
        The scaffold name in the assembly ("scaffold_start_stop" becomes "scaffold").
        The start position of the blast hit in the assembly "start het region + start blast hit"
        The stop position of the blast hit in the assemly "stop het region + stop blast hit"
    """
    name_split = region_name.split("_")
    # If the original scaffold name contained underscores, we put them back with join
    scaffold_name = "_".join(name_split[:-2])
    scaffold_start = int(name_split[-2])
    scaffold_end = int(name_split[-1])

    scaffold_region_start = int(region_start) + scaffold_start
    scaffold_region_end = int(region_stop) + scaffold_end

    return scaffold_name, scaffold_region_start, scaffold_region_end, len(assembly_dict[scaffold_name])


def check_file(filename):
    """
    Checks if "filename" exists and is a file.
    Returns:
        True if file exists and is a file.
        False if filename==None or is not a file.
    """
    file_ok = True
    error_mssg = ""
    if(filename == None):
        error_mssg = "Error: file '"+filename+"' is missing."
        file_ok = False
    else:
        if not os.path.isfile(filename):
            error_mssg = "Error: '"+filename+"' is not a file."
            file_ok = False
    return file_ok, error_mssg

#=================================================================
#                             Input checks                       =
#=================================================================
# Only 5 arguments needed, but length of 6 because sys.argv[0] = name of the script
if(len(sys.argv) < 6):
    usage()
    print("5 arguments expected, only "+str(len(sys.argv)-1)+ " found")
    sys.exit(2)

blast_filename = sys.argv[1]
blast_identity_threshold = float(sys.argv[2])
blast_length_threshold = float(sys.argv[3])
assembly = sys.argv[4]
output_folder = sys.argv[5]


file_ok, error_mssg =  check_file(blast_filename)
if not(file_ok):
    print("Blast output file could not be found. See error message below for more details:")
    print(error_mssg)
    sys.exit(2)

file_ok, error_mssg = check_file(assembly)
if not(file_ok):
    print("Assembly fasta could not be found. See error message below for more details:")
    print(error_mssg)
    sys.exit(2)

#=================================================================
#                          Main                                  =
#=================================================================
print("Creating dictionnary for the assembly...")
assembly_dict = SeqIO.index(assembly, "fasta")
print("Done !")

if not os.path.isdir(output_folder):
    print("Creating the output directory...")
    process = Popen(["mkdir", output_folder], stdout=PIPE)
    process.wait()

print("\nFiltering blast results...")
toRemovehap1_name = output_folder+"/toRemoveFromhap1.bed"
process = Popen(["touch", toRemovehap1_name], stdout=PIPE)
process.wait()
toRemovehap2_name = output_folder+"/toRemoveFromhap2.bed"
process = Popen(["touch", toRemovehap2_name], stdout=PIPE)
process.wait()

to_keep = list()
with open(blast_filename, "r") as all_blasts:
    # For each duplicated pair, we write one to hap1 and the other to hap2
    with open(toRemovehap1_name, "a") as hap1, open(toRemovehap2_name, "a") as hap2:
        for blast in all_blasts:
            tabs = blast.split("\t")
            if((float(tabs[2]) >= blast_identity_threshold) and
                (float(tabs[3]) >= blast_length_threshold)):
                # If "S1 S2" is a valid blast hit, we don't want to eliminate S2 S1 later (double removal),
                # So we add the reversed pair "S2;S1" in to_keep
                pair = tabs[0]+";"+tabs[1]
                reversed_pair = tabs[1]+";"+tabs[0]
                if(pair not in to_keep):
                    # The QUERY should always be "start stop", but if there is a inverted alignment, the SUBJECT can be "stop   start".
                    # So we get the start by taking the min and the stop with the max
                    query_start = int(tabs[6])
                    query_end = int(tabs[7])

                    subject_start = min(int(tabs[8]), int(tabs[9]))
                    subject_end = max(int(tabs[8]), int(tabs[9]))

                    # We need to get the original positions, since "query_start" and "query_end" correspond to the het regions positions and not the scaffold position.
                    q_scaffold_name, q_scaffold_start, q_scaffold_end, q_length = get_assembly_coordinates_from_blast_results(tabs[0],query_start, query_end, assembly_dict)
                    s_scaffold_name, s_scaffold_start, s_scaffold_end, s_length = get_assembly_coordinates_from_blast_results(tabs[1], subject_start, subject_end, assembly_dict)

                    # Instead of randomly removing heterozygous regions, we remove regions from the smallest contig in hap1 and from the longest in hap2.
                    # This will reduce the chance of cutting a gene in half and hence downgrading hap1 assembly quality. Hap1's N50 will also be less affected.
                    if(q_length < s_length):
                        hap1.write(q_scaffold_name+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\n")
                        hap2.write(s_scaffold_name+"\t"+str(s_scaffold_start)+"\t"+str(s_scaffold_end)+"\n")
                    else:
                        hap1.write(s_scaffold_name+"\t"+str(s_scaffold_start)+"\t"+str(s_scaffold_end)+"\n")
                        hap2.write(q_scaffold_name+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\n")

                    to_keep.append(reversed_pair)
print("Done !")
