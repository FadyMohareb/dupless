#!/usr/bin/python

# Main DupLess script

# Dependencies:
#   python v2.7 or higher
#   samtools v1.9 or higher (important for the "-o" parameter)
#   bedtools v2.27 (lower version should now work)
#   blastn v2.6.0+
#   pandas, numpy, matplotlib, multiprocessing, getopt, biopython, sys, os, subprocess
#   sed and awk

import getopt
import subprocess
import sys
import os

import detect_het_regions_from_coverage as dh
import detect_duplicates_from_het_regions as dd
import processing as proc
import utils_dupless as ud


global VERSION
VERSION = "1.0.0"

def print_version():
    """Print the version.
    """
    global VERSION
    print("DupLess v" + VERSION)


def usage():
    """Print the usage.
    """
    print("\npython DupLess.py -t [nb_threads] -b [coverage.bed] -a [assembly.fasta] -w [window_size] -c [expected_coverage] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]")
    print("\nRequired:")
    print("     -a/--assembly               The assembly corresponding to the bed coverage in fasta format.")
    print("")
    print("     -b/--bed_cov                The bed file containing the coverage at each base (can be generated with 'bedtools genomecov').")
    print("                                 /!\ If using paired end reads: make sure that you set the -w or -l option higher than the insert size,")
    print("                                     to avoid false positives due to coverage drop at the ends of contigs (because of unaligned mates).")
    print("\nOptional:")
    print("     -t/--nThreads               The number of threads (default 10)")
    print("     -o/--out_folder             The output folder (default './DupLess_out/')")
    print("")
    print("     -c/--expected_cov           The expected read coverage for the homozygous regions. The homozygosity / heterozygosity will be determined based on this value.")
    print("                                 You can determine the value to use by plotting the coverage distribution. It should correspond to the homozygous peak")
    print("                                 If no value is given, it will be based on the mode of the coverage distribution (not reliable if high heterozygosity).")
    print("")
    print("     -w/--window_size            The size of the windows in basepairs (default: 1000)")
    print("                                 The value of the coverage for each window will be the median of the coverage at each base.")
    print("                                 All the windows classified as 'heterozygous' will be considered for the detection of duplication.")
    print("")
    print("     -g/--bed_gaps               A bed file containing the gaps along the genome. If given, the graphs will contain a grey background where the gaps are.")
    print("")
    print("     -i/--blast_identity         The minimum percentage of identity between the het regions to consider them duplicates (default: 90, range 0 to 100).")
    print("     -l/--blast_length           The blast alignments with a length lower than this threshold will be filtered (default=300).")
    print("")
    print("     -n/--no_plot                Skip the creation of all the plots")
    print("\nSkipping part of pipeline:")
    print("     -s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed identifying the heterozygous regions:")
    print("                                      python DupLess.py -s [het_regions_bed] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]")
    print("")
    print("     -f/--filter_blast_only      Skip the detection of the heterozygous regions AND the pairwise alignments. If so, you must provide a blast ouput with -oufmt 6:")
    print("                                      python DupLess.py -f [blast_output] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]")
    print("\nOther:")
    print("     -h/--help                   Print the usage and help and exit.")
    print("     -v/--version                Print the version and exit.")


#=================================================================
#                           GetOpt                               =
#=================================================================
window_size =  1000             # The coverage of each window will be based on the median of all the coverages inside the window.
coverage_bed = None             # Bed file with the coverage value for each position. Can be produced with "bedtools coverage".
assembly_name = None            # Assembly in fasta format, used to extract the het regions and also check the scaffold lengths.
expected_coverage = None        # Any window with 0 < coverage < expected_cov/1.5 will be considered as heterozygous.
gaps_bed = None                 # Optional. Used to draw gaps as grey bars on the coverage graphs.
output_folder = "./DupLess_out/"
nbThreads = 10
blast_identity_threshold = 90   # Two regions will be considered duplicated if...
blast_length_threshold = 300      # ...these two blast thresholds are met (min identity and min length).
skip_het_dect = False           # Possibility to skip the first step (het detection) if bed of heterozygous regions is provided.
het_bed = None                  # Bed defining the heterozygous region (Created by DupLess, or given by the user if skip_het_dect = T).
skip_blast = False              # Possibility to skip the "het detection" and "pairwise blasting" and just filter the blast results.
blast_output = None             # Default output file for blast results (given by the user if skip_blast = T)
skip_plot = False               # Skip the generation of the coverage plots

try:
    opts, args = getopt.getopt(sys.argv[1:], "t:w:b:a:c:g:o:s:f:i:l:nhv", ["nThreads=", "window_size=", "bed_cov=", "assembly=", "expected_cov=", "bed_gaps=", 
                                                                            "out_folder=", "skip_het_detection=", "filter_blast_only=",
                                                                            "blast_identity=", "blast_length=", "no_plot", "help", "version"])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

for o,a in opts:
    if o in ("-t", "--nThreads"):
        nbThreads = int(a)
    elif o in ("-w", "--window_size"):
        window_size = int(a)
    elif o in ("-b", "--bed_cov"):
        coverage_bed = str(a)
    elif o in ("-a", "--assembly"):
        assembly_name = str(a)
    elif o in ("-c", "--expected_cov"):
        expected_coverage = int(a)
    elif o in ("-g", "--bed_gaps"):
        gaps_bed = str(a)
    elif o in ("-o", "--out_folder"):
        output_folder = str(a)
    elif o in ("-s", "--skip_het_dect"):
        het_bed = str(a)
        skip_het_dect = True
    elif o in ("-f", "--filter_blast_only"):
        blast_output = str(a)
        skip_het_dect = True
        skip_blast = True
    elif o in ("-i", "--blast_identity"):
        blast_identity_threshold = int(a)
    elif o in ("-l", "--blast_length"):
        blast_length_threshold = float(a)
    elif o in ("-n", "--no_plot"):
        skip_plot = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit(1)
    elif o in ("-v", "--version"):
        print_version()
        sys.exit(1)
    else:
        assert False, "Unhandled option !"


# If we do not skip the het detection step:
if not skip_het_dect:
    # Then we need the coverage bed
    file_ok, error_mssg = ud.check_file(coverage_bed)
    if not file_ok:
        print("Error with option -b/--bed_cov: "+error_mssg)
        usage()
        sys.exit(2)

# Exit if output folder already exists (to avoid overwriting an already existing project)
if(os.path.isdir(output_folder)):
    print("\nFolder '"+output_folder+"' already exists, stopping now...\n")
    sys.exit(2)

file_ok, error_mssg = ud.check_file(assembly_name)
if not file_ok:
    print("Error with option -a/--assembly: "+error_mssg)
    usage()
    sys.exit(2)

if(window_size <= 0):
    print("The window size can not be lower than 0 (-w/--window_size option).\n")
    usage()
    sys.exit(2)

if(nbThreads <= 0):
    print("The number of threads can not be lower than 0 (-t/--nThreads option).\n")
    usage()
    sys.exit(2)

if((blast_identity_threshold < 0) or (blast_identity_threshold > 100)):
    print("The blast identity treshold (-i/--blast_identity option) can not be higher than 100 or lower than 0. Current value: "+str(blast_identity_threshold)+"\n")
    usage()
    sys.exit(2)

if((blast_length_threshold < 0)):
    print("The blast coverage treshold (-l/--blast_length option) can not be lower than 0. Current value: "+str(blast_length_threshold)+"\n")
    usage()
    sys.exit(2)


#=================================================================
#                          Main                                  =
#=================================================================
# Creating the output folder architecture
# individuals_beds/  contains the bed files describing the het. regions for each sequence.
# invidual_blasts/   contains the blast results for each het. region.
# graphs/            contains the coverage graphs for each sequence.
# temp/              contains temp file for blast.
# deduplicated/      contains the results of DupLess: deduplicated.fasta and discarded.fasta
for folder in [output_folder, output_folder+"/individual_beds", output_folder+"/graphs", output_folder+"/individual_blasts", output_folder+"/temp", output_folder+"/deduplicated"]:
    try:
        pr = subprocess.Popen(["mkdir", folder], shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, "mkdir "+folder)
    except Exception as e:
        print("Error during mkdir "+folder)
        print("Exception:"+str(e))
        sys.exit()


# Indexing the assembly, needed later on for extraction of het regions
# Also a good way to check if samtools exists at the start of the script
ud.index_fasta_file(assembly_name)

#==================================================
#   Detection of the heterozygous regions
#==================================================
if not skip_het_dect:
    # If the user does not skip the het dect then we need the coverage bed file
    # We check if it exists
    file_ok, error_mssg = ud.check_file(coverage_bed)
    if file_ok:
        # Launch the bed and graph creation for heterozygous regions, detection based on coverage values.
        het_bed = dh.detect_het_regions(coverage_bed, gaps_bed, expected_coverage, window_size, output_folder, nbThreads, skip_plot)
    else:
        print("Error with the coverage bed file: "+error_mssg)
        usage()
        sys.exit(2)


#==================================================
#   Pairwise alignment of heterozygous regions
#==================================================
if not skip_blast:
    # Check if the het_bed file exists (wether it has been created by the step before or just given by the user)
    file_ok, error_mssg = ud.check_file(het_bed)
    if file_ok:
        # Launch pairwise blast comparisons between the detected heterozygous regions to detect duplication
        blast_output = dd.detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads)
    else:
        print("Error with the heterozygous bed file: "+error_mssg)
        usage()
        sys.exit(2)


#==================================================
#   Filtering the blast results
#==================================================
# Check if the blast output file exists (wether it has been created by the step before or just given by the user)
file_ok, error_mssg = ud.check_file(blast_output)
if file_ok:
    # Filter the blasts by identity and length
    toRemoveBed, discardedBed = dd.filter_blast_results(blast_output, blast_identity_threshold, blast_length_threshold, assembly_name, output_folder)
else:
    print("Error with the blast output file: "+error_mssg)
    usage()
    sys.exit(2)


#==================================================
#   Generating the output files
#==================================================
deduplicated_assembly = output_folder+"/deduplicated/deduplicated_assembly.fasta"
discarded_assembly = output_folder+"/deduplicated/discarded.fasta"

print("Generating the deduplicated fasta files from the blast results...")
proc.remove_duplications_assembly(deduplicated_assembly, assembly_name, toRemoveBed, output_folder)
proc.generate_discarded_fasta(assembly_name, discardedBed, discarded_assembly)

# If we skip blast, no intermediate files created, no need to clean
if not skip_blast:
    # Cleaning the intermediate files:
    ud.remove_file(output_folder+"/All_Blasts_region_coord.tab")
    ud.remove_file(output_folder+"/assembly_HET_ONLY.fa")
    ud.remove_file(output_folder+"/assembly_HET_ONLY.fa.fai")

print("Done !\n")
print("Deduplicated assembly generated in: " + deduplicated_assembly)
print("Discarded sequences in: " + discarded_assembly)
