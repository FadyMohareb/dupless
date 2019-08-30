#!/usr/bin/python

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
from Bio import SeqIO
import pandas as pd
import re

import detect_het_regions_from_coverage as dh
import detect_duplicates_from_het_regions as dd
import utils_dupless as ud


global VERSION
VERSION = "1.0.0"

def print_version():
    """Print the version.
    """
    global VERSION
    print("DupLess v"+VERSION)


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



def phasing(toRemoveBed, blast_output, outname, assembly_name, output_folder):
    """Reading the heterozygous regions blast tab and the toRemoveBed to locate to what subject contigs the queries (recorded in the Bed) have matched.
    Therefore, we will generate a second bed containing subject sequences instead of query sequences, and the removal of these will be the second haplotype.
    """
    file_ok, error_mssg =  ud.check_file(blast_output)
    if not(file_ok):
        print("Blast output file could not be found. See error message below for more details:")
        print(error_mssg)
        sys.exit(2)

    file_ok, error_mssg = ud.check_file(assembly_name)
    if not(file_ok):
        print("Assembly fasta could not be found. See error message below for more details:")
        print(error_mssg)
        sys.exit(2)
    
    toRemoveBed1 = output_folder+"/toDiscard_hap1.bed"
    process = subprocess.Popen(["touch", toRemoveBed1], stdout=subprocess.PIPE)
    process.wait()    
    toRemoveBed2 = output_folder+"/toDiscard_hap2.bed"
    process = subprocess.Popen(["touch", toRemoveBed2], stdout=subprocess.PIPE)
    process.wait()
    repetitiveBed = output_folder+"/repetitive.bed"
    process = subprocess.Popen(["touch", repetitiveBed], stdout=subprocess.PIPE)
    process.wait()
    
    # Needed to quickly get the lengths of the sequences
    assembly_dict = SeqIO.index(assembly_name, "fasta")
    
    to_keep = list()
    blast_df = pd.DataFrame(columns=['query_name','subject_name','query_start','query_end','subject_start','subject_end','q_length','s_length','pair','reversed_pair','subject_is_reversed'])
    with open(blast_output, "r") as het_blasts:
        # For each duplicated pair, we write the blast values to a dataframe
        for blast in het_blasts:
            tabs = blast.split("\t")
            # tabs[2] contains the blast identity
            # tabs[3] contains the length of the blast hit
            # For the output format 6 (which is used by DupLess)
            if((float(tabs[2]) >= blast_identity_threshold) and (float(tabs[3]) >= blast_length_threshold)):
                query_name = tabs[0]
                subject_name = tabs[1]
                # The QUERY should always be "start stop", but if there is a inverted alignment, the SUBJECT can be "stop   start".
                # So we get the start by taking the min and the stop with the max
                query_start = int(tabs[6])
                query_end = int(tabs[7])
                subject_start = min(int(tabs[8]), int(tabs[9]))
                subject_end = max(int(tabs[8]), int(tabs[9]))
                subject_is_reversed = True if int(tabs[8])>int(tabs[9]) else False
              
                q_length = len(assembly_dict[query_name])
                s_length = len(assembly_dict[subject_name])

                # If "S1 S2" is a valid blast hit, we don't want to eliminate "S2 S1" later (double removal),
                # so we will add the reversed pair "S2_start_stop;S1_start_stop" in to_keep.
                # We add the start and stop in case there are several hits between different regions on the same two scaffolds.
                region1 = query_name+"_"+str(query_start)+"_"+str(query_end)
                region2 = subject_name+"_"+str(subject_start)+"_"+str(subject_end)
                pair = region1+";"+region2
                reversed_pair = region2+";"+region1
                
                #Adding all parsed blast variables to a pandas dataframe row
                blast_df.loc[len(blast_df)]=[query_name,subject_name,query_start,query_end,subject_start,subject_end,q_length,s_length,pair,reversed_pair,subject_is_reversed]
    
    
    
    #Using kat sect to count kmers to find repetitive regions in the original assembly
    popen_args=["kat","sect","-E","-F","-t",str(nbThreads),"-M",str(4),assembly_name,assembly_name, "-v"]
    process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
    process.wait()
    #Reading the fasta with the repetitive regions and extracting their sequence locations
    repetitive_region_headers=[x for x in open("kat-sect-repetitive.fa").readlines() if x.startswith(">")]
    repetitive_df = pd.DataFrame(columns=['contig','start','end'])
    for x in repetitive_region_headers:
      contig=x.split("__")[0][1:]
      contig_positions = re.search('pos:(.*)_cov', x).group(1)
      start=int(contig_positions.split(":")[0])
      end=int(contig_positions.split(":")[1])
      repetitive_df.loc[len(repetitive_df)]=[contig,start,end]
    
    
    
    
    
    
    ########################################## 
    #must also check for the position in the contig of the hit, if its start/end or not
    #Checking for duplicate rows based on the name, start and end of the query:(maybe should also check if there are internal blasts and much more, there is a lot of detail here
    multipleHits = blast_df[blast_df.duplicated(['query_name','query_start','query_end'])] 
    #if there have been hits it is time to treat them accordingly
    if not multipleHits.empty:
      print('there are multiple hits')
    ############################################################################   
        
    
    
    
    
    
    
    
    #Writing duplicated bed files, and eliminating the reverse pair of the written pairs so as to not appear twice
    for index, row in blast_df.iterrows():
        with open(toRemoveBed1, "a") as toRemove1Handle,open(toRemoveBed2, "a") as toRemove2Handle:
            if(row['pair'] not in to_keep):
                toRemove1Handle.write(str(row['query_name'])+"\t"+str(row['query_start'])+"\t"+str(row['query_end'])+"\n")#"random" haplotype generation based on query/subject
                toRemove2Handle.write(str(row['subject_name'])+"\t"+str(row['subject_start'])+"\t"+str(row['subject_end'])+"\n")#without looking at lengths
                to_keep.append(row['reversed_pair'])                        
    #Writing bed file of the kat sect repetitive regions fasta  
    for index, row in repetitive_df.iterrows():
        with open(repetitiveBed, "a") as repetitiveHandle:
          repetitiveHandle.write(str(row['contig'])+"\t"+str(row['start'])+"\t"+str(row['end'])+"\n")
     
    #Subtracting repetitive features from the discard beds         
    popen_args=["bedtools","subtract","-a",toRemoveBed1,"-b",repetitiveBed]
    with open("hap1.bed","wb") as out:
      process = subprocess.Popen(popen_args, stdout=out)   
    popen_args=["bedtools","subtract","-a",toRemoveBed2,"-b",repetitiveBed]
    with open("hap2.bed","wb") as out:
      process = subprocess.Popen(popen_args, stdout=out)
          
                
    ################################################################################    
    remove_duplications_assembly("hap1.fasta", assembly_name, "hap1.bed", output_folder)
    remove_duplications_assembly("hap2.fasta", assembly_name, "hap2.bed", output_folder)  
    
    

def remove_duplications_assembly(outname, assembly_name, bedname, output_folder):
    """From an assembly and a bed of regions to remove, create a fasta with duplications removed.
    The regions in the bed are masked as "$" with "bedtools maskfasta", then sed is used to replace "$" by "".
    Also use "awk" to have each sequence on one line (to avoid creating empty lines with sed).
    """
    # Temp files used by sed and awk, removed at the end of the function
    fasta_masked = output_folder+"/deduplicated/temp_masked.fasta"
    fasta_masked_oneLine = output_folder+"/deduplicated/temp_masked_oneLine.fasta"

    # Mask the duplicate regions with "$"
    cmd_mask = ["bedtools", "maskfasta", "-fi", assembly_name, "-fo", fasta_masked, "-bed", bedname, "-mc", "$"]
    print("\t"+ " ".join(cmd_mask))
    try:
        pr = subprocess.Popen(cmd_mask, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd_mask))
    except Exception as e:
        print("Error for: " + " ".join(cmd_mask))
        print("Exception:"+str(e))
        sys.exit()

    # Make the fasta as one sequence per line, needed for sed below:
    # else we have "empty/truncated" lines where we remove sequences longer than fasta line length.
    ud.make_fasta_one_line(fasta_masked, fasta_masked_oneLine)

    # Replace the "$" which marks the duplicate regions by an empty string (to remove them)
    # The backup extension is required for Mac OS ("-i'backup' -e" should work on both GNU and Mac)
    cmd_sed = "sed -i'.dupless_sed_backup' -e 's/\$//g' "+fasta_masked_oneLine
    print("\t"+cmd_sed)
    try:
        pr = subprocess.Popen(cmd_sed, shell=True, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, cmd_sed)
    except Exception as e:
        print("Error for: " + cmd_sed)
        print("Exception:"+str(e))
        sys.exit()
    
    # Cleaning temp files:
    # move the fasta_masked_oneLine to the deduplicated.fasta
    mv_cmd = ["mv", fasta_masked_oneLine, outname]
    try:
        pr = subprocess.Popen(mv_cmd, shell=False)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(mv_cmd))
    except Exception as e:
        print("Error for: "+ " ".join(mv_cmd))
        print("Exception:"+str(e))
        sys.exit()
    # remove the temp file fasta_masked
    ud.remove_file(fasta_masked)
    # remove the backup made by sed
    ud.remove_file(fasta_masked_oneLine+".dupless_sed_backup")


def generate_discarded_fasta(assembly_name, region_file, discarded_fasta):
    """Write the sequences that have been removed to a discarded.fasta file.
    The bed file should correspond to "discarded.region"
    """
    with open(discarded_fasta, "w") as discarded_handle:
        # Extract the duplicated regions
        cmd_faidx = ["samtools", "faidx", assembly_name, "-r", region_file]
        print("\t"+ " ".join(cmd_faidx))
        try:
            pr = subprocess.Popen(cmd_faidx, shell=False, stdout=discarded_handle)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd_faidx))
        except Exception as e:
            print("Error for: " + " ".join(cmd_faidx))
            print("Exception:"+str(e))
            sys.exit()

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
    
    #######discardedBed is not a Bed, it is a regions I believe
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
remove_duplications_assembly(deduplicated_assembly, assembly_name, toRemoveBed, output_folder)
generate_discarded_fasta(assembly_name, discardedBed, discarded_assembly)
phasing(toRemoveBed, blast_output,"hap2.fasta", assembly_name, output_folder)



# If we skip blast, no intermediate files created, no need to clean
if not skip_blast:
    # Cleaning the intermediate files:
    ud.remove_file(output_folder+"/All_Blasts_region_coord.tab")
    ud.remove_file(output_folder+"/assembly_HET_ONLY.fa")
    ud.remove_file(output_folder+"/assembly_HET_ONLY.fa.fai")

print("Done !\n")
print("Deduplicated assembly generated in: " + deduplicated_assembly)
print("Discarded sequences in: " + discarded_assembly)
