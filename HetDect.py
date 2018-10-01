#!/bin/python

# TO DO :
#       Add legend for coverage and coverage/2 lines on plot
#       Do a better expected coverage calculation : detect two (add option for ploidy ?) peaks and select one that is freq*2 from the other
#       Remove previous folders: /!\ individual beds are APPENDED if ran again

# /!\ Values next to gaps are "dragged down" by the 0 coverage of the gaps and could be misidentified as heterozygotes.

# Realign reads to duplicated regions to try to find misassemblies ?

# Dependencies:
# bedtools
# samtools 1.9 or higher
# blastn
# pandas, numpy, matplotlib, multiprocessing

import getopt
import subprocess
import sys
import os

import Bio.SeqIO as SeqIO
from Bio.SeqIO import FastaIO

import detect_het_regions_from_coverage as dh
import detect_duplicates_from_het_regions as dd

def usage():
    """
    Prints the usage.
    """
    print("\npython HetDect.py -t [nb_threads] -w [window_size] -b [coverage.bed] -a [assembly.fasta] -c [expected_coverage] -g [gaps.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]\n")
    print("\nOptions:\n")
    print("     -t/--nThreads               The number of threads (default 20).")
    print("     -o/--out_folder             The output folder (default is the current directory).")
    print("")
    print("     -b/--bed_cov                The bed file containing the coverage for each position (can be generated with bedtools genomecov).")
    print("     -a/--assembly               The assembly corresponding to the bed coverage in fasta format.")
    print("     -g/--bed_gaps               A bed file contaning the gaps along the genome. If given, the graphs will contain a grey background where the gaps are.")
    print("     -w/--window_size            The size of the window. The value of the read coverage will be the median of the values inside each window (default: 1000).")
    print("     -c/--expected_cov           The expected read coverage along the genome. The homozygosity / heterozygosity will be determined based on this value.")
    print("                                 If no value is given, it will be based on the mode of the coverage distribution (not reliable if high heterozygosity).")
    print("     -i/--blast_identity         The minimum percentage of identity between the het region and the blast hit to consider it valid (default: 90, range 0 to 100).")
    print("     -l/--blast_length           The minimum length for the blast hit to be considered as valid (default=0).")
    print("")
    print("     -s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed with the heterozygous regions positions:")
    print("                                     python HetDect.py -t [nb_threads] -a [assembly.fasta] -s [het_regions.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]")
    print("")
    #print("     -p/--plot_histogram         Only plot the coverage histogram (useful to determine the homozygous peak  and so the expected_coverage)")
    #print("                                     python HetDect.py -p -b [coverage.bed]")
    #print("")
    print("     -h/--help                   Print the usage and help.")


def check_file_with_option(filename, option):
    """
    Checks if "filename" exists and is a file.
    Option is used to display an informative error message.
    Returns:
        True if file exists and is a file.
        False if filename==None or is not a file (Also prints an error message).
    """
    file_ok = True
    if(filename == None):
        print("\nOption "+option+" is missing.\n")
        file_ok = False
    else:
        if not os.path.isfile(filename):
            print("\nValue for option "+option+" is not a file.\n")
            file_ok = False
    return file_ok


def make_haplotype(hapname, assembly_name, bedname, output_folder):
    """
    From an assembly fasta and a bed of regions to remove, creates an haplotype fasta.
    Uses "bedtools maskfasta" and "sed" to remove the regions.
    """
    fasta_masked = output_folder+"haplotypes/temp_masked.fasta"
    fasta_masked_oneLine = output_folder+"haplotypes/temp_masked_oneLine.fasta"

    cmd_mask = ["bedtools", "maskfasta", "-fi", assembly_name, "-fo", fasta_masked, "-bed", bedname, "-mc", "$"]
    print(" ".join(cmd_mask))
    process = subprocess.Popen(cmd_mask, stdout=subprocess.PIPE)
    process.wait()

    # Transform to single line fasta to avoid empty lines after sed step (fasta-2line format)
    # Indeed if a region is longer than the fasta wrapping (usually 80 caracters), then the fasta will contain empty lines.
    fasta_masked_handle = SeqIO.parse(fasta_masked, "fasta")
    with open(fasta_masked_oneLine, "w") as fasta_masked_oneLine_handle:
        for seq_record in fasta_masked_handle:
            SeqIO.write(seq_record, fasta_masked_oneLine_handle, "fasta-2line")
  
    cmd_sed = "sed -i 's/\$//g' "+fasta_masked_oneLine
    print(cmd_sed)
    process = subprocess.Popen(cmd_sed, shell=True, stdout=subprocess.PIPE)
    process.wait()

    process = subprocess.Popen(["mv", fasta_masked_oneLine, hapname])
    process.wait()

    process = subprocess.Popen(["rm", fasta_masked])
    process.wait()

    print("Haplotype generated in "+hapname)

#=================================================================
#                           GetOpt                               =
#=================================================================
window_size =  1000         # The coverage of each window will be based on the median of the coverages inside the window.
coverage_bed = None         # Bed file with the coverage value for each position. Can be produced with "bedtools coverage".
assembly_name = None        # Assembly in fasta format, used to extract the het regions and also check the scaffold lengths.
expected_coverage = None    # Any window with coverage < expected_cov/1.5 will be considered as heterozygous.
gaps_bed = None             # Optional. Draw gaps as grey bars on the graphs.
output_folder = "./"
nbThreads = 20
blast_identity_threshold = 90   # Two regions will be considered duplicated if...
blast_length_threshold = 0      # these two blast thresholds are met.
het_bed = None                  # Created by the script. Bed defining the heterozygous region.
skip_het_dect = False           # Possibility to skip the first step (het detection) which is time consuming.

try:
    opts, args = getopt.getopt(sys.argv[1:], "t:w:b:a:c:g:o:s:i:l:h", ["nThreads=", "window_size=", "bed_cov=", "assembly=", "expected_cov=", "bed_gaps=", "out_folder=", "skip_het_detection=", "blast_identity=", "blast_length=", "help"])
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
    elif o in ("-i", "--blast_identity"):
        blast_identity_threshold = int(a)
    elif o in ("-l", "--blast_length"):
        blast_length_threshold = float(a)
    elif o in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        assert False, "Unhandled option !"


# If we do not skip the het step, then we need the coverage bed
if not skip_het_dect:
    if not check_file_with_option(coverage_bed, "-b/--bed_cov"):
        usage()
        sys.exit(2)

if not check_file_with_option(assembly_name, "-a/--assembly"):
    usage()
    sys.exit(2)

if(window_size <= 0):
    print("The window size can not be lower than 0.\n")
    usage()
    sys.exit(2)

if(nbThreads <= 0):
    print("The number of threads can not be lower than 0.\n")
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
HetDect_folder = sys.path[0]

for folder in [output_folder, output_folder+"/individual_beds", output_folder+"/graphs", output_folder+"/individual_blasts", output_folder+"/temp", output_folder+"/haplotypes"]:
    process = subprocess.Popen(["mkdir", folder], stdout=subprocess.PIPE)
    process.wait()


if not skip_het_dect:
    # Launch the bed and graph creation for heterozygous regions, detection based on coverage values.
    het_bed = dh.detect_het_regions(coverage_bed, gaps_bed, expected_coverage, window_size, output_folder, nbThreads)

if check_file_with_option(het_bed, "-s/--skip_het_dect"):
    # Launch pairwise blast comparison between the detected heterozygous regions to remove duplication
    dd.detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads, HetDect_folder)

    # Filter the blasts by identity and length.
    print("Filtering blast results with "+str(blast_identity_threshold)+"% identity and min length of "+str(blast_length_threshold)+" bp :")
    cmd_filter = ["python", HetDect_folder+"/filter_blast_results.py", output_folder+"/All_Blasts_scaffolds_coord.tab", str(blast_identity_threshold), str(blast_length_threshold), assembly_name, output_folder]
    print(" ".join(cmd_filter))
    process = subprocess.Popen(cmd_filter, stdout=subprocess.PIPE)
    process.wait()
    print("Blast filtered !\n")

    # Create the haplotype from the bed files resulting from blast filtration.
    print("Generating the haplotype fasta files from the blast results...")
    make_haplotype(output_folder+"/haplotypes/haplotype1.fasta", assembly_name, output_folder+"/toRemoveFromhap1.bed", output_folder)
    make_haplotype(output_folder+"/haplotypes/haplotype2.fasta", assembly_name, output_folder+"/toRemoveFromhap2.bed", output_folder)
    print("Haplotypes created !\n")
else:
    usage()
    sys.exit(2)

print("Done !\n")
