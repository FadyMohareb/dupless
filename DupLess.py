#!/usr/bin/python

# TO DO :
#       Do a better expected coverage calculation : detect two peaks and select one that is freq*2 from the other (can be an issue if more than two peaks)
#       If long reads: align them to duplicated regions to check for misassemblies
#       Correct coverage for edge effect on the end of the contigs (if paired ends) ?
#       Add error trapping = check Popen output number and exit if error
#       Have logs (running and error log from communicate()[1]), with command used to launch DupLess
#       Add threshold to remove whole contig if duplicated regions covers > threshold% of total length
#       Add checks for samtools, bedtools, blastn, awk and sed...
#       Add checks for het bed file format

# Dependencies:
# bedtools
# samtools 1.9 or higher (important for the "-o" parameter)
# blastn
# pandas, numpy, matplotlib, multiprocessing, getopt, biopython
# sed and awk

# For matplotlib.pyplot: needs python-tk: sudo apt-get install python-tk

import getopt
import subprocess
import sys
import os

import detect_het_regions_from_coverage as dh
import detect_duplicates_from_het_regions as dd
import utils_dupless as ud


def usage():
    """
    Prints the usage.
    """
    print("\npython DupLess.py -t [nb_threads] -w [window_size] -b [coverage.bed] -a [assembly.fasta] -c [expected_coverage] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]\n")
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
    print("     -n/--no_plot                Skip the creation of all the plots")
    print("")
    print("     -s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed with the heterozygous regions positions:")
    print("                                     python DupLess.py -t [nb_threads] -a [assembly.fasta] -s [het_regions.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]")
    print("")
    print("     -h/--help                   Print the usage and help.")


def make_haplotype(hapname, assembly_name, bedname, output_folder):
    """
    From an assembly fasta and a bed of regions to remove, creates an haplotype fasta.
    Uses "bedtools maskfasta" and "sed" to remove the regions.
    Also need "awk" to have each sequence on one line.
    """
    fasta_masked = output_folder+"/haplotypes/temp_masked.fasta"
    fasta_masked_oneLine = output_folder+"/haplotypes/temp_masked_oneLine.fasta"

    # Mask the duplicate regions with a "$"
    cmd_mask = ["bedtools", "maskfasta", "-fi", assembly_name, "-fo", fasta_masked, "-bed", bedname, "-mc", "$"]
    print("\t"+ " ".join(cmd_mask))
    try:
        pr = subprocess.Popen(cmd_mask, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd_mask))
    except:
        print("Error for: " + " ".join(cmd_mask))
        print(sys.exc_info()[0])
        sys.exit()

    # Transform to single line fasta to avoid empty lines after sed step
    # Indeed if a region is longer than the fasta wrapping (usually 80 caracters), then the fasta will contain empty lines.
    # awk from: https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file, to avoid error with biopython: "fasta-2line"
    with open(fasta_masked_oneLine, "w") as fasta_masked_oneLine_handle:
        cmd_oneLine = "awk \'/^>/{print s? s\"\\n\"$0:$0;s=\"\";next}{s=s sprintf(\"%s\",$0)}END{if(s)print s}\' "+fasta_masked
        try:
            pr = subprocess.Popen(cmd_oneLine, shell=True, stdout=fasta_masked_oneLine_handle)
            pr.communicate()
            ud.check_return_code(pr.returncode, cmd_oneLine)
        except:
            print("Error for: " + " ".join(cmd_oneLine))
            print(sys.exc_info()[0])
            sys.exit()

    # Replace the "$" which marks the duplicate regions by an empty string (to remove them)
    cmd_sed = "sed -i 's/\$//g' "+fasta_masked_oneLine
    print("\t"+cmd_sed)
    try:
        pr = subprocess.Popen(cmd_sed, shell=True, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, cmd_sed)
    except:
        print("Error for: " + " ".join(cmd_sed))
        print(sys.exc_info()[0])
        sys.exit()
    
    # remove the temp file fasta_masked
    # move the fasta_masked_oneLine to the haplotype1 or 2 fasta
    try:
        pr = subprocess.Popen(["mv", fasta_masked_oneLine, hapname], shell=False)
        pr.communicate()
        ud.check_return_code(pr.returncode, "mv "+fasta_masked_oneLine+" "+hapname)
    except:
        print("Error for: mv "+ fasta_masked_oneLine + " " + hapname)
        print(sys.exc_info()[0])
        sys.exit()
    try:
        pr = subprocess.Popen(["rm", fasta_masked], shell=False)
        pr.communicate()
        ud.check_return_code(pr.returncode, "rm "+fasta_masked)
    except:
        print("Error for: rm " + fasta_masked)
        print(sys.exc_info()[0])
        sys.exit()


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
blast_length_threshold = 0      # ...these two blast thresholds are met (min identity and min length).
het_bed = None                  # Created by the script. Bed defining the heterozygous region. Can be used to try more than one blast thresholds without rerunning everything
skip_het_dect = False           # Possibility to skip the first step (het detection) which is time consuming.
skip_plot = False               # Skip the generation of the coverage plots

try:
    opts, args = getopt.getopt(sys.argv[1:], "t:w:b:a:c:g:o:s:i:l:nh", ["nThreads=", "window_size=", "bed_cov=", "assembly=", "expected_cov=", "bed_gaps=", 
                                                                        "out_folder=", "skip_het_detection=", "blast_identity=", "blast_length=", "no_plot", "help"])
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
    elif o in ("-n", "--no_plot"):
        skip_plot = True
    elif o in ("-h", "--help"):
        usage()
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
    # And we stop if folder already exists (to avoid overwriting an already existing project)
    if(os.path.isdir(output_folder)):
        print("\nFolder '"+output_folder+"' already exists, stopping now...\n")
        sys.exit(2)

file_ok, error_mssg = ud.check_file(assembly_name)
if not file_ok:
    print("Error with option -a/--assembly: "+error_mssg)
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
# Get the path to dupless, useful to call subscripts
DupLess_folder = sys.path[0]

for folder in [output_folder, output_folder+"/individual_beds", output_folder+"/graphs", output_folder+"/individual_blasts", output_folder+"/temp", output_folder+"/haplotypes"]:
    try:
        pr = subprocess.Popen(["mkdir", folder], shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, "mkdir "+folder)
    except:
        print("Error during mkdir "+folder)
        print(sys.exc_info()[0])
        sys.exit()

# Indexing the fasta, needed later on for extraction of het regions
# Also a good way to check if samtools exists at the start of the script
ud.index_fasta_file(assembly_name)


if not skip_het_dect:
    # Launch the bed and graph creation for heterozygous regions, detection based on coverage values.
    het_bed = dh.detect_het_regions(coverage_bed, gaps_bed, expected_coverage, window_size, output_folder, nbThreads, skip_plot)

# Check if the het_bed file exists (wether it has been created by the step before or just given by the user)
file_ok, error_mssg = ud.check_file(het_bed)
if file_ok:
    # Launch pairwise blast comparisons between the detected heterozygous regions to detect duplication
    dd.detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads, DupLess_folder)

    # Filter the blasts by identity and length.
    print("Filtering blast results with "+str(blast_identity_threshold)+"% identity and min length of "+str(blast_length_threshold)+" bp :")
    cmd_filter = ["python", DupLess_folder+"/filter_blast_results.py", output_folder+"/All_Blasts_scaffolds_coord.tab", str(blast_identity_threshold), str(blast_length_threshold), assembly_name, output_folder]
    print(" ".join(cmd_filter))
    try:
        pr = subprocess.Popen(cmd_filter, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd_filter))
    except:
        print("Error for: " + " ".join(cmd_filter))
        print(sys.exc_info()[0])
        sys.exit()  
    print("Blast filtered !\n")

    # Create the haplotype from the bed files resulting from blast filtration.
    print("Generating the haplotype fasta files from the blast results...")
    make_haplotype(output_folder+"/haplotypes/haplotype1.fasta", assembly_name, output_folder+"/toRemoveFromhap1.bed", output_folder)
    print("Haplotype 1 generated in :" + output_folder+"/haplotypes/haplotype1.fasta")
    make_haplotype(output_folder+"/haplotypes/haplotype2.fasta", assembly_name, output_folder+"/toRemoveFromhap2.bed", output_folder)
    print("Haplotype 2 generated in :" + output_folder+"/haplotypes/haplotype2.fasta")
else:
    print("Error with the heterozygous bed: "+error_mssg)
    usage()
    sys.exit(2)

print("Done !\n")