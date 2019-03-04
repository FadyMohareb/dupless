#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib
# non-interactive backend to avoid error : "XIO:  fatal IO error 25 (Inappropriate ioctl for device)"
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
import sys
import subprocess
from multiprocessing import Pool

import utils_dupless as ud

#==============================
#       Functions
#==============================
def get_colors(classifications):
    """
    Takes a list of classifications ("het", "hom", anythingElse) and
    Returns:
        The correspond list of colors ("red", "green", "purple").
    """
    colors = []
    for cl in classifications:
        if(cl == 'het'):
            colors.append("#d62728")
        elif(cl == 'hom'):
            colors.append("#5ac739")
        else:
            colors.append("#440a7f")
    return colors



def classify_median(value, expected_coverage):
    """
    Takes a value and an expected coverage and
    Returns:
        A classification ("het", "hom", "outlier") according to the distance to the expected_coverage.
    """
    cl = ''
    if(value != 0 and
       value <= expected_coverage/1.5):
        cl ='het'
    elif(value > expected_coverage/1.5 and
         value < expected_coverage*1.5):
        cl = 'hom'
    else:
        cl = 'outlier'
    return cl



def create_plot_coverage(starts, medians, classifications, contig_coverage, contig_name, gaps_dataframe, window_size, output_folder):
    """
    Generates the coverage plot from a list of start positions and a list of coverages (generally corresponding to a contig/scaffold).
    Saves it under output_folder/"filename".
    """
    filename = contig_name+".png"

    len_pos = len(starts)
    len_medians = len(medians)
    len_class = len(classifications)
    colors = []

    if(len_pos == len_medians and len_class == len_pos):
        colors = get_colors(classifications)

        plt.figure(figsize=(20,16), dpi=80)
        plt.plot(starts, medians, color='k', linewidth=0.5)
        plt.scatter(x=starts, y=medians, c=colors, edgecolor='none')
        # Plot horizontal line for coverage and coverage/2
        plt.plot([0, starts[-1]], [contig_coverage, contig_coverage], 'k-', linestyle='--', lw=1)
        plt.plot([0, starts[-1]], [contig_coverage/2, contig_coverage/2], 'k-', linestyle='--', lw=1)
        plt.xlabel('Position (window size='+str(window_size)+")")
        plt.ylabel('Median of read coverage for each window\n(expected coverage='+str(contig_coverage)+')')
        plt.title(str(contig_name))
        plt.ylim(0,contig_coverage*2)

        #Create custom artists
        hetArtist = plt.Line2D((0,1),(0,0), color='#d62728', marker='o')
        homArtist = plt.Line2D((0,1),(0,0), color='#5ac739', marker='o')
        outArtist = plt.Line2D((0,1),(0,0), color = '#440a7f', marker='o')
        #Create legend from custom artist/label lists
        plt.legend([hetArtist,homArtist, outArtist],['Heterozygous', 'Homozygous', 'Outlier'])

        if(gaps_dataframe.empty == False):
            gaps_df_contig = gaps_dataframe[gaps_dataframe['contig'].isin([contig_name])]
            if(len(gaps_df_contig['contig']) > 0):
                for i in range(0, len(gaps_df_contig['contig'])):
                    plt.axvspan(int(gaps_df_contig['start'].iloc[i]), int(gaps_df_contig['stop'].iloc[i]), facecolor='0.2', alpha=0.5)

        plt.savefig(output_folder+"/graphs/"+filename)
        plt.close()

    else:
        print("ERROR: there are "+len_pos+" positions for "+len_medians+" coverage values and "+len_class+" classifications !")
        with open(output_folder+"/"+contig_name+"_PLOT_ERROR.log", "w") as error_file:
                error_file.write("ERROR: there are "+len_pos+" positions for "+len_medians+" coverage values and "+len_class+" classifications !")



def get_windows_medians_contig(contig_df, expected_cov, window_size):
    """
    Generate the list of positions of the list of coverages for the contig "contig_name".
    Returns::
        A list of positions and a list of medians classifications. ('hom'=homozygous, 'het'=heterozygous, 'outlier'=possible outlier)
    """
    # Getting the last position of the contig and the median coverage on this contig (to detect het regions, we assume that homozygous is the dominant)
    last_contig_position = contig_df['position'].iloc[-1] - 1

    previous_position = 0
    current_position = 0
    window_medians = []
    starts = []
    stops = []
    classifications = []

    # For every window, we get the start, stop, median coverage and classification (hom,het or other)
    while(current_position < last_contig_position-window_size):
        previous_position = current_position
        current_position = previous_position + window_size
        starts.append(int(previous_position))
        stops.append(int(current_position))
        window = contig_df[(contig_df['position'] > previous_position) & (contig_df['position'] < current_position)]
        window_medians.append(float(window['coverage'].median()))
        classifications.append(classify_median(window_medians[-1], expected_cov))

    starts.append(int(current_position))
    stops.append(int(last_contig_position))
    window = contig_df[(contig_df['position'] > current_position) & (contig_df['position'] < last_contig_position)]
    window_medians.append(float(window['coverage'].median()))
    classifications.append(classify_median(window_medians[-1], expected_cov))

    return starts, stops, window_medians, classifications



def create_bed_het_regions(contig_name, starts, stops, classifications, output_folder):
    """
    Creates a bed file containing the regions detected as 'het' in the "classifications" list.
    The regions comes from the "positions" list. The two lists should correspond (positions[i] should refer to classifications[i])
    """
    i = 0
    start = 0

    with open(output_folder+"/individual_beds/"+contig_name+'HET.bed', 'a') as het_bed_file:
        if(len(starts) == len(classifications)):
            while(i < len(classifications)):
                if(classifications[i] == 'het'):
                    start = starts[i]
                    while(i < len(classifications)-1 and classifications[i+1] == 'het'):
                        i = i + 1
                    stop = stops[i]
                    # Name of the scaffold is "chr_start_stop", to avoid issues with faidx later, that has issues with "chr:start-stop" Ids.
                    name = contig_name+"_"+str(start)+"_"+str(stop)
                    het_bed_file.write(contig_name+"\t"+str(start)+"\t"+str(stop)+"\t"+name+"\n")
                i = i + 1
        else:
            print("Different lengths for starts and classifications:"+str(len(starts))+" and "+str(len(classifications)))
            with open(output_folder+"/"+contig_name+"_BED_ERROR.log", "w") as error_file:
                error_file.write("Different lengths for starts and classifications:"+str(len(starts))+" and "+str(len(classifications)))


# genome_expected_coverage is coming from the user "-c" or fom computation of the genome mode
# Computation of the genome mode does not work if Heterozygous peak > Homozygous peak
def coverage_histogram(contig_coverages, genome_expected_coverage, output_folder):
    """
    Takes a list of coveragesself.
    Saves the histogram of coverage under output_folder/"Histogram_coverage.png".
    """
    plt.figure(figsize=(20,16), dpi=80)
    plt.hist(contig_coverages, genome_expected_coverage, facecolor='green', alpha=0.75)
    plt.grid(True)
    plt.axvline(x=genome_expected_coverage, color='r', linestyle='--', lw=2)
    plt.xlabel("Coverage (x fold), up to expected coverage*2 only.")
    plt.ylabel('Frequency')
    plt.title("Coverage distribution, zero coverage regions ignored (Red bar is the expected coverage)")
    plt.savefig(output_folder+"/Histogram_coverage.png")
    plt.close()



def detect_genome_mode(bed_dataframe):
    """
    Detects the genome mode from the coverage values. Will be false if het peak > hom peak.
    Returns:
        The detected genome mode.
    """
    genome_mode = None
    # Remove values = 0 => if the scaffold possess mostly of gaps it can false the calculation of the mode
    bed_dataframe_filtered = bed_dataframe[bed_dataframe['coverage'] != 0]
    modes = bed_dataframe_filtered['coverage'].mode()
    if(len(modes) == 1):
            genome_mode = int(modes.iloc[0])
    else:
        print("ERROR: More than one mode found !!!")
        sys.exit(2)
    return genome_mode



# Try and catch needed here due to bug in python with multiprocessing
# ctrl-c makes the child unjoinable and pool creates new workers instead of exiting
def process_contig(contig_name):
    global BED_DF
    global GENOME_MODE
    global WINDOW_SIZE
    global OUTPUT_FOLDER
    global GAPS_DF
    global SKIP_PLOT

    # try here because this funcion is used with pool
    # the except KeyboardInterrupt will allow a gracefull exit with ctrl-c
    try:
        contig_df = BED_DF[BED_DF['contig'].isin([contig_name])]
        starts, stops, window_medians, classifications = get_windows_medians_contig(contig_df, GENOME_MODE, WINDOW_SIZE)
        if not(SKIP_PLOT):
            create_plot_coverage(starts, window_medians, classifications, GENOME_MODE, contig_name, GAPS_DF, WINDOW_SIZE, OUTPUT_FOLDER)
        create_bed_het_regions(contig_name, starts, stops, classifications, OUTPUT_FOLDER)
    except KeyboardInterrupt:
        return


BED_DF = pd.DataFrame()
GENOME_MODE = 0
WINDOW_SIZE = 0
OUTPUT_FOLDER = 0
GAPS_DF = pd.DataFrame()
#==========================================
# "Main" function that uses the other ones
#==========================================
def detect_het_regions(coverage_bed, gaps_bed, genome_mode, window_size, output_folder, nbThreads, skip_plot):
    """
    Reads the coverage bed file and produces bed and graphs of the heterozygous regions.
    Returns:
        The name of the bed files containing all the heterozygous regions.
    """
    global BED_DF
    global GENOME_MODE
    global WINDOW_SIZE
    global OUTPUT_FOLDER
    global GAPS_DF
    global SKIP_PLOT

    GENOME_MODE = genome_mode
    WINDOW_SIZE = window_size
    OUTPUT_FOLDER = output_folder
    SKIP_PLOT = skip_plot

    #GAPS_DF = pd.DataFrame()
    if(gaps_bed != None):
        GAPS_DF = pd.read_csv(gaps_bed, sep='\t', index_col=False, names=['contig', 'start', 'stop'], header=None)

    print("\nReading coverage bed, this can take a while...")
    BED_DF = pd.read_csv(coverage_bed, sep='\t', index_col=False, names=['contig', 'position', 'coverage'], header=None)
    print("Done !\n")

    # If user mode not set by user, then compute it
    # Use the mode from the whole genome to avoid wrong calculation on mostly het scaffolds.
    if(GENOME_MODE == None):
        GENOME_MODE = detect_genome_mode(BED_DF)

    if(SKIP_PLOT == False):
        # Creates the histogram of the coverage, to visually check the chosen expected coverage (genome_mode)
        bed_DF_hist = BED_DF[(BED_DF['coverage'] <= GENOME_MODE*2) & (BED_DF['coverage'] != 0)]
        coverage_histogram(bed_DF_hist['coverage'], GENOME_MODE, OUTPUT_FOLDER)

    print("The mode is :"+str(GENOME_MODE))

    print("Processing the contigs... (Creating graphs and bed files)")
    # Multiprocessing, we process each contig in parallel
    try:
        pool = Pool(nbThreads)
        pool.map(process_contig, BED_DF['contig'].unique())
        pool.close()
        pool.join()
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
        sys.exit()
    except:
        print("Error during the processing of the contigs")
        pool.terminate()
        pool.join()
        sys.exit()
        
    print("Contigs processed !\n")

    concat_bed_name = OUTPUT_FOLDER+"/Heterozygous_regions_ALL.bed"
    print("Concatenating the bed files to "+concat_bed_name+" ...")
    with open(concat_bed_name, "w") as concat_bed:
        # Combination of find and cat with "+" to avoid issue of "argument list too long"
        cmd = ["find", OUTPUT_FOLDER+"/individual_beds/", "-maxdepth", "1", "-type", "f", "-exec", "cat", "{}", "+"]
        try:
            pr = subprocess.Popen(cmd, shell=False, stdout=concat_bed)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd))
        except:
            print("Error for: " + " ".join(cmd))
            print(sys.exc_info()[0])
            sys.exit()  
    print("Bed files concatenated to: "+concat_bed_name+"\n")
    return concat_bed_name
