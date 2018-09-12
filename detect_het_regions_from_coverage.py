#!/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import subprocess

#==============================
#       Functions
#==============================
def get_color(classifications):
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
        colors = get_color(classifications)

        plt.figure(figsize=(20,16), dpi=80)
        plt.plot(starts, medians, color='k', linewidth=0.5)
        plt.scatter(x=starts, y=medians, c=colors, edgecolor='none')
        # Plot horizontal line for coverage and coverage/2
        plt.plot([0, starts[-1]], [contig_coverage, contig_coverage], 'k-', linestyle='--', lw=1)
        plt.plot([0, starts[-1]], [contig_coverage/2, contig_coverage/2], 'k-', linestyle='--', lw=1)
        plt.xlabel('Position (window size='+str(window_size)+")")
        plt.ylabel('Median read coverage')
        plt.title(str(contig_name))
        plt.ylim(0,contig_coverage*2)

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
    print("Length= "+str(last_contig_position))

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



def coverage_histogram(contig_coverages, genome_mode, output_folder):
    """
    Takes a list of coveragesself.
    Saves the histogram of coverage under output_folder/"Histogram_coverage.png".
    """
    plt.figure(figsize=(20,16), dpi=80)
    plt.hist(contig_coverages, genome_mode, facecolor='green', alpha=0.75)
    plt.grid(True)
    plt.axvline(x=genome_mode, color='r', linestyle='--', lw=2)
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



#==========================================
# "Main" function that uses the other ones
#==========================================
def detect_het_regions(coverage_bed, gaps_bed, genome_mode, window_size, output_folder):
    """
    Reads the coverage bed file and produces bed and graphs of the heterozygous regions.
    Returns:
        The name of the bed files containing all the heterozygous regions.
    """
    gaps_dataframe = pd.DataFrame()
    if(gaps_bed != None):
        gaps_dataframe = pd.read_csv(gaps_bed, sep='\t', index_col=False, names=['contig', 'start', 'stop'], header=None)

    print("Reading coverage bed, this can take a while...")
    bed_dataframe = pd.read_csv(coverage_bed, sep='\t', index_col=False, names=['contig', 'position', 'coverage'], header=None)
    print("Done !")

    # If user mode not set by user, then compute it
    # Use the mode from the whole genome to avoid wrong calculation on mostly het scaffolds.
    if(genome_mode == None):
        genome_mode = detect_genome_mode(bed_dataframe)

    # Creates the histogram of the coverage, to visually check the chosen expected coverage (genome_mode)
    bed_DF_hist = bed_dataframe[(bed_dataframe['coverage'] <= genome_mode*2) & (bed_dataframe['coverage'] != 0)]
    coverage_histogram(bed_DF_hist['coverage'], genome_mode, output_folder)

    print("The mode is :"+str(genome_mode))

    # Then process each scaffold => get regions of het and output a figure and a bed.
    for contig_name in bed_dataframe['contig'].unique():
        # We extract the contig information from the bed file. /!\ the dataframe for the contig subset still uses the row index from the original dataframe.
        print("processing contig: "+contig_name)
        contig_df = bed_dataframe[bed_dataframe['contig'].isin([contig_name])]

        starts, stops, window_medians, classifications = get_windows_medians_contig(contig_df, genome_mode, window_size)
        create_plot_coverage(starts, window_medians, classifications, genome_mode, contig_name, gaps_dataframe, window_size, output_folder)
        create_bed_het_regions(contig_name, starts, stops, classifications, output_folder)


    concat_bed_name = output_folder+"/Heterozygous_regions_ALL.bed"
    print("Concatenating the bed files to "+concat_bed_name+" ...")
    with open(concat_bed_name, "w") as concat_bed:
        process = subprocess.Popen("cat "+output_folder+"/individual_beds/*.bed", shell = True, stdout=concat_bed)
        process.wait()
    print("Done !")
    return concat_bed_name
