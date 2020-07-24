#!/usr/bin/python

# Part of DupLess dealing with the detection of duplications from the heterozygous
# regions detected by the previous step (detect_het_regions_from_coverage)

import subprocess
from multiprocessing import Pool
from Bio import SeqIO
import os
import sys

import utils_dupless as ud


def get_lengths_fasta(fasta):
    """Take a fasta file name.
    Returns:
        A list containing the length(s) of the sequence(s)
    """
    lengths = []
    fasta_handle = SeqIO.parse(fasta, "fasta")

    for record in fasta_handle:
        lengths.append(len(record.seq))

    return lengths


# Need to use global name for multiprocessing
# We can use map only with a function and a list of values (here the reigon ID)
def filter_fasta_file(ID):
    """Remove the sequence corresponding to ID from HET_FASTA_NAME and writes the result to "filtered_fasta".
    Needed to make a blast database (to avoid self hits)
    """
    global HET_FASTA_NAME
    global OUTPUT_FOLDER

    # The fasta with the sequence correspond to "ID" removed, it will be used to make the blast database
    filtered_fasta = OUTPUT_FOLDER+"/temp/Het_regions_filtered_"+ID+".fasta"

    file_ok, error_mssg = ud.check_file(HET_FASTA_NAME)
    if not(file_ok):
        print(error_mssg)
        sys.exit(2)

    old_fasta = SeqIO.parse(HET_FASTA_NAME, "fasta")
    ID = ID.rstrip("\n")
    ID_found = False

    # In a loop, we write all the sequences except the one with the ID
    with open(filtered_fasta, "w") as new_fasta:
        for seq_record in old_fasta:
            # If match, do nothing, else write the sequence
            if(str(seq_record.name) != ID):
                SeqIO.write(seq_record, new_fasta, "fasta")
            else:
                ID_found = True

    if not ID_found:
        print("Warning: ID '"+str(ID)+"' not found.")
    
    return ID_found



def extract_and_blast(region, het_fasta, output_folder):
    """Create commands for blasting a region (format "name:start-stop") against all the other het regions (needs a fasta file with all the het regions: "het_fasta").
    This will be use with multiprocessing
    Returns the commands (strings) to:
        Extract the het region from the het fasta.
        Create the blast DB from the filtered het fasta.
        Blast the het region against the blast DB.
    """
    region = region.rstrip("\n")
    region_temp_fasta = output_folder+"/temp/"+region+"_temp.fasta"
    region_temp_blastdb = output_folder+"/temp/"+region+"_temp.blastdb"
    # We remove the sequence from the list of het sequences to avoid self hit with blast.
    filtered_fasta = output_folder+"/temp/Het_regions_filtered_"+region+".fasta"
    
    # List of commands that will be run in parallel
    extract = "samtools faidx "+het_fasta+" "+region+" -o "+region_temp_fasta
    makeBlastDB = "makeblastdb -in "+filtered_fasta+" -input_type fasta -out "+region_temp_blastdb+" -dbtype nucl"
    megablast = "blastn -task megablast -db "+region_temp_blastdb+" -query "+region_temp_fasta+" -outfmt 6 -out "+output_folder+"/individual_blasts/"+region+".tab -max_target_seqs 1 -max_hsps 1"

    return extract, makeBlastDB, megablast



def run_cmd(cmd):
    """Used by the multiprocessing.Pool to run commands in parallel
    """
    # Here the try and expect are used to resolve a bug in python, when using multiprocessing and trying to ctrl-c:
    # Subprocesses do not return anything and the pool just creates the next worker, so the ctrl-c does not stop everything
    # We need to catch the keyboard exception here AND in the main thread to exit gracefully
    try:
        subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    except KeyboardInterrupt:
        return


def get_assembly_coordinates_from_blast_results(region_name, blast_start, blast_stop):
    """Translate the blast coordinates (which are relative to the het. regions) to scaffold coordinates.
    The blast hits are relative to the het region coordinates: "region_name start stop" (region_name is in the format: "scaffold_start_stop").
    We do this because we want to remove sequences according to the scaffolds coordinates.
    Used by: "convert_region_coord_to_scaffold_coord"
    Returns:
        The scaffold name in the assembly ("scaffold_start_stop" becomes "scaffold").
        The start position of the blast hit in the assembly "start het region + start blast hit"
        The stop position of the blast hit in the assembly "start het region + blast hit length"
        Basically, we just "shift" the blast coordinates by "region start" to get the scaffold coordinates !
    """
    # We split the region by underscores (name is under the format: "scaffold_regionStart_regionStop")
    name_split = region_name.split("_")
    # If the original scaffold name contained underscores, we put them back with join 
    # The last 2 elements of the split are the start and end positions of the region
    scaffold_name = "_".join(name_split[:-2])
    region_start = int(name_split[-2])
    #region_end = int(name_split[-1])

    # The part to remove from the scaffold is from: "region start + blast start" to "region start + blast length"
    scaffold_blast_start = region_start + int(blast_start)
    scaffold_blast_end = region_start + int(blast_stop)

    return scaffold_name, scaffold_blast_start, scaffold_blast_end



def convert_region_coord_to_scaffold_coord(region_blasts_name, output_folder):
    """Blast hits coordinates are relatives to the heterozygous regions
    But we need the scaffolds coordinates to remove from reference fasta
    This function converts the region coordinates to scaffold coordinates
    """
    scaffold_coord_blasts_name = output_folder+"/All_Blasts_scaffolds_coord.tab"

    print("Converting Blast results in region coordinates to scaffold coordinates...")
    with open(region_blasts_name, "r") as blasts_regions:
        with open(scaffold_coord_blasts_name, "w") as blasts_scaffolds:
            for blast in blasts_regions:
                tabs = blast.split("\t")
                    
                # The QUERY should always be "start stop", but if there is a inverted alignment, the SUBJECT can be "stop   start".
                # So we get the start by taking the min and the stop with the max
                query_start = int(tabs[6])
                query_end = int(tabs[7])

                subject_start = min(int(tabs[8]), int(tabs[9]))
                subject_end = max(int(tabs[8]), int(tabs[9]))

                # We need to get the original positions, since "query_start" and "query_end" correspond to the het regions positions and not the scaffold position.
                q_scaffold_name, q_scaffold_start, q_scaffold_end = get_assembly_coordinates_from_blast_results(tabs[0],query_start, query_end)
                s_scaffold_name, s_scaffold_start, s_scaffold_end = get_assembly_coordinates_from_blast_results(tabs[1], subject_start, subject_end)

                # We want to keep the order of start and end for the subject if there is an reverse alignment.
                if(int(tabs[8]) < int(tabs[9])):
                    blasts_scaffolds.write(q_scaffold_name+"\t"+s_scaffold_name+"\t"+tabs[2]+"\t"+tabs[3]+"\t"+tabs[4]+"\t"+tabs[5]+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\t"+str(s_scaffold_start)+"\t"+str(s_scaffold_end)+"\t"+tabs[10]+"\t"+tabs[11])
                else:
                    blasts_scaffolds.write(q_scaffold_name+"\t"+s_scaffold_name+"\t"+tabs[2]+"\t"+tabs[3]+"\t"+tabs[4]+"\t"+tabs[5]+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\t"+str(s_scaffold_end)+"\t"+str(s_scaffold_start)+"\t"+tabs[10]+"\t"+tabs[11])
    return scaffold_coord_blasts_name



def extract_heterozygous_regions(assembly, heterozygous_bed, output_folder):
    """Create a fasta file containing all the heterozygous regions from the "heterozygous_bed" file
    The bed file comes from the previous step of DupLess (or given by the user). Used for pairwise blast.
    Returns:
        The name of the fasta file with the heterozygous regions
    """
    global HET_FASTA_NAME

    cmd = ["bedtools", "getfasta", "-fi", assembly, "-bed", heterozygous_bed, "-name", "-fo", HET_FASTA_NAME]
    try:
        pr = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd))
    except Exception as e:
        print("Error for: " + " ".join(cmd))
        print("Exception:"+str(e))
        sys.exit()


def concatenate_blast_results(output_folder):
    """Concatenate all the individual blasts into one file
    Return:
        The name of the file with concatenated blast results
    """
    region_blasts_name = output_folder+"/All_Blasts_region_coord.tab"
  
    print("Concatenating the blast results to "+region_blasts_name+" ...")
    with open(region_blasts_name, "w") as blasts_regions:
        # Combination of find and cat with "+" to avoid issue of "argument list too long"
        cmd = ["find", output_folder+"/individual_blasts/", "-maxdepth", "1", "-type", "f", "-exec", "cat", "{}", "+"]
        try:
            pr = subprocess.Popen(cmd, shell=False, stdout=blasts_regions)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd))
        except Exception as e:
            print("Error for: " + " ".join(cmd))
            print("Exception:"+str(e))
            sys.exit() 

    return region_blasts_name


def do_blast_parallel(threads, output_folder):
    """Do the blast in parallel with a Pool.
    The blast results are written to "output_folder/individual_blasts/"
    """
    global HET_FASTA_NAME

    filter_cmds = list()
    extract_cmds = list()       # Extract the region
    makeBlastDB_cmds = list()   # Create a copy of the heterozygous regions fasta without the region and create a database
    megablast_cmds = list()     # Blast the region against the database
    n_process = 0

    pl = Pool(threads)

    # We do each step by batch to parallelise, all the extractions, then all the makeBlastDB, then all the megablasts
    het_fasta = SeqIO.parse(HET_FASTA_NAME, "fasta")
    for seq_record in het_fasta:
        # Create the list of commands to process 1 region:
        # 1) filter region from the fasta (to make blastdb and avoid self hits)
        # 2) extract the region from the fasta
        # 3) make the blast database form filtered fasta
        # 4) command to launch megablast between region and blast database
        extract, makeBlastDB, megablast = extract_and_blast(seq_record.name, HET_FASTA_NAME, output_folder)
        filter_cmds.append(seq_record.name)
        extract_cmds.append(extract)
        makeBlastDB_cmds.append(makeBlastDB)
        megablast_cmds.append(megablast)
        n_process = n_process + 1

        # If we get more commands than nbThreads specified, we launch them in parallel
        if n_process >= int(threads):
            # Bug in python: if ctrl-c during multiprocessing, the children becomes unjoinable
            # To resolve this, have to add a keyboardInterrupt exception to both the main thread AND the subprocessess
            try:
                id_founds = pl.map(filter_fasta_file, filter_cmds)
                pl.map(run_cmd, extract_cmds)
                pl.map(run_cmd, makeBlastDB_cmds)
                pl.map(run_cmd, megablast_cmds)
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pl.terminate()
                pl.join()
                sys.exit()
            except Exception as e:
                print("Error during the processing of the contigs:")
                print("Exception: "+str(e))
                pl.terminate()
                pl.join()
                sys.exit()

            n_process = 0
            filter_cmds = list()
            extract_cmds = list()
            makeBlastDB_cmds = list()
            megablast_cmds = list()

            # We empty the temp folder as we go along, this avoid creating a temp folder with a lot of files
            ud.empty_folder(output_folder+"/temp/")

    # Do the remaining commands (the ones remaining because "n_process <= int(nbThreads)")
    try:
        id_founds = pl.map(filter_fasta_file, filter_cmds)
        pl.map(run_cmd, extract_cmds)
        pl.map(run_cmd, makeBlastDB_cmds)
        pl.map(run_cmd, megablast_cmds)
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pl.terminate()
        pl.join()
        sys.exit()
    except Exception as e:
        print("Error during the processing of the contigs:")
        print("Exception: "+str(e))
        pl.terminate()
        pl.join()
        sys.exit()
    
    ud.empty_folder(output_folder+"/temp/")
    pl.close()
    pl.join()



# Used from "DupLess.py", used after "detect_dupl_regions()" 
def filter_blast_results(blast_filename, blast_identity_threshold, blast_length_threshold, assembly, output_folder):
    """Filter the blast results on identity and length.
    Write the valid blast hits to "duplications_to_remove.bed" and "discarded.regions"
    These files will contain the regions to remove from the assembly:
        - "duplications_to_remove.bed": in bed format
        - "discarded.regions": in samtools regions format (to use it later with "samtools faidx -r discarded.regions")
    """

    file_ok, error_mssg =  ud.check_file(blast_filename)
    if not(file_ok):
        print("Blast output file could not be found. See error message below for more details:")
        print(error_mssg)
        sys.exit(2)

    file_ok, error_mssg = ud.check_file(assembly)
    if not(file_ok):
        print("Assembly fasta could not be found. See error message below for more details:")
        print(error_mssg)
        sys.exit(2)

    # Needed to quickly get the lengths of the sequences
    print("Creating dictionary for the assembly...")
    assembly_dict = SeqIO.index(assembly, "fasta")

    print("\nFiltering blast results...")
    toRemove_name = output_folder+"/toDiscard.bed"
    process = subprocess.Popen(["touch", toRemove_name], stdout=subprocess.PIPE)
    process.wait()

    discarded_name = output_folder+"/toDiscard.regions"
    process = subprocess.Popen(["touch", discarded_name], stdout=subprocess.PIPE)
    process.wait()

    to_keep = list()
    with open(blast_filename, "r") as all_blasts:
        # For each duplicated pair, we write the coordinates to the ".bed" and ".regions" files
        with open(toRemove_name, "a") as toRemoveHandle, open(discarded_name, "a") as discar_handle:
            for blast in all_blasts:
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

                    q_length = len(assembly_dict[query_name])
                    s_length = len(assembly_dict[subject_name])

                    # If "S1 S2" is a valid blast hit, we don't want to eliminate "S2 S1" later (double removal),
                    # So we add the reversed pair "S2_start_stop;S1_start_stop" in to_keep
                    # We add the start and stop in case there are several hits between different regions on the same two scaffolds.
                    region1 = query_name+"_"+str(query_start)+"_"+str(query_end)
                    region2 = subject_name+"_"+str(subject_start)+"_"+str(subject_end)
                    pair = region1+";"+region2
                    reversed_pair = region2+";"+region1

                    if(pair not in to_keep):
                        # Instead of randomly removing heterozygous regions, we remove regions from the smallest contig, this avoid introducing too many misassemblies,
                        # indeed, sometimes the duplication is a small contig that maps to part of a big one: it is better to remove a small contig that removing the middle of a large one.
                        # This will reduce the chance of cutting a gene in half and hence downgrading the resulting assembly quality.
                        if(q_length < s_length):
                            toRemoveHandle.write(query_name+"\t"+str(query_start)+"\t"+str(query_end)+"\n")    
                            discar_handle.write(query_name+":"+str(query_start)+"-"+str(query_end)+"\n")
                        else:
                            toRemoveHandle.write(subject_name+"\t"+str(subject_start)+"\t"+str(subject_end)+"\n")
                            discar_handle.write(subject_name+":"+str(subject_start)+"-"+str(subject_end)+"\n")

                        to_keep.append(reversed_pair)

    return(toRemove_name, discarded_name)


def detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads):
    """Main script to launch the comparison between the het regions.
    Writes the results of Blast in: "output_folder/individual_blasts/"
    Return:
        The name of the file with the concatenated blast outputs.
    """
    global HET_FASTA_NAME
    HET_FASTA_NAME = output_folder+"/assembly_HET_ONLY.fa"

    global OUTPUT_FOLDER
    OUTPUT_FOLDER = output_folder

    # Extract all the heterozygous regions from the assembly and create a fasta file with it.
    extract_heterozygous_regions(assembly_name, het_bed, output_folder)

    # bedtools before v2.27(?) adds the region to the fasta header during getfasta
    # Checks if bedtools < 2.27 and remove 
    if(ud.check_old_bedtools_version()):
        # The backup extension is required for Mac OS ("-i 'backup' -e" should make it work on both GNU and Mac)
        cmd = "sed -i'.dupless_sed_backup' -e 's/::.*//g' " + HET_FASTA_NAME
        print("bedtools is older than v2.27, removing the part behind the '::'\n\t"+cmd)
        try:
            pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            pr.communicate()
            ud.check_return_code(pr.returncode, cmd)
            ud.remove_file(output_folder+"/assembly_HET_ONLY.fa.dupless_sed_backup")
        except Exception as e:
            print("Error for: " + cmd)
            print("Exception: "+str(e))
            sys.exit()

    # Index fasta of heterozygous regions, needed to extract sequences for blast
    ud.index_fasta_file(HET_FASTA_NAME)

    # Do the blast in parallel, each blast writes its results to "output_folder/individual_blasts/region_start_stop.tab"
    print("Blasting each heterozygous regions against the others, this could take a while....")
    do_blast_parallel(nbThreads, output_folder)
    print("Blast done !\n")
    
    # Concatenate the blast results and convert blast coordinates from het regions to coordinates on scaffolds
    region_blasts_name = concatenate_blast_results(output_folder)
    print("Blast files concatenated to: "+region_blasts_name)

    # We need to output the scaffold coordinates for the next step
    # The positions from the blast outputs are relative to the het regions
    blast_output = convert_region_coord_to_scaffold_coord(region_blasts_name, output_folder)
    print("Blast files with scaffold coordinates written to: "+blast_output+"\n")

    return blast_output
    