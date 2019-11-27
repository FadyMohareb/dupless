#!/usr/bin/python

# Dependencies:
#   python v2.7 or higher
#   samtools v1.9 or higher (important for the "-o" parameter)
#   bedtools v2.27 (lower version should now work)
#   blastn v2.6.0+
#   pandas, numpy, matplotlib, multiprocessing, getopt, biopython, sys, os, subprocess
#   sed and awk

import fileinput
import getopt
import subprocess
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd
import re
import json
import time
import numpy as np

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


def writeKmerIDs(kmerFile):#not being used
  kmerID=1
  with fileinput.FileInput(kmerFile, inplace=True, backup='.bak') as fileK:
    for line in fileK:
      print(line.replace(">", ">"+str(kmerID)+"_"), end='')
      kmerID+=0.5
      
def kmerReadAssemblyBlasting(rowDf,rowDf2,readsDf,assembly,readFastaDF,numberKmer):   
  #function that receives the blast output row of kmer vs Assembly for both contigs,      
  assemblySeqs = SeqIO.parse(open(assembly),'fasta')
  #writing sequence of the involved assembly contigs to fasta
  assemblySeqID=rowDf.tolist()[1]
  startKmerContig1=rowDf.tolist()[8]#subject start of kmers vs Contigs, contig 1
  with open("assemblyContig1.fasta", "w") as f:
    for record in assemblySeqs:
      if (record.id==assemblySeqID):
        SeqIO.write(record, f, "fasta")
         
  assemblySeqID=rowDf2.tolist()[1]
  startKmerContig2=rowDf2.tolist()[8]#subject start of kmers vs Contigs, contig 2 
  assemblySeqs = SeqIO.parse(open(assembly),'fasta')
  with open("assemblyContig2.fasta", "w") as f2:
    print("check")  
    for record in assemblySeqs:
      if (record.id==assemblySeqID):
        print("contig found")
        SeqIO.write(record, f2, "fasta")
              
        
  #writing sequences of reads to fasta  
  readsRowsDf=readsDf.loc[readsDf['qseqid'] == rowDf.tolist()[0]]  
  with open("readsCompare.fasta", "w") as f:
    for index,row in readsRowsDf.iterrows():
      readSeqID=str(int(row.tolist()[1] ))
      subsetReadDf=readFastaDF.loc[readFastaDF['read_id'] == readSeqID]
      for index,row in subsetReadDf.iterrows():
        record = SeqRecord(Seq(row.tolist()[1],IUPAC.unambiguous_dna),id=row.tolist()[0])
        SeqIO.write(record, f, "fasta")
        
  #assemble a consensus sequence of the reads
  velveth = ["velveth","output","21","-fasta","readsCompare.fasta"]
  velvetg = ["velvetg","output"]
  process = subprocess.Popen(velveth, stdout=subprocess.PIPE)
  process.wait()
  process = subprocess.Popen(velvetg, stdout=subprocess.PIPE)
  process.wait()
  #must read in the generated consensus contigs, and if there is more than one, only keep the one that has the kmer, then the longest 
  velvetSeqs = SeqIO.index("output/contigs.fa",'fasta')
  if len(list(velvetSeqs))>1:
    kmerSeqs = SeqIO.index("removalKmers.fasta",'fasta')
    with open("kmerSequence.fasta", "w") as f:
      record = SeqRecord(Seq(str(kmerSeqs[str(numberKmer)].seq),IUPAC.unambiguous_dna),id=str(numberKmer))
      SeqIO.write(record, f, "fasta")
    makeBlastDB = ["makeblastdb", "-in","output/contigs.fa","-input_type","fasta","-out","velvetContigs","-dbtype","nucl"]
    blast = ["blastn","-word_size",str(21),"-db","velvetContigs","-query","kmerSequence.fasta","-outfmt", str(6),"-out","kmerVsVelvetContigs.tab","-num_threads",str(nbThreads)]
    #kmer 32 did not hit the contigs, let's hope it is resolved by being more lenient in the alignments (reducing kmer size of the blast)
    process = subprocess.Popen(makeBlastDB, stdout=subprocess.PIPE)
    process.wait()
    process = subprocess.Popen(blast, stdout=subprocess.PIPE)
    process.wait()
    velvetContigsDF=pd.read_csv("kmerVsVelvetContigs.tab", sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']).sort_values(by=['length'])
    print("Velvet generated more than one contig")
    print(velvetContigsDF)
    with open("newconsensus.fasta", "w") as f:
        record = SeqRecord(Seq(str(velvetSeqs[str(velvetContigsDF.iloc[0]['sseqid'])].seq),IUPAC.unambiguous_dna),id=str(velvetContigsDF.iloc[0]['sseqid']))
        SeqIO.write(record, f, "fasta")
    blast1 = ["blastn","-word_size",str(27),"-db","contig2","-query","newconsensus.fasta","-outfmt", str(6),"-out",("consensusContig2-"+str(numberKmer)+".tab"),"-num_threads",str(nbThreads)]
    blast2 = ["blastn","-word_size",str(27),"-db","contig1","-query","newconsensus.fasta","-outfmt", str(6),"-out",("consensusContig1-"+str(numberKmer)+".tab"),"-num_threads",str(nbThreads)]
    if len(velvetContigsDF)==1:#make a fasta with the only (or longest) contig that matched the kmer and blast that fasta
      print(" but only one hit the kmer sequence")
    else:
      print(" and several hit the kmer sequence, so we keep the longest velvet contig")
      
  #there are cases in which there is no hit between the consensus and the contig, might have to reduce word size too  
  else: #we blast the original velvet contig file
    print("Velvet generated one contig")
    blast1 = ["blastn","-word_size",str(27),"-db","contig2","-query","output/contigs.fa","-outfmt", str(6),"-out",("consensusContig2-"+str(numberKmer)+".tab"),"-num_threads",str(nbThreads)]
    blast2 = ["blastn","-word_size",str(27),"-db","contig1","-query","output/contigs.fa","-outfmt", str(6),"-out",("consensusContig1-"+str(numberKmer)+".tab"),"-num_threads",str(nbThreads)]
    
  #blast the consensus sequence vs the contigs  
  makeBlastDB = ["makeblastdb", "-in","assemblyContig2.fasta","-input_type","fasta","-out","contig2","-dbtype","nucl"]
  process = subprocess.Popen(makeBlastDB, stdout=subprocess.PIPE)
  process.wait()
  process = subprocess.Popen(blast1, stdout=subprocess.PIPE)
  process.wait()
  makeBlastDB = ["makeblastdb", "-in","assemblyContig1.fasta","-input_type","fasta","-out","contig1","-dbtype","nucl"]
  process = subprocess.Popen(makeBlastDB, stdout=subprocess.PIPE)
  process.wait()
  process = subprocess.Popen(blast2, stdout=subprocess.PIPE)
  process.wait()
   ######################################################################must blast specifically close to the kmer position originally, as it is possible to assemble a consensus sequence that has repetitive regions: i can just blast and check if the start is between the start and end, or just extract from rowDF
  
  #read tabs into dfs, the tabs should only be of one velvet contig now, and we will keep the longest hit, only passing the contig ID and positions##################################     #must check this further, there will be several hits most surely in some cases 
  contig1DF=pd.read_csv(("consensusContig1-"+str(numberKmer)+".tab"), sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']).sort_values(by=['length'])
  contig2DF=pd.read_csv(("consensusContig1-"+str(numberKmer)+".tab"), sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']).sort_values(by=['length'])
  print(contig1DF.iloc[0,[1,8,9]].tolist())
  print(contig2DF.iloc[0,[1,8,9]].tolist())
  return (contig1DF.iloc[0,[1,8,9]].tolist(),contig2DF.iloc[0,[1,8,9]].tolist())
    
def phasing(toRemoveBed, blast_output, outname, assembly_name, output_folder):
    """Reading the heterozygous regions blast tab and the toRemoveBed to locate to what subject contigs the queries (recorded in the Bed) have matched.
    Therefore, we will generate a second bed containing subject sequences instead of query sequences, and the removal of these will be the second haplotype.
    """
    skipall=True
    if (skipall==False):
      #Secondary check as it is not carried out in the current tests previously
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
      
      #Beginning a timecount and creating empty bed files for later use
      start_time = time.time()
      toRemoveBed1 = output_folder+"/toDiscard_hap1.bed"
      process = subprocess.Popen(["touch", toRemoveBed1], stdout=subprocess.PIPE)
      process.wait()    
      toRemoveBed2 = output_folder+"/toDiscard_hap2.bed"
      process = subprocess.Popen(["touch", toRemoveBed2], stdout=subprocess.PIPE)
      process.wait()
      repetitiveBed = output_folder+"/repetitive.bed"
      process = subprocess.Popen(["touch", repetitiveBed], stdout=subprocess.PIPE)
      process.wait()
      kmerBed = output_folder+"/kmerRemoval.bed"
      process = subprocess.Popen(["touch", kmerBed], stdout=subprocess.PIPE)
      process.wait()
      
      # Needed to quickly get the lengths of the sequences
      assembly_dict = SeqIO.index(assembly_name, "fasta")
      
      #Reading the blast output and parsing everything into a pandas Dataframe
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
      

    
    #lets decide using kat sect (0) or kat comp (1) for testing purposes
    kat=1
    if (kat==0):
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
          #Writing bed file of the kat sect repetitive regions fasta  
      for index, row in repetitive_df.iterrows():
        with open(repetitiveBed, "a") as repetitiveHandle:
          repetitiveHandle.write(str(row['contig'])+"\t"+str(row['start'])+"\t"+str(row['end'])+"\n")
          
          
          
    elif (kat==1):
      skipfornow=True
      if (skipfornow==False):
        #using kat hist to create a json distanalysis file with the frequency of the homozygous peak
        popen_args=["kat","hist","sim_dup_reads100X.fa"]
        process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
        #process.wait()      
        with open("kat.hist.dist_analysis.json") as json_file:
          data = json.load(json_file)
          hom_peak=data['peaks'][int(data['hom_peak']['index']-1)]
          #print(hom_peak)####################################################or allow user input#################################################
          maxMultiplicity=int(round(float(hom_peak["mean_freq"])+float(hom_peak["stddev"])))
          minMultiplicity=int(round(float(hom_peak["mean_freq"])-float(hom_peak["stddev"])))
      
        #Using jellyfish counts and dump to get the 2x kmers of the assembly and make a dictionary with keys as sequences and values as multiplicities (all 2) 
        popen_args=["jellyfish","count","-m",str(27),"-C","-t",str(nbThreads),"-s","100M",assembly_name,"-o","jellyCountAssembly2x.jf"]
        process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
        process.wait()
        popen_args=["jellyfish","dump","-U",str(2),"-L",str(2),"-o","assembly2xKmers.fasta","jellycountAssembly2x.jf"]
        process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
        process.wait()
        assemblyKmers = SeqIO.parse(open("assembly2xKmers.fasta"),'fasta')
        assemblyKmerDict={}
        for seq in assemblyKmers:
          assemblyKmerDict[str(seq.seq)]=seq.id
        
        
  
        #Getting the homozygous peak kmers of the reads and creating another dictionary like the previous one
        popen_args=["jellyfish","count","-m",str(27),"-C","-t",str(nbThreads),"-s","100M","sim_dup_reads100X.fa","-o","jellyCountReadsAll.jf"]
        process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
        process.wait()
        popen_args=["jellyfish","dump","-o","readsKmers.fasta","-U",str(maxMultiplicity),"-L",str(minMultiplicity),"jellyCountReadsAll.jf"]
        process = subprocess.Popen(popen_args, stdout=subprocess.PIPE)
        process.wait()
        readKmers = SeqIO.parse(open("readsKmers.fasta"),'fasta')
        readKmerDict={}
        for seq in readKmers:
          readKmerDict[str(seq.seq)]=seq.id
        
        #dictionary comprehension to iterate and obtain the kmer sequences of the assembly that have been found in the homozygous peak of the reads, the keys are          the sequences
        #and the values are the multiplicities in the reads
        kmersCompared={key:readKmerDict[key] for key in assemblyKmerDict if key in readKmerDict}
        
        #dict comprehension to remove all the kmers that do not a multiplicity in the homozygous peak (double check, probably not necessary)
        removalKmers={key:kmersCompared[key] for key in kmersCompared if int(kmersCompared[key])>=minMultiplicity and int(kmersCompared[key])<=maxMultiplicity}
        print(removalKmers)      
        
        #Now to write a fasta with the previous dictionary 
        with open("removalKmers.fasta", "w") as f:
          kmerID=1
          for seq in removalKmers:
            print(str(seq))
            record = SeqRecord(Seq(seq,IUPAC.unambiguous_dna),id=str(kmerID))
            SeqIO.write(record, f, "fasta") 
            kmerID+=1
        
        
      perform=False  #blast to find the position of the kmers in the assembly and the reads(jellyfish dump does not allow this to be done previously)
      if (perform==True):         
        makeBlastDB = ["makeblastdb", "-in",assembly_name,"-input_type","fasta","-out","assemblyDB","-dbtype","nucl"]
        blast = ["blastn","-word_size",str(27),"-db","assemblyDB","-query","removalKmers.fasta","-outfmt", str(6),"-out","kmersRemovalAssembly.tab","-num_threads",str(nbThreads)]#"-max_target_seqs",str(1),"-max_hsps",str(1)
        process = subprocess.Popen(makeBlastDB, stdout=subprocess.PIPE)
        process.wait()
        process = subprocess.Popen(blast, stdout=subprocess.PIPE)
        process.wait()        
        makeBlastDB = ["makeblastdb", "-in","sim_dup_reads100X.fa","-input_type","fasta","-out","readsDB","-dbtype","nucl"]
        blast = ["blastn","-word_size",str(27),"-db","readsDB","-query","removalKmers.fasta","-outfmt", str(6),"-out","kmersRemovalReads.tab","-num_threads",str(nbThreads)]#"-max_target_seqs",str(1),"-max_hsps",str(1)
        process = subprocess.Popen(makeBlastDB, stdout=subprocess.PIPE)
        process.wait()
        process = subprocess.Popen(blast, stdout=subprocess.PIPE)
        process.wait()
        
        #also, blast reads vs assembly directly
        blast = ["blastn","-word_size",str(27),"-db","readsDB","-query",assembly_name,"-outfmt", str(6),"-out","bigAssemblyVsReads.tab","-num_threads",str(nbThreads)]#"-max_target_seqs",str(1),"-max_hsps",str(1)
        process = subprocess.Popen(blast, stdout=subprocess.PIPE)
        process.wait()
        
        
      ##########time to add the hard part, compare the reads of the kmers to the assembly of the kmers, blast them one by one.... 
      #first we make dataframes with all the results
      kmersAssemblyDF=pd.read_csv("kmersRemovalAssembly.tab", sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
      #print(kmersAssemblyDF)
      kmersReadsDF=pd.read_csv("kmersRemovalReads.tab", sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
      assemblyReadsDF=pd.read_csv("bigAssemblyVsReads.tab", sep='\t',header=None, names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
      
      #we count the time it took to get here and start counting again
      #elapsed_time  = time.time()- start_time    
      #print("Elapsed time for %s: %s" % ("steps before dataframe with all reads",str(elapsed_time)))
      start_time = time.time()
      
           
      #store all reads into a pandas dataframe for faster data retrieval    
      with open("sim_dup_reads100X.fa", "r") as fasta:  
        readIds = []
        readSequences = []
        for record in SeqIO.parse(fasta, 'fasta'): 
          readIds.append(str(record.id))
          readSequences.append(str(record.seq))   
        dictDF={'read_id':readIds,'sequence':readSequences}        
        readFastaDf=pd.DataFrame(dictDF) 
      #how long did the dataframe creation take?    
      elapsed_time  = time.time()- start_time    
      print("Elapsed time for %s: %s" % ("creating dataframe with all the reads",str(elapsed_time)))
      
      
      
      
      
          
          
          
      #iterate through each row of the dataframe kmersAssemblyDF to find the position of the kmer in the contig 1, find the reads at that position, and blast them to contig 2
      numberOfKmersIteration=len(kmersAssemblyDF.index)
      enditeration=numberOfKmersIteration-50# im doing this here just to stop the iterations quickly, the number is the amount of iterations performed
      #print(str(numberOfKmersIteration))                5 and 6 have multiple velvet contigs   
      
      
      
      
      
      #iterate every two positions, sending in the second row every time
      contigs1RowsList=[]
      contigs2RowsList=[]
      numberKmer=1
      # convert to itertuples, much faster supposedly
      for index, row in kmersAssemblyDF.iterrows():#convert to apply() for the df if possible, probably not
        print(index)
        row=kmersAssemblyDF.iloc[index]
        row2=kmersAssemblyDF.iloc[index+1]
        if (index % 2 == 0):
          #the positions of the duplications of the contigs are arbitrarily separated into lists
          contig1row,contig2row=kmerReadAssemblyBlasting(row,row2,kmersReadsDF,assembly_name,readFastaDf,numberKmer)
          contigs1RowsList.append(contig1row)
          contigs2RowsList.append(contig2row)
          numberKmer=numberKmer+1
          print("individual blasting")
          numberOfKmersIteration=numberOfKmersIteration-1
          print(str(numberOfKmersIteration))  
          if (numberOfKmersIteration==enditeration):
            break
      print(contigs1RowsList)
      with open('listfile1.txt', 'w') as f:
        for listitem in contigs1RowsList:
          f.write('%s\n' % listitem)
      with open('listfile2.txt', 'w') as f:
        for listitem in contigs2RowsList:
          f.write('%s\n' % listitem)
      #






















            #Writing bed file of the kmers to remove  
    for index, row in kmersAssemblyDF.iterrows():
      with open(kmerBed, "a") as kmerHandle:
        kmerHandle.write(str(row['contig'])+"\t"+str(row['start'])+"\t"+str(row['end'])+"\n")
      
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
  
    if (kat==0):####it would be merge, not subtract i believe
      #Subtracting repetitive features from the discard beds         
      popen_args=["bedtools","subtract","-a",toRemoveBed1,"-b",repetitiveBed]
      with open("hap1.bed","wb") as out:
        process = subprocess.Popen(popen_args, stdout=out) 
        process.wait()  
      popen_args=["bedtools","subtract","-a",toRemoveBed2,"-b",repetitiveBed]
      with open("hap2.bed","wb") as out:
        process = subprocess.Popen(popen_args, stdout=out)
        process.wait()
    elif (kat==1):
      #Merging duplicated kmer features with the discard beds  
      with open("sortedRemoveBed1.bed","wb") as out:
        popen_args=["sort","-k1,1","-k2,2n",toRemoveBed1]  
        process = subprocess.Popen(popen_args, stdout=out) 
        process.wait() 
      with open("sortedRemoveBed2.bed","wb") as out:
        popen_args=["sort","-k1,1","-k2,2n",toRemoveBed2]  
        process = subprocess.Popen(popen_args, stdout=out) 
        process.wait() 
      with open("sortedKmerBed.bed","wb") as out:
        popen_args=["sort","-k1,1","-k2,2n",kmerBed]  
        process = subprocess.Popen(popen_args, stdout=out) 
        process.wait()      
             
      popen_args=["multiIntersectBed","-i","sortedRemoveBed1.bed","sortedKmerBed.bed"]
      with open("hap1.bed","wb") as out:
        process = subprocess.Popen(popen_args, stdout=out) 
        process.wait()  
      popen_args=["multiIntersectBed","-i","sortedRemoveBed2.bed","sortedKmerBed.bed"]
      with open("hap2.bed","wb") as out:
        process = subprocess.Popen(popen_args, stdout=out)
        process.wait()      
                
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
