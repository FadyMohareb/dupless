# DupLess v1.0.0 (Duplication Less):

## A duplication removal tool for heterozygous genomes based on read coverage and pairwise alignments.

Most of the currently available assemblers are designed to work on highly inbred, homozygous species and are treating differing haplotypes as separate contigs. However, inbreeding is not always an option and attempts to assemble a highly heterozygous species often results in a heavily duplicated assembly.
For these cases, we created "DupLess", a tool capable of detecting and removing the duplicated regions issued from heterozygosity in diploid genomes.

---

## Dependencies

DupLess is supported on Mac and Linux.

You will need to have the following dependencies:

**Note:** The following python packages are already built-in from python2.7 and do not need to be installed: getopt, subprocess, multiprocessing, sys and os. Moreover, **awk** and **sed** should also be available on most systems.

- **Python v2 or v3**
- **The following python packages:** [numpy](http://www.numpy.org/ "Numpy Homepage"), [pandas](https://pandas.pydata.org/ "Pandas Homepage"), [biopython](https://biopython.org/ "biopython Homepage"), [matplotlib.pyplot](https://matplotlib.org/ "Matplotlib Homepage"), getopt, subprocess, multiprocessing, sys, os.
- [samtools v1.9](http://www.htslib.org/ "samtools Homepage") or higher (/!\ DupLess will not work with version prior to 1.9, as it needs the "-o" parameter)
- [bedtools v2.27.1](https://bedtools.readthedocs.io/en/latest/ "Bedtools Homepage") or higher (lower versions should also work now, but only v2.26 has been tested)
- [blastn v2.6.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download "Blast download page") or higher
- **awk** and **sed**


## Installation

DupLess itself is a collection of python scripts, no installation is needed. You just have to clone the repository (or directly download the python files) and run "python DupLess.py" to use it.
```
     git clone https://github.com/MCorentin/DupLess
     cd DupLess
     python DupLess.py --help
```

To install samtools 1.9:
```
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	tar -vxjf samtools-1.9.tar.bz2
	cd samtools-1.9
	make
```
To install samtools you may also have to install HTSlib (cf: https://www.biostars.org/p/328831/)

DupLess needs to have the dependencies (samtools ,bedtools, blastn, awk and sed) available in the $PATH to work.
You can run the following command to add a tool to the PATH variable:
```
	export PATH=/path/to/tool_folder/:$PATH
```

## Testing the installation:

The "test_data" folder contains two files:

 - **AT_duplicated.fa**: this is a subset of chr3 of *Arabidopsis thaliana* with artificially induced duplications (15 sequences, see below for explanations)
 - **AT_duplicated_simReads.sorted.coverage.gz**: the coverage file for these 3 contigs, based on simulated reads **you will need to unzip this file to run DupLess**.

To test if your DupLess installation works you can run the following command (~30 min):
```
     cd test_data/
     gunzip AT_duplicated_simReads.sorted.coverage.gz
     python DupLess.py -t 20 -o dupless_AT -b AT_duplicated_simReads.sorted.coverage -a AT_duplicated.fa -w 250 -c 50 -i 95 -l 500
```

The 15 duplicated sequences should be removed or heavily truncated, you can filter the fasta by length to remove remaining artifacts.

The test dataset was created with the following pipeline (to simulate duplication due to heterozygosity):

 1. Extraction of chr3 of *Arabidopsis thaliana* assembly (ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/arabidopsis_thaliana/dna/).
 2. Creation of contig assembly by splitting chr3 everytime 2 or more consecutive "N" appeared.
 3. Creation of mutated version of contig assembly with BBmap's "mutate.sh" (mutation rate: 3%).
 4. Extraction of 15 regions from the mutated genome (~33 kbp in total) and concatenation of these 15 sequences with the contig assembly (to induce duplication).
 5. Simulation of reads with 25x coverage on the whole assembly.
 6. Simulation of reads with 25x coverage only on non-duplicated regions.

After this pipeline **AT_duplicated.fa** contained 15 duplicated regions, with 25x coverage. The non-duplicated regions had ~50x coverage. DupLess was run on this dataset and managed to remove 95% of the induced duplicated sequences.

The list of duplicated contigs is available in: "test_data/additional_data/mutated_list.txt".

---


## Input Files:

**Required**

- The assembly in fasta format.

- A bed file with the coverage value for each base of the assembly. (See below for instructions on how to create this file: "Running DupLess on your own assembly:").

**Optional**

- A bed file containing the gaps coordinates in the assembly. If provided, they will be represented as grey bars on the coverage graphs.

- If you wish to skip the detection of heterozygous regions based on the coverage, you can directly input a bed file with the regions to consider for duplication. (This file is also produced during DupLess first step)


## Usage

python DupLess.py -t [nb_threads] -b [coverage.bed] -a [assembly.fasta] -w [window_size] -c [expected_coverage] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

**Required:**
     -a/--assembly               The assembly corresponding to the bed coverage in fasta format.

     -b/--bed_cov                The bed file containing the coverage at each base (can be generated with 'bedtools genomecov').
                                 /!\ If using paired end reads: make sure that you set the -w or -l option higher than the insert size,
                                     to avoid false positives due to coverage drop at the ends of contigs (because of unaligned mates).

**Optional:**
     -t/--nThreads               The number of threads (default 10) 
     -o/--out_folder             The output folder (default './DupLess_out/')

     -c/--expected_cov           The expected read coverage for the homozygous regions. The homozygosity / heterozygosity will be determined based on this value.
                                 You can determine the value to use by plotting the coverage distribution. It should correspond to the homozygous peak
                                 If no value is given, it will be based on the mode of the coverage distribution (not reliable if high heterozygosity).

     -w/--window_size            The size of the window in basepairs (default: 1000)
                                 The value of the coverage for each window will be the median of the coverage at each base.
                                 All the windows classified as 'heterozygous' will be considered for the detection of duplication.

     -g/--bed_gaps               A bed file containing the gaps along the genome. If given, the graphs will contain a grey background where the gaps are.

     -i/--blast_identity         The minimum percentage of identity between the het regions to consider them duplicates (default: 90, range 0 to 100).
     -l/--blast_length           The blast alignments with a length lower than this threshold will be filtered (default=0).

     -n/--no_plot                Skip the creation of all the plots

**Skipping part of pipeline:**
     -s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed identifying the heterozygous regions:
                                      python DupLess.py -s [het_regions_bed] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]

     -f/--filter_blast_only      Skip the detection of the heterozygous regions AND the pairwise alignments. If so, you must provide a blast ouput with -oufmt 6:
                                      python DupLess.py -f [blast_output] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]

**Other:**
     -h/--help                   Print the usage and help and exit.
     -v/--version                Print the version and exit.


## Output files

- Two fasta files containing the different versions of the deduplicated assembly. (under "output_folder/haplotypes/")
- A bed file with the identified heterozygous regions, useful if one wants to explore the regions in more details ("output_folder/Heterozygous_regions_ALL.bed").
- A histogram of the coverage distribution, to help the user decide the expected coverage value (see below).
- Graphs of the coverage along each sequences of the assembly (see below).
- The results of the blast between the heterozygous regions ("output_folder/All_Blasts_scaffolds_coord.tab").

**We recommand filtering the resulting deduplicated assembly by length as DupLess does not remove entire contigs, so some very small contigs may be present in the output**

![alt text](https://bitbucket.org/MCorentin/hetdect/src/master/figures/Histogram_coverage.png "Histogram of coverage")

![alt text](https://bitbucket.org/MCorentin/hetdect/src/master/figures/Super_scaffold_1.png "Graph of coverage along a sequence")

---


## Running DupLess on your own assembly:

### How to generate the coverage bed file ("-b/--bed_cov" option):

You need to generate a file with the coverage value at each position (format: "sequence_name   position  coverage"). You can use any pipeline you want to generate this file. We used the following pipeline to generate the coverage files durnig our testing (/!\ do not forget the "-d" option for "genomecov" as DupLess needs the coverage on every bases):
```
     bwa index genome.fa
     bwa mem -t 20 genome.fa read_1.fq read_2.fq | samtools view -b -@ 20 -o genome_reads.bam
     samtools sort -@ 20 genome_reads.bam > genome_reads.sorted.bam
     bedtools genomecov -ibam genome_reads.sorted.bam -d > genome_reads.coverage
```

You can then use "genome_reads.coverage" with the "-b/--bed_cov" parameter. 

### How to choose the right value for the expected coverage ("-c/--expected_cov" option):

The expected coverage should be the coverage corresponding to the homozygous regions. To choose it you can plot the coverage distribution from the coverage file.

You can use R:
```
cov <- read.table("illumina.coverage")
hist(cov$V3)
```
If your genome is heterozygous, you should obtain two peaks (see graph below):
![alt text](https://bitbucket.org/MCorentin/hetdect/src/master/figures/Histogram_coverage_R.png "Histogram of coverage")

The second peak corresponds to the coverage on the homozygous regions, and the value on the x-axis for the maximum of this peak corresponds to the homozygous coverage. DupLess also generates a histogram of the coverage.

### How to choose the right value for the window length ("-w/--window_size" option):

DupLess works with a window-based approach: each sequence is splitted into windows of a certain size. Each window will be catergorized as "homozygous", "heterozygous" or "outlier" depending on the coverage median for this window.

Larger window sizes will decrease the running time but may decrease the sensitivity of DupLess. Smaller window sizes will increase the sensitivity but also possibly the number false positives, this may be counterbalanced by choosing a higher value for the blast length threshold (-l/--blast_length). The optimal value depends on the fragmentation of your assembly and your objectives: removing as many duplications as possible (small window size) or reducing the number of false positives (large window size).

/!\ **When using paired end reads** (when creating the coverage bed file), a drop in coverage is expected at the extreme ends of the sequences. This is due to the fact that only one read of the pair align to these regions. We recommand setting -w or -l higher than the insert size when using paired end reads to avoid false positives.

### Trying different blast thresholds:

We have implemented two options to allow the user to skip parts of the pipeline:

 - **-s/--skip_het_detection**: to avoid the detection of the heterozygous regions. You need to give a bed file with the coordinates of the regions to blast (it can be from a previous run of DupLess: "Heterozygous_regions_ALL.bed").

```
     python DupLess.py -s [het_regions_bed] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]
```

 - **-f/--filter_blast_only"**: an option to try different blast thresholds, this will modify the sensitivity/specificity of DupLess. You need to give a blast output with "-outfmt 6" (it can be from a previous run of DupLess: "All_Blasts_scaffolds_coord.tab").

```
     python DupLess.py -f [blast_output] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]
```

---


## How are the duplication detected/removed?

DupLess workflow is composed of two main steps:

 1. The detection, based on the coverage, of heterozygous regions in the assembly.

 2. The detection of duplicates among the heterozygous regions, based on their sequences similarity (using megablast).

For the first step, DupLess processes each sequence by splitting it into windows of size defined by the "-w/--window_size" option. Then the median coverage of each window is calculated based on the coverage at each base. The window is classified in three categories depending on its median value (EC = Expected Coverage):

 - Heterozygous if:  "0 < median <=  EC / 1.5"
 - Homozygous if:    "EC < median < EC * 1.5"
 - Outlier if:       "median = 0 **OR** median >= EC * 1.5"

Only the heterozygous regions are considered for later analysis and consecutives heterozygous windows are merged together.

The second step of the analysis is the pairwise megablast alignment of the heterozygous regions. Each region is aligned against all the others, only the best hit is retained for each region. After all the hits have been found, they are filtered by identity (-i/--blast_identity) and length (-l/--blast_length). The aligned pairs that are left after the filtering are the regions to remove from the assembly. DupLess does not remove the whole region, only the part that aligned.

For each aligned pair DupLess removes one from *haplotype1.fasta* and the other from *haplotype2.fasta*, so no genomic data is lost. The one removed for *haplotype1.fasta* is always the one on the smaller sequence (contig/scaffold) of the pair, this is done to reduce the possible misassemblies introduced by DupLess, indeed aligned pairs are mostly one large sequence that align to a much smaller and almost only heterozygous sequence, see figure below. Hence, *haplotype1.fasta* is expected to be of better quality than *haplotype2.fasta*.

![alt text](small_contig.png "Small contig, will be removed from haplotype1.fasta")

Output files are produced all along DupLess pipeline, so that the user can explore in more details how the duplications have been removed:
 
 - The heterozygous regions coordinates are written to a bed file.
 - The blast results are written to a tsv file.
 - The classification of the windows for each sequence are plotted, with a color code for each category (green=homozygous, red=heterozygous, purple=outlier) 

---

## Troubleshooting:

**Q:** haplotype1.fasta is the same as the assembly despite the "toRemoveFromhap1.tsv" being not empty.

**A:** Check if your assembly does not contain Windows new line characters: "\r". You can use "cat -v file", the hidden characters will appear as "^M". To resolve this issue you can run "tr -d '\r' < assembly.fasta > assembly_corrected.fasta"


## Future work:

- Support long reads to realign to duplicated regions and detect misassemblies
- Possibly: add Mummer to check how the regions around the duplications align to each other.
- Add an option to remove whole contigs if the blast hit span > threshold% of the total length.
- Improve speed.
- Flag regions with half the coverage but no blast hits.


## Performances:

For a genome of 1Gbp, DupLess requires ~50 Gb of RAM.

For running time, the bottleneck is the pairwise blasting of heterozygous regions, an assembly with a lot of regions will increase the running time significantly.
