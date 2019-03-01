# DupLess v1.0.0 (Duplication Less):

## Deduplication for assemblies of heterozygous genomes based on read coverage and blast pairwise alignments.

Most of the currently available assemblers are designed to work on highly inbred, homozygous species and are treating differing haplotypes as separate contigs. However, inbreeding is not always an option and attempts to assemble a highly heterozygous species often results in a heavily duplicated assembly.
For these cases, we created "DupLess", a tool capable of detecting and removing the duplicated regions issued from heterozygosity in diploid genomes.

DupLess workflow is composed of two main steps:

 1. The detection, based on the coverage, of heterozygous regions in the assembly.

 2. The detection of duplicates among the heterozygous regions, based on their sequences similarity (using blast).

---

## Dependencies

DupLess is supported on Mac and Linux.

You will need to have the following dependancies:

**Note:** The following python packages are already built-in from python2.7 and do not need to be installed: getopt, subprocess, multiprocessing, sys and os. Moreover, **awk** and **sed** should also be available on most systems.

- **Python v2 or v3**
- **The following python packages:** [numpy](http://www.numpy.org/ "Numpy Homepage"), [pandas](https://pandas.pydata.org/ "Pandas Homepage"), [biopython](https://biopython.org/ "biopython Homepage"), [matplotlib.pyplot](https://matplotlib.org/ "Matplotlib Homepage"), getopt, subprocess, multiprocessing, sys, os.
- [samtools v1.9](http://www.htslib.org/ "samtools Homepage") or higher (/!\ DupLess will not work with version prior to 1.9, as it needs the "-o" parameter)
- [bedtools v2.27.1](https://bedtools.readthedocs.io/en/latest/ "Bedtools Homepage") or higher (lower versions should also work now, but only v2.26 has been tested)
- [blastn v2.6.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download "Blast download page") or higher
- **awk** and **sed**


## Installation

DupLess is a collection of python scripts, no installation is needed. You just have to clone the repository (or directly download the python files) and run "python DupLess.py" to use it.
```
     git clone https://github.com/MCorentin/DupLess
     cd DupLess
     python DupLess.py --help
```

To install samtools 1.9 (version 1.9 is not yet available from "apt-get install"):
```
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	tar -vxjf samtools-1.9.tar.bz2
	cd samtools-1.9
	make
```
To install samtools you may also have to install HTSlib (cf: https://www.biostars.org/p/328831/)

DupLess needs to have the tools in the PATH to work.
You can run the following command to add a tool to the PATH variable.
```
	export PATH=/path/to/tool_folder/:$PATH
```

## Testing the installation:

The "test_data" folder contains to files:
 - genome.fasta: a subset of S. chilense, containing only 3 contigs.
 - illumina.coverage: the coverage file for these 3 contigs.

To test if your DupLess installation works you can run the following command:
```
     /usr/bin/python2.7 /home/corentin/git_scripts/hetdect/DupLess.py -t 5 -o test_dupless -b illumina.coverage -a genome.fasta -w 500 -c 80 -i 80 -l 100
```

You can compare your output to the one under "exemple_output".

---


## Input Files:

**Required**

- The assembly in fasta format.

- A bed file with the coverage value for each base of the assembly. (See below for instructions on how to create this file).

**Optional**

- A bed file containing the gaps coordinates in the assembly. If provided, they will be represented as grey bars on the coverage graphs.

- If you wish to skip the detection of heterozygous regions based on the coverage, you can directly input a bed file with the regions to consider for duplication. (This file is also produced during DupLess first step)


---

## Usage

	python DupLess.py -t [nb_threads] -w [window_size] -b [coverage.bed] -a [assembly.fasta] -c [expected_coverage] -g [gaps.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

**Options:**

     -t/--nThreads               The number of threads to use (default 10)
     -o/--out_folder             The output folder (default is './DupLess_out/').

     -b/--bed_cov                REQUIRED: The bed file containing the coverage for each position (can be generated with bedtools genomecov).
     -a/--assembly               REQUIRED: The assembly corresponding to the bed coverage in fasta format.
     -g/--bed_gaps               A bed file contaning the gaps along the genome. If given, the graphs will contain a grey background where the gaps are.
     
     -w/--window_size            The size of the window. The value of the read coverage will be the median of the values inside each window (default: 1000).
     -c/--expected_cov           The expected read coverage along the genome. The homozygosity / heterozygosity will be determined based on this value. You can assess this value by plot the coverage distribution.
                                 If no value is given, it will be based on the mode of the coverage distribution (not reliable if high heterozygosity).
     -i/--blast_identity         The minimum percentage of identity between the het region and the blast hit to consider it valid (default: 90, range 0 to 100).
     -l/--blast_length           The minimum length for the blast hit to be considered as valid (default=0).


**Other:**

	-n/--no_plot		        Skip the creation of all plots.

	-s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed with the heterozygous regions positions:
                                     python DupLess.py -s [het_regions.bed] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

	-h/--help                   Print the usage and help and exit.
     -v/--version                Print the version and exit.


## Output

- **Two fasta files containing the different versions of the deduplicated assembly.** (under "output_folder/haplotypes/")
- A bed file with the identified heterozygous regions, useful if one wants to explore the regions in more details ("output_folder/Heterozygous_regions_ALL.bed").
- A histogram of the coverage distribution, to help the user decide the expected coverage value (see below).
- Graphs of the coverage along each sequences of the assembly (see below).
- The results of the blast between the heterozygous regions ("output_folder/All_Blasts_scaffolds_coord.tab").

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
If your genome is heterozygous you should obtain two peaks (see graph below):
![alt text](https://bitbucket.org/MCorentin/hetdect/src/master/figures/Histogram_coverage_R.png "Histogram of coverage")

The second peak corresponds to the homozygous regions, and the value on the x-axis for the maximum of this peak corresponds to the homozyous coverage.

### Trying different blast thresholds:

You may want to try different blast thresholds as this will modify the sensitivity of DupLess. We implemented the "-s/--skip_het_detection" parameter for this purpose, so you do not have to run the whole pipeline again. Or if you already have a file with the position of the sequences where you expect duplications.

You can use the bed file produced by a previous run of DupLess with the "-s/--skip_het_detection" option:
```
     python DupLess.py -s Heterozygous_regions_ALL.bed -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [new_output_folder]
```

---

## Future work

- Support long reads to realign to duplicated regions and detect misassemblies
- Possibly: add Mummer to check how the regions around the duplications align to each other.
- Add an option to remove whole contigs if the blast hit span > threshold% of the total length.
- Improve speed.
- Flag regions with half the coverage but no blast hits.


## Performances

For a genome of 1Gbp, DupLess requires ~50 Gb of RAM.

For running time, the bottleneck is the pairwise blasting of heterozygous regions, an assembly with a lot of regions will increase the running time significantly.
