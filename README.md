express-full-genome
==================

Expand eXpress from RNA-seq to any genome wide high throughput sequencing assay


##Setup

Install python (v2.7), pysam, scipy, numpy and BioPython (python modules), samtools, bowtie, express

Create environmental variables: $BOWTIE_INDEXES and $EXPRESS_FILES

##Steps to run pipeline 

**Bin the entire genome to create the target sequence**

BinnedFasta.py takes the following arguments
* -g: reference genome fasta file
* -l: length of the sequence
* -b: size of the bin

Example:

`` BinnedFasta.py -g mm9.fa -l 50 -b 200 ``

**Map reads to reference genome**

Map the sequences to the original reference genome and save it in SAM/BAM. File sizes can grow quickly when allowing for multiply aligned sequences so BAM files are suggested instead.

Example: 

`` bzip2 -dc reads.fastq.bz2 | bowtie -c -q -v 2 -k 100 -S mm9 - | samtools view -bS - > aligned.k100.bam ``

**Map reference genome alignment coordinates to the target sequence alignment**

BinMapping.py takes the following arguments
* -g: reference genome fasta file
* -l: length of the sequence
* -b: size of the bin
* -r: aligned BAM file
* -o: prefix of the converted BAM file

Example:

`` BinMapping.py -g mm9.fa -r aligned.k100.bam -l 50 -b 200 -o aligned.200bp.k100 ``

**Pass into express**

eXpress will take the converted BAM file and calculate the likelihood for each alignment with relation to the local bin abundances in the target sequences.

Example:

``express -o output_dir --output-align-samp $EXPRESS_FILES/eXpress_200bp_50.mm9.fa aligned.200bp.k100_converted.bam``

**Stream above steps**

To avoid creating multiple large BAM files, it is possible to stream the mapping into eXpress and then save a significantly smaller final file. 

``bzip2 -dc reads.fastq.bz2 | bowtie -c -q -v 2 -aS mm9 - | BinMapping.py -g mm9.fa -l 50 -b 200 -r - -o - | samtools view -h - | express -o --output-dir --output-align-samp -B 1 $EXPRESS_FILES/eXpress_200bp_50.mm9.fa``

**Map target sequence alignments back to the reference genome**

eXpress2wiggle.py takes the following arguments
* -g: reference genome prefix (e.g. mm9)
* -r: eXpress outputted alignment file
* -o: name for all outputted files
* -bo: do you want to output the converted BAM file (y/n)

If you are using a shape based algorithm to identify enriched regions within the genome, there is no need to use the -bo y option. Wiggle and BedGraph files are generated automatically.

Example:

``eXpress2wiggle.py -g mm9 -r output_dir/hits.1.samp.bam -bo y -o Aligned.200bp.samp_wB``

###Disclaimer
This code is still in beta form and some bugs exist. 

