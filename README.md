README.txt
Authors: Arthur Hsu, Jan Schroeder, Anthony T Papenfuss 
Date: 27/3/2013

A driver script (developed in Python2.7) - Socrates, can be found in this directory.
For a test example run with the script, download the test data from 
http://bioinf.wehi.edu.au/socrates/test_data.tar.gz, extract in the script directory and execute:
./Socrates all --botie2_db test_data/data/bowtie2_db/ecoli1665_bt2 test_data/data/random_ecoli_63_s.bam

Additional documentation is currently being developed. 

The Socrates package contains several Java programs (developed in JDK1.6).
Each program is designed to process data for a specific stages of analysis.

A Python driver script - "Socrates", for retaining cross-platform compatibility of Java,
is included to, is included to execute the programs. 

To use Socrates without the driver script, Java class path needs to be set:

socrates=`dirname $0`/
libs=${socrates}lib/sam-1.77.jar:${socrates}lib/commons-lang3-3.1.jar:${socrates}lib/commons-cli-1.2.jar:${socrates}lib/picard-1.85.jar:${socrates}lib/snappy-java-1.0.3-rc3.jar
java -Xmx4g -cp ${socrates}bin:$libs net.wehi.socrates.[PROG] [OPTIONS]`

where PROG is one of BamStratifier, RealignmentBAM, RealignmentClustering and AnnotatePairedClusters.

While default values have worked satisfactorily in our simulated and real
cancer genome sequencing datasets, users should set program parameters in the
driver scripts appropriately for their own data. Full lists of program
parameters are provided in the following sections, together with a discussion
on the impact of changing them where applicable.


1.1. Preprocess BAM File 

Stratifies the original BAM file into...

usage: Socrates preprocess [options] alignment bam

-b,	--base-quality <score>	Minimum average base quality score of soft clipped sequence  [default: 5] 

-h,	--help 			print this message 

-k,	--keep-duplicate	keep duplicate reads [default: false]

-l,	--long-sc-len <length>	Length threshold of long soft-clip [default: 25 (bp)]

-p,	--percent-id <pid> 	Minimum alignment percent identity to reference [default: 95 (%)]

-q,	--min-mapq <mapq> 	Minimum alignments mapq [default: 5]

-t,	--threads <threads> 	Number of threads to use [default: 1]

-v,	--verbose		be verbose of progress

    
Minimum base quality option:
 
A reasonable threshold helps removing low quality
soft clips that could lead to erroneous breakpoint calls.

Long soft clip length: 

Studies have shown that longer the sequences the more
likely they can be uniquely placed in a genome. In an early study, it is
demonstrated that while percentage of unique mapping improves with increasing
read length, the rate of gain di- minishes past 25nt ( 80% at 25nt and 90% at
40nt). If value for this parameter is too low, many non-unique soft clips will
be produced and impact on system requirement, processing time and reliability
of results downstream. On the other hand, too high the value results in low
number of long soft clips and hence risk of missing breakpoints.
Percent identity: We often observe higher-than-expected base mismatch rate for
reads in satellite, centromeric and telomeric regions where correctness of
alignments can be con- tentious. Minimum percent identity threshold, which is
equivalent to maximum allowable mismatch rate, can greatly reduce these
erroneous alignments.

Minimum mapping quality:
 
Higher mapping quality, while may not guarantee
unique align- ment, is sufficient to exclude multi-mapping anchor alignments
from further analysis for Bowtie2 and BWA aligned reads.


1.2. Process the re-alignment BAM file

usage: Socrates realignment [options] input_bam output_bam

input_bam 	Re-aligned soft clip BAM file. Use “-” to accept input from stdin

output_bam 	Output re-alignment BAM with anchor info merged

anchor_info	Anchor info file produced by BAMStratifier


This program merges soft clip re-alignment BAM file with anchor alignment
information. The program has built-in sorting mechanism and therefore can take
unsorted, raw re- alignment output from aligner. While the program accepts
input BAM file from standard input channel, this requires more system memory
for buffering.


1.3. Predict rearrangements

usage: Socrates predict [options] realigned_sc_bam short_sc_bam metrics_file

-f,	--flank <flank> 		Size of flank for promiscuity filter [default: 50 (bp)]

-h,	--help 				print this message 

-i,	--ideal-only 			Use only proper pair 5’ SC and anomalous pair 3’ SC [default: false]

-l,	--long-sc-len <length> 		Length threshold of long soft-clip [default: 25 (bp)]

-m,	--promiscuity <threshold>	Exclude cluster if more than [promiscuity] clusters within [flank]bp of a break

-p,	--percent-id <pid>		Minimum realignment percent identity to reference [default: 95 (%)]

-q,	--min-mapq <mapq> 		Minimum realignments mapq [default: 2]

-c,	--short-sc-cluster 		Search for short soft clip cluster support for unpaired clusters 

-s,	--max-support <support> 	Maximum realignment support to search for short SC cluster [default: 30]

-t,	--threads <threads> 		Number of threads to use [default: 3]

-v,	--verbose			be verbose of progress
     
The rearrangement predictor is the main part of the algorithm. It clusters the 
split reads, and then pairs clusters to form the output. 
There are two files generated: the paired and unpiared outputs. The paired output
contains the best results, while the unpaired contains any soft-clips that have 
realigned anywhere else in the process re-alignments stage of the algorithm.
The unpaired results are included for completeness as they can contain some
useful information, but overall this output is comprised of false positives 
due to mapping errors and other artefacts. 
The paired output contains various columns of information that describe the 
location of the break point and the level of support. The column labelled 
'BP_condition' describes the nature of the fusion event. It can take any of 
five values:
1. Blunt-end joining: the most straight forward case of a clean join (none
of the below).
2. Micro-homology: Xbp homology found! (XXX): the two joined regions are 
identical for X bases across the break. Therefore the true location of the 
breakpoint is only known within those boundaries.
3. Inserted sequence: XXX: There is a short bit of sequence inserted in 
between the two loci of the fusion. The sequence is either untemplated or 
from somewhere else in the genome (but too short to map).
4. unequal distances of realigned breakpoint to anchor breakpoint: X v Y: 
In this case the realignment and anchor loci of the two paired clusters do 
not support the exact same coordinate for a fusion (|X-Y| indicates the 
difference). This is usually due to mis-mappings.
5. Unequal inserted sequence: XXX v Y: an insert occurs as above, but 
Socrates was unable to determine the exact sequence. One of the values 
should contain the correct sequence.
6. short SC cluster: In most experiments the most prevalent type, yet the
least trustworthy. Only one side of the breakpoint is supported by realiged
split reads, the other by short bits of soft-clipped sequence only. This sort
of cluster pairing makes Socrates very sensitive, but introduces false 
positives.
  
1.4. Annotating rearrangements

usage: Socrates annotate [options] socrates_paired_cluster_output

-n,	--normal <normal>	Socrates paired breakpoint calls for normal sample 

-r,	--repeatmask <file>	UCSC repeat masker track file in BED format, Tabix indexed.


