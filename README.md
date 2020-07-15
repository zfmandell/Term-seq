# Term-seq
Directions and associated custom python scripts to obtain total 3' ends from raw .fastq.gz file obtained by conducting Term-seq on the Illumina Sequencing Platform (single-end sequencing)

requires in total :

(these are the versions of each package that I used, but other versions likely work too):

Python 3.8.3 (with associated modules: biopython v. 1.77, numpy v. 1.19.0, scipy v. 1.5.1), Trimmomatic v. 0.38, Cutadapt v. 1.16, BWA-MEM v.0.7.12-r1034, samtools v. 0.1.19-44428cd, bedtools v. 2.26.0

to run :

1) run trimmomatic in SE mode
java -jar (trimmomatic location) SE -phred33 (input file name .fastq.gz) (output file name .fastq.gz) ILLUMINACLIP:(Illumina adaptor sequences .fasta):2:30:7 MINLEN:20

This yields RNA-seq dataset

2) run cutadapt
cutadapt -g TAGCTCATCGACTGGATCTCAGTGTCTCATT -O 8 -m 20 --discard-untrimmed --info-file=(information file name .txt) -o (output reads file name .fq) (input reads file name, output of step 1 .fastq.gz)

This yields all reads containing custom oligo sequence at 5' end of R1

3) run bwa-mem -> sam to bam 

bwa mem (reference genome .fasta file) (custom oligo containing reads, output of step 2 .fq) | samtools view -b - > (output aligned reads name .bam)

This yields sorted .bam file of mapped custom oligo containing reads

4) run samtools view
samtools view -b -F 4 -f 16 (sorted mapped reads, output of step 3 .bam) > (mapped reads, split by strand (fwd) fwd.bam)
samtools view -b -F 4 -F 16 (sorted mapped reads, output of step 3 .bam) > (mapped reads, split by strand (rev) rev.bam)

this will yield all reads mapped to fwd strand and rev strand

5) samtools sort
samtools sort -o (output file name, sorted and mapped by strand (fwd) sorted.fwd.bam) (mapped reads split by strand (fwd) fwd.bam)
samtools sort -o (output file name, sorted and mapped by strand (rev) sorted.rev.bam) (mapped reads split by strand (rev) rev.bam)

This will sort by genomic coordinate the .bam files that were split by strand in step 4 

6) run bedtools
bedtools  genomecov -d  -ibam (name of sorted and split bam file, output of step 5 (fwd) sorted.fwd.bam)  -g (reference genome .fasta file used in step 2)  > (output file name fwd.cov)
bedtools  genomecov -d  -ibam (name of sorted and split bam file, output of step 5 (rev) sorted.rev.bam)  -g (reference genome .fasta file used in step 2)  > (output file name rev.cov)

This will get the coverage at each position in the genome

7) run positional_slopy.py (custom python script #1)
positional_slope.py plus (+optional -window_size x) (run in directory with fwd strand .cov files)
positional_slope.py minus (+optional -window_size x) (run in directory with rev strand .cov files)

This will yield genome-wide Cv .bedgraph files for each strand

8) run local_peak.py (custom python script #2)
local_peak.py plus (+optional -min_delta x) (+optional -fasta_out <true/false>) (+optional -fasta <.fasta file used for alignment>) (+optional -fasta_window x) (run in directory with fwd strand .bedgraph files, output of step 7)
local_peak.py minus (+optional -min_delta x) (+optional -fasta_out <true/false>) (+optional -fasta <.fasta file used for alignment>) (+optional -fasta_window x) (run in directory with rev strand .bedgraph files, output of step 7)

This will yield all 3' ends that have Cv values above threshold (-min_delta, default = 10) for each strand

9) run rev_inverse.py
rev_inverse.py (rev strand called 3' ends, output of step 8 rev.bedgraph) (output file inverse.rev.bedgraph)

This will yield inversed (Cv values will be negative) .bedgraph file for rev strand 3' ends, used for IGV viewing

10) run cat
cat (fwd strand called 3' end file, output of step 8) (inverse rev strand called 3' end file, output of step 9) > (final 3' end file name .bedgraph)

This will yield final 3' end file for analysis and viewing on IGV







