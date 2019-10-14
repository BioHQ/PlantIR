# PlantIR
Detection of intron retention of palnts


Step1: Building reference file for PlantIR


Example:
perl PlantIR_reference_v1.2.pl Arabidopsis_thaliana.TAIR10.44.gff3 Ara_PlantIR_reference.txt


Step2: Preparing SAM file for PlantIR


1.Raw RNA-Seq data is downloaded from NCBI SRA and .sra file is coversed to .fastq file using fastq-dumpt of SRA Toolkit
Example:
fastq-dump SRR4048211.sra

2. Adaptor removal and quanlity filtering are performed by Trimmomatic
Example:
java -jar trimmomatic-0.33.jar SE -phred33 -threads 10 SRR4048211.fastq SRR4048211_trim.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36&>SRR4048211_trim_log.txt

3. Sequencing reads are mapped to reference genome by HISAT2
Example:
extract_splice_sites.py Arabidopsis_thaliana.TAIR10.44.gff3 > TAIR10.ss
hisat2-build -p 10 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10_genome_index
hisat2 --known-splicesite-infile TAIR10.ss--dta -t -p 10 -x TAIR10_genome_index -U SRR4048211_trim.fastq -S SRR4048211_accepted_hits.sam &> SRR4048211_alignment_summary.txt


Step3: Running PlantIR


perl PlantIR_main.pl Ara_PlantIR_reference.txt [read_length] [thread] [input file (SAM)] [prefix of output file] [minimum of read counts mapped to junctions] [threshold of PSI]

Example:
perl PlantIR_main.pl Ara_PlantIR_reference.txt 100 10 SRR4048211_accepted_hits.sam SRR4048211 8 0.1

If multithreading is not supported by your Perl:

perl PlantIR_main_single.pl Ara_PlantIR_reference.txt [read_length] [input file (SAM)] [prefix of output file] [minimum of read counts mapped to junctions] [threshold of PSI]

Example:
perl PlantIR_main_single.pl Ara_PlantIR_reference.txt 100 SRR4048211_accepted_hits.sam SRR4048211 8 0.1
