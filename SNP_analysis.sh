#aligns the fastq files to a reference genome
module load BWA/0.7.17-foss-2016b
bwa mem -M -t 16 /lustre1/gwc32007/crypto_fastq_files/crypto_iowaIIfullgenome.fasta /lustre1/gwc32007/crypto_fastq_files/UKP4/SRR6147581_1.fastq.gz /lustre1/gwc32007/crypto_fastq_files/UKP4/SRR6147581_2.fastq.gz > UKP4_aln.sam

#converts the .sam file to a more condensed .bam file
module load SAMtools/1.6-foss-2016b
samtools view -bt ~/crypto_fastq_files/crypto_iowaIIfullgenome.fasta -o UKP3_aln.bam ~/fastq_alignments/UKP3_aln.sam

#sorts the alignments from smallest to largest
samtools sort /lustre1/gwc32007/fastq_alignments/UKP3_aln.bam -o /lustre1/gwc32007/fastq_alignments/UKP3_aln.sorted.bam

#marks the duplicate reads in the alignment

#PBS -S /bin/bash
#PBS -N j_picard
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=25g

cd $PBS_O_WORKDIR

module load picard/2.16.0-Java-1.8.0_144
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates I=/lustre1/gwc32007/fastq_alignments/UKP3_aln.sorted.bam O=/lustre1/gwc32007/fastq_alignments/UKP3_aln.marked_duplicates.bam M=mark
ed_dup_metrics.txt REMOVE_DUPLICATES=false

#replaces read groups in the bam file with one read group

#PBS -S /bin/bash
#PBS -N j_picard
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=25g

cd $PBS_O_WORKDIR

module load picard/2.16.0-Java-1.8.0_144
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=/lustre1/gwc32007/fastq_alignments/UKP3_aln.marked_duplicates.bam O=/lustre1/gwc32007/fastq_alignments/UKP3_aln.dupl.read.
sort.bam SORT_ORDER=coordinate RGLB=lib1 RGPL=illumina RGSM=20 RGPU=unit1 VALIDATION_STRINGENCY=LENIENT

#indexing the reference

#PBS -S /bin/bash
#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR

module load SAMtools/1.6-foss-2016b

time samtools faidx /lustre1/gwc32007/crypto_fastq_files/crypto_iowaIIfullgenome.fasta

time samtools index /lustre1/gwc32007/fastq_alignments/UKP3_aln.dupl.read.sort.bam

#creating dictionary for reference

module load picard/2.16.0-Java-1.8.0_144
#PBS -S /bin/bash
#PBS -N j_picard
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=25g

cd $PBS_O_WORKDIR

module load picard/2.16.0-Java-1.8.0_144
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar CreateSequenceDictionary R=/lustre1/gwc32007/crypto_fastq_files/crypto_iowaIIfullgenome.fasta O=/lustre1/gwc32007/crypto_fastq_files/crypto_iowaIIf
ullgenome.fasta.dict
