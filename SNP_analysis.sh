#indexes reference genome

#PBS -S /bin/bash
#PBS -N j_bwa
#PBS -l walltime=128:00:00
#PBS -l nodes=1:ppn=1:AMD
#PBS -q batch

cd $PBS_O_WORKDIR

module load BWA/0.7.17-foss-2016b
time bwa index /scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta

#aligns reads to reference genome

#PBS -S /bin/bash
#PBS -N j_bwa
#PBS -l walltime=128:00:00
#PBS -l nodes=1:ppn=1:AMD
#PBS -q batch

cd $PBS_O_WORKDIR

module load BWA/0.7.17-foss-2016b
time bwa mem -M -t 16 /scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta /scratch/gwc32007/crypto_fastq/UKP2/fastq#1 /scratch/gwc32007/crypto_fastq/UKP2/fastq#2 > /scratch/gwc32007/fastq_alignments/UKP2_aln.sam

#converts the .sam file to a more condensed .bam file

#PBS -S /bin/bash
#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR

module load SAMtools/1.6-foss-2016b
time samtools view -bt /scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta -o /scratch/gwc32007/fastq_alignments/UKP2_aln.bam /scratch/gwc32007/fastq_alignments/UKP2_aln.sam

#sorts the alignments from smallest to largest

#PBS -S /bin/bash
#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR

module load SAMtools/1.6-foss-2016b

time samtools sort /scratch/gwc32007/fastq_alignments/UKP2_aln.bam -o /scratch/gwc32007/fastq_alignments/UKP2_aln.sorted.bam

#marks the duplicate reads in the alignment

#PBS -S /bin/bash
#PBS -N j_picard
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=25g

cd $PBS_O_WORKDIR

module load picard/2.16.0-Java-1.8.0_144
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates I=/scratch/gwc32007/fastq_alignments/UKP2_aln.sorted.bam O=/scratch/gwc32007/fastq_alignments/UKP2_aln.sorted_duplicates.bam M=marke
d_dup_metrics.txt REMOVE_DUPLICATES=false

#replaces read groups in the bam file with one read group

#PBS -S /bin/bash
#PBS -N j_picard
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=25g

cd $PBS_O_WORKDIR

module load picard/2.16.0-Java-1.8.0_144
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=/scratch/gwc32007/fastq_alignments/UKP2_aln.sorted_duplicates.bam O=/scratch/gwc32007/fastq_alignments/UKP2_aln.dupl.read.
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

time samtools faidx /scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta

time samtools index /scratch/gwc32007/fastq_alignments/UKP2_aln.dupl.read.sort.bam

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
time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar CreateSequenceDictionary R=/scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta O=/scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.dict

#find SNP calls in the alignments

#PBS -S /bin/bash
#PBS -N j_gatk
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l mem=2gb

module load GATK/3.8-0-Java-1.8.0_144

java -Xmx4g -jar /usr/local/apps/eb/GATK/3.8-0-Java-1.8.0_144/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scratch/gwc32007/crypto_genomes/cparvum_iowaII_genome.fasta -I /scratch/gwc32007/fastq_alignments/UKP2_aln.dupl.read.sort.bam -o /scratch/gwc32007/SNPGenie/gvcf_CURO/UKP2_output.g.vcf.g
z -ERC GVCF
