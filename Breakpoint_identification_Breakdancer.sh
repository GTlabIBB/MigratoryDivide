### Breakdancer pipeline ###

# First we perform read mapping against the reference genome 
WD="/proj/uppstore2017185/b2014034_nobackup/Aurora/PSMC_cardui/bam_files/"
REF_DIR="/proj/uppstore2017185/b2014034_nobackup/Aurora/Assemblies/Darwin_tree_of_life/"
SAMPLEID="SAMPLE_A"
SAMPLELETTER="A"
DATA="/proj/uppstore2017185/b2014034/private/vanessa_cardui_project/raw_data/OEB_Insect_Genome_Sequencing/$SAMPLEID/TRUSEQ_DNA_PCR_FREE_LIBRARY_$SAMPLEID/"

cd $TMPDIR

# Remove adapters from reads. We have four lanes of sequencing
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L001_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L001_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L001_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L001_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L002_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L002_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L002_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L002_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L003_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L003_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L003_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L003_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L004_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L004_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L004_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L004_R2_001.fastq.gz

wait
 
# Mapping
bwa mem -t 10 -M $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna -R  "@RG\tLB:Lib1\tID:1\tSM:SAMPLE_B\tPL:ILLUMINA" $TMPDIR/TS_A_S1_L001_R1_001_trimmed.fastq.gz $TMPDIR/TS_A_S1_L001_R2_001_trimmed.fastq.gz | samtools view -bt $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna > $TMPDIR/$SAMPLEID\_1.bam
samtools sort -o $TMPDIR/$SAMPLEID\_1_sorted.bam $TMPDIR/$SAMPLEID\_1.bam
samtools index $TMPDIR/$SAMPLEID\_1_sorted.bam

bwa mem -t 10 -M $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna -R  "@RG\tLB:Lib1\tID:1\tSM:SAMPLE_B\tPL:ILLUMINA" $TMPDIR/TS_A_S1_L002_R1_001_trimmed.fastq.gz $TMPDIR/TS_A_S1_L002_R2_001_trimmed.fastq.gz | samtools view -bt $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna > $TMPDIR/$SAMPLEID\_2.bam
samtools sort -o $TMPDIR/$SAMPLEID\_2_sorted.bam $TMPDIR/$SAMPLEID\_2.bam
samtools index $TMPDIR/$SAMPLEID\_2_sorted.bam

bwa mem -t 10 -M $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna -R  "@RG\tLB:Lib1\tID:1\tSM:SAMPLE_B\tPL:ILLUMINA" $TMPDIR/TS_A_S1_L003_R1_001_trimmed.fastq.gz $TMPDIR/TS_A_S1_L003_R2_001_trimmed.fastq.gz | samtools view -bt $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna > $TMPDIR/$SAMPLEID\_3.bam
samtools sort -o $TMPDIR/$SAMPLEID\_3_sorted.bam $TMPDIR/$SAMPLEID\_3.bam
samtools index $TMPDIR/$SAMPLEID\_3_sorted.bam

bwa mem -t 10 -M $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna -R  "@RG\tLB:Lib1\tID:1\tSM:SAMPLE_B\tPL:ILLUMINA" $TMPDIR/TS_A_S1_L004_R1_001_trimmed.fastq.gz $TMPDIR/TS_A_S1_L004_R2_001_trimmed.fastq.gz | samtools view -bt $REF_DIR/GCA_905220365.1_ilVanCard2.1_genomic.fna > $TMPDIR/$SAMPLEID\_4.bam
samtools sort -o $TMPDIR/$SAMPLEID\_4_sorted.bam $TMPDIR/$SAMPLEID\_4.bam
samtools index $TMPDIR/$SAMPLEID\_4_sorted.bam

# Merge the aligment files into a single file
samtools merge -cf $WD/$SAMPLEID\_dtol.bam $TMPDIR/$SAMPLEID\_1_sorted.bam $TMPDIR/$SAMPLEID\_2_sorted.bam $TMPDIR/$SAMPLEID\_3_sorted.bam $TMPDIR/$SAMPLEID\_4_sorted.bam

# Index the final alignment file
samtools index $WD/$SAMPLEID\_dtol.bam


# Subset the inversion region
bedtools intersect -sorted -a $WD/$SAMPLEID\_dtol.bam -b chr8_inverted_region.bed > sample_A_chr8inv.bam

# Create the stats file needed for the execution of Breakdancer
perl /usr/local/bin/breakdancer/perl/bam2cfg.pl -g -h sample_A_chr8inv.bam > SAMPLE_A_chr8inv.cfg

# Run Breakdancer
/usr/local/bin/breakdancer/build/bin/breakdancer-max -o $ref SAMPLE_A_chr8inv.cfg > SAMPLE_A_breakd_output