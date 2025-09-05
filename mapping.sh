#!/bin/bash

########################################################
## Vanessa cardui samples: Map, sort and index bam files
########################################################

module load bioinfo-tools samtools bwa bcftools vcftools java cutadapt

ulimit -c unlimited

WD="/proj/uppstore2017185/b2014034_nobackup/Aurora/PSMC_cardui/bam_files/"
REF_DIR="/proj/uppstore2017185/b2014034_nobackup/Aurora/Assemblies/Darwin_tree_of_life/"
SAMPLEID="SAMPLE_A"
SAMPLELETTER="A"
DATA="/proj/uppstore2017185/b2014034/private/vanessa_cardui_project/raw_data/OEB_Insect_Genome_Sequencing/$SAMPLEID/TRUSEQ_DNA_PCR_FREE_LIBRARY_$SAMPLEID/"


cd $TMPDIR

cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L001_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L001_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L001_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L001_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L002_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L002_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L002_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L002_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L003_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L003_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L003_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L003_R2_001.fastq.gz &
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -o $TMPDIR/TS_$SAMPLELETTER\_S1_L004_R1_001_trimmed.fastq.gz -p $TMPDIR/TS_$SAMPLELETTER\_S1_L004_R2_001_trimmed.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L004_R1_001.fastq.gz $DATA/TS_$SAMPLELETTER\_S1_L004_R2_001.fastq.gz

wait


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

samtools merge -cf $WD/$SAMPLEID\_dtol_sorted.bam $TMPDIR/$SAMPLEID\_1_sorted.bam $TMPDIR/$SAMPLEID\_2_sorted.bam $TMPDIR/$SAMPLEID\_3_sorted.bam $TMPDIR/$SAMPLEID\_4_sorted.bam

samtools index $WD/$SAMPLEID\_dtol_sorted.bam


#########################################################
## Vanessa kershawi sample: Map, sort and index bam files
#########################################################

module load bioinfo-tools Stampy samtools BEDTools seqtk bcftools gvcftools QualiMap
module load java/sun_jdk1.8.0_151

ref="/proj/uppstore2017185/b2014034_nobackup/Aurora/dnds/reference_genomes/GCA_905220365.1_ilVanCard2.1_genomic_chroms"
wd="/proj/uppstore2017185/b2014034_nobackup/Aurora/dnds/00_assemblies_stampy_output"

PE1='/proj/uppstore2017185/b2014034_nobackup/Orazio/resequecing_data/P13011_106_S215_L002_R1_001.fastq.gz'
PE2='/proj/uppstore2017185/b2014034_nobackup/Orazio/resequecing_data/P13011_106_S215_L002_R2_001.fastq.gz'

cd $TMPDIR

species="kershawi"

/sw/apps/bioinfo/Stampy/1.0.32/rackham/stampy.py -g $ref -h $ref --readgroup=ID:group1,SM:sample1,PL:illumina,LB:lib1,PU:unit1 --threads=10 --substitutionrate=0.0406 -o $wd/$species.sam -M $PE1 $PE2 
cp $TMPDIR/$species.sam $wd/$species.sam

samtools view -bt $ref.fa.fai -o $wd/$species.bam $wd/$species.sam 
samtools sort -o $wd/$species\_sorted.bam $wd/$species.bam
samtools index $wd/$species\_sorted.bam
qualimap bamqc -bam $wd/$species\_sorted.bam -outdir $wd -outfile $species\_qualimap.pdf -nt 10 --java-mem-size=5G


##########################################################
# Common to all samples: MarkDuplicates and IndelRealigner
##########################################################


## MarkDuplicates with GATK

java -Xmx32g -jar /sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar ValidateSamFile I=$wd/$species\_sorted.bam
java -Xmx32g -jar /sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar CreateSequenceDictionary REFERENCE=$ref.fa OUTPUT=$ref.fa.dict 
java -Xmx32g -jar /sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar MarkDuplicates INPUT=$wd/$species\_sorted.bam OUTPUT=$wd/$species.bam.dedup.bam METRICS_FILE=$wd/$species.bam.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
samtools index $wd/$species.bam.dedup.bam

## IndelRealigner with GATK

java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I $wd/$species.bam.dedup.bam -R $ref.fa -T RealignerTargetCreator -o $wd/$species.intervals
java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -I $wd/$species.bam.dedup.bam -R $ref.fa -T IndelRealigner --filter_bases_not_stored -o $wd/$species\_realigned.bam -targetIntervals $wd/$species.intervals
samtools index $wd/$species\_realigned.bam
qualimap bamqc -bam $wd/$species\_realigned.bam -outdir $wd -outfile $species\_realigned_qualimap.pdf -nt 10 --java-mem-size=5G
