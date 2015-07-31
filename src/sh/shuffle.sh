#!/bin/bash

#SBATCH --job-name shuffle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output /tmp/shuffle.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

THEN="$(date)"

# make output directory
outdir="output/shuffle-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	gff2bed version,$(convert2bed --version)
	wgsim version,$(wgsim 2>&1 | sed '3q;d')
	bedtools version,$(bedtools --version)
	tophat2 version,$(tophat2 --version)
	bowtie2-build version,$(bowtie2-build --version | head -n 1)
_EOF_

# 1. convert gffs to BED -----------------------------------------------------------

echo -e "[ $(date): Downloading micro RNA, repeat mask, tRNA and mRNA files ]"

# download osa.gff3 miRBase miRNAs
wget --output-document $outdir/osa.gff3 ftp://mirbase.org/pub/mirbase/CURRENT/genomes/osa.gff3 &

# download rice_osa1r7_rm.gff3 repeat mask
wget --output-document $outdir/rice_osa1r7_rm.gff3.gz ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/rice_osa1r7_rm.gff3.gz &

# download tRNA/rRNA annotation from Rap-DB

# Rap-DB is down, for now use a symlink in the data/ folder to a version downloaded earlier (22/5/15)
#wget --output-document $outdir/irgsp1_rRNA_tRNA.gff.gz http://rapdb.dna.affrc.go.jp/download/archive/irgsp1_rRNA_tRNA.gff.gz &
cp data/irgsp1_rRNA_tRNA.gff $outdir/irgsp1_rRNA_tRNA.gff

# gunzip archives after downloading
wait
gunzip $outdir/rice_osa1r7_rm.gff3.gz 
#gunzip $outdir/irgsp1_rRNA_tRNA.gff.gz

# extract "gene" features from gtf
echo -e "[ $(date): Extracting genes from gff and converting files to bed9 ]"
grep "gene"  data/genome/Osativa_204_v7.0.gene_exons.gff3 | gff2bed > $outdir/Osativa_204_v7.0.gene.bed

# remove genes with wonky coordinates from the tRNA/rRNA database 
sed '282d' $outdir/irgsp1_rRNA_tRNA.gff | sed '537d' | sed '913d' | gff2bed > $outdir/irgsp1_rRNA_tRNA.bed

# convert to bed
gff2bed < $outdir/rice_osa1r7_rm.gff3 > $outdir/rice_osa1r7_rm.bed
gff2bed < $outdir/osa.gff3 > $outdir/osa.bed

# convert bed9 output to bed6 (call R script with location of outdir)
echo -e "[ $(date): Converting bed9 to bed6 ]"
src/R/fixGffBed.R $outdir

# 2. generate bedfile for repeat features -----------------------------------------------------------
echo -e "[ $(date): Generating bedfile from plant repeat database ]"

# download oryza repeat database
wget --output-document $outdir/TIGR_Oryza_Repeats.v3.3_0_0.fsa ftp://ftp.plantbiology.msu.edu/pub/data/TIGR_Plant_Repeats/TIGR_Oryza_Repeats.v3.3

# simulate 1 million reads from plant repeat database
echo -e "[ $(date): Simulating reads from plant repeat database ]"
wgsim -e 0 -1 50 -2 1 -r 0 -R 0 -X 0 -d 0 -s 0 	$outdir/TIGR_Oryza_Repeats.v3.3_0_0.fsa $outdir/TIGR_Oryza_Repeats.wgsim.fq /dev/null

# map simulated reads to reference genome with tophat (might need to build bowtie2 index)
if [[ ! -e "data/genome/Osativa_204_v7.0.1.bt2" ]]; then
	echo -e "[ $(date): Need to build a bowtie2 index ]"
	bowtie2-build data/genome/Osativa_204_v7.0.fa data/genome/Osativa_204_v7.0
fi
echo -e "[ $(date): Mapping simulated reads ]"
tophat2 --output-dir $outdir/tophat data/genome/Osativa_204_v7.0 $outdir/TIGR_Oryza_Repeats.wgsim.fq

# convert tophat accepted hits to bed format
bedtools bamtobed -i $outdir/tophat/accepted_hits.bam > $outdir/accepted_hits.bed6

# 3. combine features and add of slop ----------------------------------------------------------------

# count chromosome lengths
echo -e "[ $(date): Getting chromosome lengths ]"
src/py/seqLength.py data/genome/Osativa_204_v7.0.fa > $outdir/genome.bed

# add slop
echo -e "[ $(date): Adding slop to gtf ]"
bedtools slop -b 1000 -i $outdir/Osativa_204_v7.0.gene.bed6 -g $outdir/genome.bed > $outdir/Osativa_204_v7.0.gene.slop.bed6

# mask uslopped gtf in tmp
mv $outdir/Osativa_204_v7.0.gene.bed6 $outdir/Osativa_204_v7.0.gene.tmp

bed6files=("$outdir/*.bed6")

cat <<- _EOF_
	[ $(date): Merging bed files ]
	$(for item in $bed6files; do echo $item; done)
_EOF_

sort-bed $bed6files	| bedtools merge -i stdin > $outdir/excl.combined.bed

# 4. convert the combined befile to GTF and shuffle ----------------------------------------------------------------

# make a GTF file from the combined bedfile (takes ages)
echo -e "[ $(date): Generating dummy GTF from combined bedfile ]"
src/R/generateCdsLengthsForShuffle.R $outdir

# *finally* shuffle
echo -e "[ $(date): Shuffling the dummy GTF ]"
bedtools shuffle -excl $outdir/excl.combined.bed -chromFirst -chrom -noOverlapping -seed 1 -i $outdir/dummyGff.gff3 -g $outdir/genome.bed > $outdir/Osativa_204_v7.shuffled.gff3

echo -e "[ $(date): Done ]"

# email output
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID started at $THEN has finished.
	Output:

	$(cat /tmp/shuffle."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)
_EOF_

mv /tmp/shuffle."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out $outdir/shuffle."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out

exit 0