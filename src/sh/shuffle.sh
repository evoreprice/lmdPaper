#!/bin/bash

#SBATCH --job-name shuffle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output /tmp/shuffle.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

# make output directory
outdir="output/shuffle-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
version="$(STAR --version)"
cat -t <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	gff2bed version,$(convert2bed --version)

_EOF_


# 1. convert gffs to BED -----------------------------------------------------------

# download osa.gff3 miRBase miRNAs

wget --output-document $outdir/osa.gff ftp://mirbase.org/pub/mirbase/CURRENT/genomes/osa.gff3 &

# download rice_osa1r7_rm.gff3 repeat mask

wget --output-document $outdir/rice_osa1r7_rm.gff3.gz ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/rice_osa1r7_rm.gff3.gz &

# download tRNA/rRNA annotation from Rap-DB (22/5/15)

wget --output-document $outdir/irgsp1_rRNA_tRNA.gff.gz http://rapdb.dna.affrc.go.jp/download/archive/irgsp1_rRNA_tRNA.gff.gz &

# gunzip archives after downloading
wait
gunzip $outdir/rice_osa1r7_rm.gff3.gz 
gunzip $outdir/irgsp1_rRNA_tRNA.gff.gz

# extract "gene" features from gtf

grep "gene"  data/genome/Osativa_204_v7.0.gene_exons.gff3 | gff2bed > $outdir/Osativa_204_v7.0.gene.bed

gff2bed < $outdir/rice_osa1r7_rm.gff3 > $outdir/rice_osa1r7_rm.bed
gff2bed < $outdir/osa.gff3 > $outdir/osa.bed

