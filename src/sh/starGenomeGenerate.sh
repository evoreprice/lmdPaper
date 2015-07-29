#!/bin/bash

#SBATCH --job-name stargg
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=6
#SBATCH --output /tmp/stargg.%N.%j.out
#SBATCH --open-mode=append

THEN="$(date)"
echo -e "[ "$(date)": Genome generation with STAR ]"

# make output directory
outdir="output/star-index-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# catch sigkill
clean_up() {
	echo -e "[ "$(date)" : Script aborted ]"
	# email output
	NOW="$(date)"
	MESSAGE="$(cat /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" aborted" tom
	Job $SLURM_JOBID submitted at $THEN was aborted.
	Concatenated stdout files:
	$MESSAGE
_EOF_
	mv /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# set parameters
genomeFastaFiles="data/genome/Osativa_204_v7.0.fa"
sjdbGTFfile="data/genome/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
sjdbGTFtagExonParentTranscript="oId"
sjdbGTFtagExonParentGene="gene_name"
sjdbOverhang=49
outFileNamePrefix="$outdir/"

# log metadata
version="$(STAR --version)"
cat -t <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$version
	genomeFastaFiles,$genomeFastaFiles
	sjdbGTFfile,$sjdbGTFfile
	sjdbOverhang,$sjdbOverhang
_EOF_

echo -e "[ "$(date)": Submitting job ]\ngenomeFastaFiles:\t\t$genomeFastaFiles\nsjdbGTFfile:\t\t\t$sjdbGTFfile\nsjdbGTFtagExonParentTranscript:\t$sjdbGTFtagExonParentTranscript\nsjdbOverhang:\t\t\t$sjdbOverhang"

cmd="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript --sjdbGTFtagExonParentGene $sjdbGTFtagExonParentGene --sjdbOverhang $sjdbOverhang --outFileNamePrefix $outFileNamePrefix"

srun --output $outdir/stargg.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &

echo -e "[ "$(date)": Waiting for jobs to finish ]"
wait
echo -e "[ "$(date)": Jobs finished, tidying up ]"

# email output
NOW="$(date)"
MESSAGE="$(cat /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID submitted at $THEN is finished.
	Concatenated stdout files:
	$MESSAGE
_EOF_

mv /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0