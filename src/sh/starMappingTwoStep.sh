#!/bin/bash

#SBATCH --job-name star2s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --output /tmp/star.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --mail-type=ALL

THEN="$(date)"
echo -e "[ "$(date)": Two-step mapping with STAR ]"

# choose the most recent STAR index
shopt -s nullglob
folders=(output/star-index*)
shopt -u nullglob
# stop if there is no STAR index
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No STAR index found ]"
	exit 1
fi
star_index_dir="${folders[-1]}"
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# choose the most recent cutadapt output
shopt -s nullglob
folders=(output/cutadapt*)
shopt -u nullglob
# stop if there is no cutadapt folder
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
cutadapt_dir="${folders[-1]}"
echo -e "[ "$(date)": Using cutadapt folder $cutadapt_dir ]"

# make today's output directory
outdir="$cutadapt_dir/STAR-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$(STAR --version)
	STAR index,$star_index_dir
	cutadapt folder,$cutadapt_dir
_EOF_

# catch SIGTERM etc
clean_up() {
	echo -e "[ "$(date)": Script aborted ]"
	# email output
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Job $SLURM_JOBID aborted" tom
	Job $SLURM_JOBID submitted at $THEN was aborted.
	
	Concatenated stdout files:

	$(cat /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out)
_EOF_
	mv /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# load genome
echo -e "[ "$(date)": Loading genome into shared memory ]"
cmd="STAR --runThreadN 6 --genomeDir $star_index_dir --genomeLoad LoadAndExit --outFileNamePrefix $outdir/gLoad."
srun --ntasks=1 --exclusive --cpus-per-task=6 $cmd

# STAR options
OPTIONS="--runThreadN 6 --genomeDir "$star_index_dir" --outSAMtype BAM Unsorted --alignIntronMax 5000 --outSJfilterReads Unique --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outSJfilterIntronMaxVsReadN 5000 --readFilesCommand zcat"

### STEP 1

# make directory for step 1
if [[ ! -d $outdir/step1 ]]; then
	mkdir -p $outdir/step1
fi

echo -e "[ "$(date)": Submitting step 1 mapping jobs ]"

# find the read files
shopt -s nullglob
fastq_files=("$cutadapt_dir/*.fastq.gz")
shopt -u nullglob
for read_file in $fastq_files
do
	n=$(basename $read_file)
	library_name=${n:0:4}
	cat <<- _EOF_
	[ $(date): Submitting STAR run ]
	library_name:   $library_name
	read_file:      $read_file	
_EOF_
	cmd="STAR $OPTIONS --genomeLoad LoadAndKeep --readFilesIn $read_file --outFileNamePrefix $outdir/step1/$library_name."
	srun --output $outdir/step1/$library_name.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &	
done

echo -e "[ "$(date)": Waiting for step 1 jobs to finish ]"
wait

echo -e "[ "$(date)": Step 1 finished. Removing index from memory ]"
srun --exclusive --ntasks=1 --cpus-per-task=6 \
	STAR --runThreadN 6 --genomeDir $star_index_dir --genomeLoad Remove --outFileNamePrefix $outdir/gRem.

### STEP 2

echo -e "[ "$(date)": Submitting step 2 mapping jobs ]"

# remove step 1 bamfile (waste of space)
rm "$outdir"/step1/*.out.bam

# get the SJ.out.tab files from step 1
sjTabs=("$outdir/step1/*.SJ.out.tab")
cat <<- _EOF_
	[ SJ.out.tab files from step 1 ]
	$(for tab in $sjTabs; do echo $tab; done)
_EOF_

for read_file in $fastq_files
do
	n=$(basename $read_file)
	library_name=${n:0:4}
	cat <<- _EOF_
	[ $(date): Submitting step 2 STAR run ]
	library_name:   $library_name
	read_file:      $read_file	
_EOF_
	cmd="STAR $OPTIONS --genomeLoad NoSharedMemory --sjdbFileChrStartEnd $sjTabs --quantMode GeneCounts --readFilesIn $read_file --outFileNamePrefix $outdir/$library_name."
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &	
done

echo -e "[ "$(date)": Waiting for step 2 jobs to finish ]"

wait
echo -e "[ "$(date)": Jobs finished. Tidying up ]"
# email output
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job $SLURM_JOBID finished" tom
	Job $SLURM_JOBID submitted at $THEN is finished.

	Job log:
	$(cat /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out)
_EOF_

mv /tmp/star."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0