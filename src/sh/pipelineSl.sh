#!/bin/bash

#SBATCH --job-name compSl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --output /tmp/compSl.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --mail-type=ALL

# catch sigkill
clean_up() {
	echo -e "[ "$(date)": Script aborted ]"
	# email output
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" aborted" tom
	Job $SLURM_JOBID submitted at $THEN was aborted.
	Output so far:
	$(cat /tmp/compSl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)
_EOF_
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# periodic updates
update_output() {
	echo -e "[ "$(date)": updating output ]"
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Pipeline job "$SLURM_JOBID" running" tom
	Job $SLURM_JOBID submitted at $THEN is running.
	Output so far:
	$(cat /tmp/compSl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)
_EOF_
}

# TOMATO

echo -e "[ "$(date)": Starting pipeline for Solanum lycopersicum ]"

# 1. Generate genome for STAR

echo -e "[ "$(date)": Genome generation with STAR ]"

# make output directory
outdir="output/madsComp/sl/star-index-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# set parameters
genomeFastaFiles="data/sl/genome/Slycopersicum_225_iTAGv2.40.fa"
sjdbGTFfile="data/sl/genome/Slycopersicum_225_iTAGv2.3.gene_exons.cuffcomp.gtf"
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

cat <<- _EOF_
[ $(date): Generating genome ]
              genomeFastaFiles:    $genomeFastaFiles
                   sjdbGTFfile:    $sjdbGTFfile
sjdbGTFtagExonParentTranscript:    $sjdbGTFtagExonParentTranscript
                  sjdbOverhang:    $sjdbOverhang
_EOF_

cmd="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript --sjdbGTFtagExonParentGene $sjdbGTFtagExonParentGene --sjdbOverhang $sjdbOverhang --outFileNamePrefix $outFileNamePrefix"

srun --output $outdir/stargg.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &

echo -e "[ "$(date)": Waiting... ]"
wait

echo -e "[ "$(date)": Genome generation finished ]"
update_output

# 2. Trim adaptors

echo -e "[ "$(date)": Converting reads to .fastq.gz ]"



echo -e "[ "$(date)": Trimming reads with cutadapt ]"

# make today's output directory
outdir="output/madsComp/sl/cutadapt-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$(cutadapt --version)
_EOF_

# parameters
adaptorFwd='TruSeq_indexed_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adaptorRev='TruSeq_universal_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
trim_qualities=20
minimum_length=50


#### THIS IS BORKED, NEED TO START AGAIN FROM HERE TO USE CONVERTED FASTQ.GZ FILES!!!

echo -e "[ "$(date)": Submitting cutadapt jobs ]"
fwd_reads=("data/reads/tomato/*L.tgz")
for readFile in $fwd_reads; do
	fFile="$(basename $readFile)"
	lib_name="${fFile:0:8}"
	rev_reads="data/reads/tomato/$(basename $readFile L.tgz)R.tgz"
	if [[ ! -e "$rev_reads" ]]; then echo -e "Error: rev_reads not found\n[ lib_name ]:\t$lib_name\n[ rev_reads ]:\t$rev_reads"; exit 1; fi
	output="$outdir/$lib_name.R1.fastq.gz"
	paired_output="$outdir/$lib_name.R2.fastq.gz"
	# print some info
	cat <<- _EOF_
	[ $(date): Running cutadapt ]
	     lib_name:    $lib_name
	    fwd_reads:    $readFile
	    rev_reads:    $rev_reads
	       output:    $output
	paired_output:    $paired_output
_EOF_
	# run cutadapt
	cmd="cutadapt -a $adaptor -A $adaptorRev --quality-cutoff $trim_qualities --minimum-length $minimum_length --output=$output --paired-output=$paired_output $readFile $rev_reads"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for adaptor trimming to finish ]"
wait

echo -e "[ "$(date)": Adaptor trimming finished ]"
update_output

# 3. Run 2-step mapping

echo -e "[ "$(date)": Two-step mapping with STAR ]"

# choose the most recent STAR index
shopt -s nullglob
folders=(output/madsComp/sl/star-index*)
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
folders=(output/madsComp/sl/cutadapt*)
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
	library_name=${n:0:5}
	cat <<- _EOF_
	[ $(date): Submitting STAR run ]
	library_name:    $library_name
	   read_file:    $read_file	
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

# remove step 1 bamfiles (waste of space)
rm $outdir/step1/*.out.bam

# get the SJ.out.tab files from step 1
sjTabs=("$outdir/step1/*.SJ.out.tab")
cat <<- _EOF_
	[ SJ.out.tab files from step 1 ]
	$(for tab in $sjTabs; do echo $tab; done)
_EOF_

for read_file in $fastq_files
do
	n=$(basename $read_file)
	library_name=${n:0:5}
	cat <<- _EOF_
	[ $(date): Submitting step 2 STAR run ]
	library_name:    $library_name
	   read_file:    $read_file	
_EOF_
	cmd="STAR $OPTIONS --genomeLoad NoSharedMemory --sjdbFileChrStartEnd $sjTabs --quantMode GeneCounts --readFilesIn $read_file --outFileNamePrefix $outdir/$library_name."
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &	
done

echo -e "[ "$(date)": Waiting for step 2 jobs to finish ]"
wait

echo -e "[ "$(date)": Step 2 jobs finished ]"

echo -e "[ "$(date)": Pipeline finished! Phew! Emailing output. ]"

# email output
cat <<- _EOF_ | mail -s "[Tom@SLURM] Pipeline job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID submitted at $THEN has finished.
	Output:

	$(cat /tmp/compSl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)
_EOF_

echo -e "[ "$(date)": Storing output ]"
echo "output/madsComp/at/compSl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out"

mv /tmp/compSl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out output/madsComp/at/

echo -e "[ "$(date)": Done, exiting ]"
exit 0
