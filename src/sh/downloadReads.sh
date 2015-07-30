#!/bin/bash

#SBATCH --job-name dlLibs
#SBATCH --ntasks 1
#SBATCH --output /tmp/dl.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --mail-type=ALL

# script for downloading publicly available reads for tomato and arabidopsis

# tomato

outdirTom="data/reads/tomato"
if [[ ! -d $outdirTom ]]; then
	mkdir -p $outdirTom
fi

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept wt_sim*.tgz ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept "wt_tm*.tgz" ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept wt_fm*.tgz ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

# arabidopsis

outdirAt="data/reads/at"
if [[ ! -d $outdirAt ]]; then
	mkdir -p $outdirAt
fi

# IM_R1
wget --output-document "$outdirAt"/IM_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351334/ERR351334.fastq.gz &

# IM_R2

wget --output-document "$outdirAt"/IM_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351333/ERR351333.fastq.gz &

# FM_R1

wget --output-document "$outdirAt"/FM_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351336/ERR351336.fastq.gz &

# FM_R2

wget --output-document "$outdirAt"/FM_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351331/ERR351331.fastq.gz &

wait

# email output
NOW="$(date)"
MESSAGE="$(cat /tmp/dl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID submitted at $THEN is finished.
	Concatenated stdout files:
	$MESSAGE
_EOF_

rm /tmp/dl."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out

exit 0
