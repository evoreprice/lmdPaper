#!/bin/bash

set -e

outdirTom="data/reads/sl"
if [[ ! -d $outdirTom ]]; then
	mkdir -p $outdirTom
fi

# handle waiting
FAIL=0
fail_wait() {
for job in $(jobs -p); do
  wait $job || let "FAIL+=1"
done
if [[ ! "$FAIL" == 0 ]]; then
  exit 1
fi
}

echo -e "[ "$(date)": Downloading tomato reads ]"

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept "wt_sim*.tgz" ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept "wt_tm*.tgz" ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

wget --no-parent --directory-prefix "$outdirTom" --recursive --no-directories \
	--accept "wt_fm*.tgz" ftp://ftp.solgenomics.net/user_requests/LippmanZ/public_releases/by_species/Solanum_lycopersicum/transcripts/

#first decompress
echo -e "[ "$(date)": Extracting reads ]"
tgzFiles=("$outdirTom/*.tgz")
FAIL=0
for tgz in $tgzFiles; do
	bn=$(basename $tgz .tgz)
	cmd="tar --gunzip --get --verbose --to-stdout --file $tgz"
	srun --ntasks=1 --cpus-per-task=1 --output=$outdirTom/$bn.fastq --error=/dev/null --exclusive $cmd &
done
fail_wait

# then recompress
echo -e "[ "$(date)": Recompressing reads with gzip ]"
FAIL=0
$(find "$outdirTom" -name "*fastq" -type f -exec bash -c 'srun --exclusive --ntasks=1 --cpus-per-task=1 gzip --best {} &' \;)
fail_wait

# then delete
tgzFiles=("$outdirTom/*.tgz")
for tgz in $tgzFiles; do
	bn=$(basename $tgz .tgz)
	if [[ -s "$outdirTom"/"$bn".fastq.gz ]]; then
		cat <<- _EOF_
		[ $(date): Done ]
		          Read file: $outdirTom/$bn.fastq.gz
		Old file (deleting): $tgz
_EOF_
		rm "$tgz"
	fi
done

exit 0
