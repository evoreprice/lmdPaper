#!/bin/bash

set -e

# catch email address (-e) and password (-p) for jgi passed to bash script
# if the jgi password has special characters this script won't work

while [ "$1" != "" ]; do
	case $1 in
		-e )	shift
				EMAIL=$1
				;;
		-p )	shift
				PASSWORD=$1
				;;
		* )		echo "Bad input"
				exit 1
	esac
	shift
done

# arabidopsis
genome_dir='data/genome/at'
genome_url='http://genome.jgi.doe.gov/Athaliana/download/_JAMO/53112a1d49607a1be0055864/Athaliana_167_TAIR9.fa.gz'
annot_url='http://genome.jgi.doe.gov/Athaliana/download/_JAMO/53112a1a49607a1be005585d/Athaliana_167_TAIR10.gene_exons.gff3.gz'
genome_file="$(basename $genome_url .fa.gz)"
annotation_file="$(basename $annot_url .gff3.gz)"

# open phytozome session
echo -e "[ "$(date)": Signing on to phytozome at JGI ]"
curl https://signon.jgi.doe.gov/signon/create --data-ascii \
	login="$EMAIL"\&password="$PASSWORD" \
	-b "$genome_dir"/cookies -c "$genome_dir"/cookies > /dev/null

# download genomes
if [ ! -d $genome_dir ]; then
	mkdir -p $genome_dir
fi

cat <<- _EOF_
	[ $(date): Downloading genome fasta ]
	$genome_url
_EOF_
curl $genome_url -b "$genome_dir"/cookies -c "$genome_dir"/cookies > $genome_dir/$genome_file.fa.gz
# make sure the file was downloaded
if [[ ! -s "$genome_dir"/"$genome_file.fa.gz" ]]; then
	echo -e "[ "$(date)": ERROR: download failed, check password? ]"
	exit 1
fi

cat <<- _EOF_
	[ $(date): Downloading annotation ]
	$annot_url
_EOF_
curl $annot_url -b "$genome_dir"/cookies -c "$genome_dir"/cookies > $genome_dir/$annotation_file.gff3.gz

echo -e "[ "$(date)": Unzipping downloads ]"
gunzip $genome_dir/$genome_file.fa.gz
gunzip $genome_dir/$annotation_file.gff3.gz

# make cuffcomp gtf
echo -e "[ "$(date)": Making GTF file with cuffcompare ]"
cuffcompare -s $genome_dir/$genome_file.fa -CG -r $genome_dir/$annotation_file.gff3 \
	-o $genome_dir/$annotation_file.cuffcomp $genome_dir/$annotation_file.gff3

# stash new gtf
mv $genome_dir/$annotation_file.cuffcomp.combined.gtf $genome_dir/gtf_final.tmp

# remove cuffcomp intermediates
rm $genome_dir/*cuffcomp*
mv $genome_dir/gtf_final.tmp $genome_dir/$annotation_file.cuffcomp.gtf

# tidy up
rm "$genome_dir"/cookies

# log metadata
cat -t <<- _EOF_ > $genome_dir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cuffcomp version,$(cuffcompare -h 2>&1 | head -n 1)
_EOF_

exit 0
