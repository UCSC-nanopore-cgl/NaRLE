#!/usr/bin/env bash
set -e

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

# Determine whether we should use fasta or fastq encoding
if [ ! -f "$1" ]; then
	echo "samtools_sort docker invocation must specify file to be sorted as first argument"
	exit 1
elif [[ $1 == *sam ]]; then
	export INPUT_PIPE="cat $1"
elif [[ $1 == *bam ]]; then
	export INPUT_PIPE="samtools view -h $1"
else
	echo "samtools_sort docker invocation must specify .sam or .bam as first argument"
	exit 1
fi

# gather time and memory statistics, save everything to 'rle.log'
echo -e "$INPUT_PIPE | /usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' samtools sort ${@:2} | samtools view -hb >/data/samtools_sort.bam\n" > /data/samtools_sort.log
eval "$INPUT_PIPE | /usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' samtools sort ${@:2} | samtools view -hb >/data/samtools_sort.bam" 2>&1 | tee -a /data/samtools_sort.log

