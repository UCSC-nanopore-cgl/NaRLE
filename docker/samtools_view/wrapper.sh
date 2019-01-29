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
if [[ $1 == view ]]; then
	echo "samtools_view docker invocation must not specify 'view' subprogram"
	exit 1
fi

# gather time and memory statistics, save everything to 'rle.log'
echo "/usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' samtools view $@ >/data/samtools_view.out\n" > /data/samtools_view.log
eval "/usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' samtools view $@ >/data/samtools_view.out" 2>&1 | tee -a /data/samtools_view.log

