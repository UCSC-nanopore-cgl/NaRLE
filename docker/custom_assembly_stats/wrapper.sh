#!/usr/bin/env bash
set -e

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

# gather time and memory statistics, save everything to log
echo "/usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' python3 /opt/nanopore_assembly_and_polishing_assessment/align_and_summarize_contigs.py $@\n" > /data/custom_assembly_stats.log
eval "/usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' python3 /opt/nanopore_assembly_and_polishing_assessment/align_and_summarize_contigs.py $@" 2>&1 | tee -a /data/custom_assembly_stats.log

