#!/usr/bin/env bash
set -e

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

# gather time and memory statistics, save everything to 'rtg.log'
echo "/usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' java $JAVA_OPTS -jar /opt/rtg_tools/rtg-tools.jar $@\n" > /data/rtg.log
eval "/usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' java $JAVA_OPTS -jar /opt/rtg_tools/rtg-tools.jar $@" 2>&1 | tee -a /data/rtg.log
