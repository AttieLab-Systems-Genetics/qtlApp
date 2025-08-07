#!/bin/bash

# Directory containing the peak summary files
PEAK_DIR="/mnt/rdrive/mkeller3/General/main_directory/annotated_peak_summaries"

# Check if the directory exists
if [ ! -d "$PEAK_DIR" ]; then
  echo "Error: Directory not found: $PEAK_DIR"
  exit 1
fi

echo "Analyzing files in: $PEAK_DIR"
echo "This may take a few moments..."

# Calculate total lines and unique LODs in a more robust, memory-efficient way
# by streaming the data instead of loading it into a variable.
{
  read total_lines
  read unique_lods
} < <(find "$PEAK_DIR" -type f -name "*.csv" -exec tail -q -n +2 {} + | tee \
    >(wc -l) \
    >(cut -d, -f3 | sort -u | wc -l) >/dev/null)

echo "========================================"
echo "      Peak Summary Analysis"
echo "========================================"
echo "Total data rows (peaks):    ${total_lines}"
echo "Unique qtl_lod values:      ${unique_lods}"
echo "========================================"

