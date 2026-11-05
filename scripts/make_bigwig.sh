#!/usr/bin/env bash
# Depth track for IGV.
#
#   usage: make_bigwig.sh <in.bam> <out.bw> <bedGraphToBigWig>
#
# Skips quietly if the UCSC binary is missing - a missing track should not
# stop the pipeline.
set -euo pipefail

IN_BAM=$1
OUT_BW=$2
BG2BW=$3

if [ ! -x "$BG2BW" ]; then
  echo "[make_bigwig] no bedGraphToBigWig, skipping $(basename "$OUT_BW")"
  exit 0
fi

BEDGRAPH="${OUT_BW%.bw}.bedgraph"
CHROM_SIZES="${OUT_BW%.bw}.chrom_sizes"

samtools idxstats "$IN_BAM" | grep -v '^\*' | cut -f1,2 > "$CHROM_SIZES"
samtools depth -a -d 0 "$IN_BAM" \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > "$BEDGRAPH"

"$BG2BW" "$BEDGRAPH" "$CHROM_SIZES" "$OUT_BW"
rm -f "$BEDGRAPH" "$CHROM_SIZES"
