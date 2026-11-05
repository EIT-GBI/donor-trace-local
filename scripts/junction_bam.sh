#!/usr/bin/env bash
# Every read with alignments on BOTH contigs, as a BAM you can open in IGV.
#
#   usage: junction_bam.sh <combined.bam> <out.bam>
#
# No filtering. These are candidates, not findings: most are reads from
# sequence that exists on both contigs, so the aligner had to split them.
# junction_reads.py is what separates those from real junctions.
#
# All alignments of a selected read are kept, supplementary included - that is
# what shows where the read actually goes.
set -euo pipefail

IN_BAM=$1
OUT_BAM=$2
NAMES="${OUT_BAM%.bam}.readnames.txt"

# -F 0x104 drops unmapped and secondary, KEEPS supplementary
samtools view -F 0x104 "$IN_BAM" recipient | cut -f1 | sort -u > "$NAMES.rec"
samtools view -F 0x104 "$IN_BAM" bac       | cut -f1 | sort -u > "$NAMES.bac"
comm -12 "$NAMES.rec" "$NAMES.bac" > "$NAMES"
rm -f "$NAMES.rec" "$NAMES.bac"

samtools view -b -N "$NAMES" "$IN_BAM" > "$OUT_BAM"
samtools index "$OUT_BAM"
