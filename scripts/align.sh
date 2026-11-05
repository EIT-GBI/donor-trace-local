#!/usr/bin/env bash
# Align one sample to one reference.
#
#   usage: align.sh <ref.fa> <platform> <threads> <out.bam> <fastq1> <fastq2>
#                   [keep_supplementary]
#
# By default supplementary (0x800) and secondary (0x100) alignments are
# dropped, so every read contributes to at most one position. Unmapped reads
# (0x4) go too: 0x4 + 0x100 + 0x800 = 0x904.
#
# Pass a 7th argument of "keep" to keep supplementary alignments (0x104 only).
# Junction detection needs them - a read spanning recipient and BAC is
# reported as one primary plus supplementary alignments, and dropping those
# deletes exactly the evidence being looked for.
set -euo pipefail

REF=$1
PLATFORM=$2
THREADS=$3
OUT_BAM=$4
FASTQ1=$5
FASTQ2=${6:-}
KEEP_SUPPLEMENTARY=${7:-}

if [ "$KEEP_SUPPLEMENTARY" = "keep" ]; then
  EXCLUDE=0x104
else
  EXCLUDE=0x904
fi

if [ "$PLATFORM" = "ont" ]; then
  minimap2 -t "$THREADS" -ax map-ont --secondary=no "$REF" "$FASTQ1"
else
  bwa mem -t "$THREADS" "$REF" "$FASTQ1" "$FASTQ2"
fi \
  | samtools view -b -F "$EXCLUDE" - \
  | samtools sort -@ 2 -o "$OUT_BAM" -

samtools index "$OUT_BAM"
