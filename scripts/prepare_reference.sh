#!/usr/bin/env bash
# Build one reference set for a (recipient, BAC) pair.
#
#   usage: prepare_reference.sh <recipient.fa> <bac.fa> <outdir> <platform>
#
# Produces three FASTAs, all indexed:
#   recipient.fa   the recipient genome alone, contig 'recipient'
#   bac.fa         the BAC alone,              contig 'bac'
#   combined.fa    both, as TWO CONTIGS: 'recipient' and 'bac'
#
# combined.fa is a plain concatenation - two contigs, no spacer, no offsets. A
# BAC position is just a position on 'bac', so synteny.bed, homology_all.bed
# and snps.bed load against any of the three references without conversion.
#
# It is used for the depth tracks and for junction detection, and for nothing
# else. Do not read allele fractions off it: it holds BOTH copies of the
# shared sequence, so a syntenic read fits each equally well, gets MAPQ 0 and
# is assigned at random. That is what recipient.fa and bac.fa are for - each
# holds only ONE copy, so there is nothing to assign.
set -euo pipefail

RECIPIENT_FA=$1
BAC_FA=$2
OUTDIR=$3
PLATFORM=$4

mkdir -p "$OUTDIR"

# Fixed contig names, so every later script can rely on them.
seqkit replace -p '.*' -r 'recipient' "$RECIPIENT_FA" > "$OUTDIR/recipient.fa"
seqkit replace -p '.*' -r 'bac'       "$BAC_FA"       > "$OUTDIR/bac.fa"
cat "$OUTDIR/recipient.fa" "$OUTDIR/bac.fa" > "$OUTDIR/combined.fa"

for FA in recipient bac combined; do
  samtools faidx "$OUTDIR/$FA.fa"
  if [ "$PLATFORM" = "illumina" ]; then
    bwa index "$OUTDIR/$FA.fa" 2> "$OUTDIR/$FA.bwa-index.log"
  fi
done

echo "[prepare_reference] $OUTDIR"
cut -f1,2 "$OUTDIR/combined.fa.fai" | sed 's/^/    /'
