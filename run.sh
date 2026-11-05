#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# donor-trace-local
#
#   bash run.sh [config.txt]
#
# For each sample: is the donor sequence there, across how much of the region,
# and is it integrated?
#
# Three tiers, three different questions, three alignments. Each tier is
# independent of the others, so they can be compared rather than assumed.
#
#   TIER 1  combined.fa - recipient and BAC as TWO CONTIGS, supplementary
#           alignments kept. Align and LOOK AT THE ALIGNMENT: depth tracks and
#           IGV. Depth is trustworthy on BAC-only sequence (is a BAC there at
#           all?) and on recipient-only sequence (coverage, breakpoints), but
#           NOT on the shared region: the cargo is ~99% identical to its target,
#           so those reads get MAPQ 0 and are split between the contigs at
#           random. Structural surprises show up here.
#
#   TIER 2  recipient.fa and bac.fa separately. Each holds only ONE copy of the
#           shared sequence, so nothing is assigned - donor-derived and
#           recipient-derived reads stack at the same position and the ALLELE
#           FREQUENCY at the diagnostic SNPs is the mixture. This is how much
#           cargo is present and across exactly which interval. Run from both
#           sides; they share no alignment step, so agreeing is a real check.
#
#   TIER 3  JUNCTION READS, off the Tier 1 BAM. Tiers 1 and 2 measure how much
#           cargo is there; neither can say whether it is INTEGRATED, because an
#           episome at one copy per cell looks identical to an integrated copy.
#           Only one molecule containing both sequences distinguishes them, and
#           a read off it aligns in pieces that land on different contigs -
#           which is why combined.fa keeps supplementary alignments.
# ---------------------------------------------------------------------------
set -euo pipefail

CONFIG=${1:-config.txt}
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$CONFIG"

# Defaults, so a config that omits an optional setting fails with a sensible
# result rather than "unbound variable" from set -u.
: "${FORCE:=0}"
: "${TRIM:=0}"
: "${THREADS:=4}"
: "${MIN_DEPTH:=10}"
: "${MIN_BASEQ:=10}"
: "${MIN_ANCHOR:=500}"
: "${MIN_SYNTENY_LEN:=1000}"
: "${SMOOTH_WINDOW:=11}"
: "${RECIPIENT_LABEL:=recipient}"
: "${BAC_LABEL:=BAC}"
: "${BEDGRAPH_TO_BIGWIG:=tools/ucsc/bedGraphToBigWig}"

case "$PLATFORM" in
  ont|illumina) ;;
  *) echo "ERROR: PLATFORM must be 'ont' or 'illumina', got '$PLATFORM'"; exit 1 ;;
esac

# Activate the uv environment that setup.sh builds; the plots need matplotlib.
# A virtualenv stores its own absolute path, so rebuild it if this folder has
# been moved or renamed. The same test covers it simply not being there yet.
export PATH="$HOME/.local/bin:$PATH"          # where the uv installer puts uv
if ! grep -qs "VIRTUAL_ENV='$HERE/.venv'" "$HERE/.venv/bin/activate"; then
  rm -rf "$HERE/.venv"
  ( cd "$HERE" && uv sync --quiet )
fi
source "$HERE/.venv/bin/activate"

# A relative BEDGRAPH_TO_BIGWIG in the config means "inside the pipeline folder".
case "$BEDGRAPH_TO_BIGWIG" in
  /*) ;;
  *) BEDGRAPH_TO_BIGWIG="$HERE/$BEDGRAPH_TO_BIGWIG" ;;
esac

REF_DIR="$OUTDIR/reference"
BAM_DIR="$OUTDIR/bam"
BW_DIR="$OUTDIR/bigwig"
ALLELE_DIR="$OUTDIR/alleles"
PLOT_DIR="$OUTDIR/plots"
JUNCTION_DIR="$OUTDIR/junctions"
TRIM_DIR="$OUTDIR/trimmed"
LOG_DIR="$OUTDIR/logs"
mkdir -p "$REF_DIR" "$BAM_DIR" "$BW_DIR" "$ALLELE_DIR" "$PLOT_DIR" \
         "$JUNCTION_DIR" "$TRIM_DIR" "$LOG_DIR"

echo "=== donor-trace-local ==="
echo "    config   $CONFIG"
echo "    platform $PLATFORM"
echo "    outdir   $OUTDIR"

# --- the sample sheet -------------------------------------------------------
# The sheet is CSV. ONT sheets have no fastq2 column. Convert both layouts to
# one tab-separated five-column form (sample, fastq1, fastq2, recipient, bac)
# so the rest of the script has a single code path.
#
# `tr -d '\r'` because a CSV saved from Excel carries CRLF line endings, and a
# trailing \r silently becomes part of the last field - here a file path, which
# then does not exist for reasons invisible in the error message.
SHEET="$OUTDIR/samplesheet.normalised.tsv"
if [ "$PLATFORM" = "ont" ]; then
  tail -n +2 "$SAMPLESHEET" | tr -d '\r' \
    | awk -F',' -v OFS='\t' '{print $1, $2, "", $3, $4}' > "$SHEET"
else
  tail -n +2 "$SAMPLESHEET" | tr -d '\r' \
    | awk -F',' -v OFS='\t' '{print $1, $2, $3, $4, $5}' > "$SHEET"
fi

# --- 1. references, one set per BAC -----------------------------------------
echo
echo "--- references ---"
# Fields come out with `cut`, not `read a b c d e`: bash treats a tab as
# whitespace, so an empty fastq2 would merge two tabs and shift every field.
cut -f4,5 "$SHEET" | sort -u | while IFS= read -r LINE; do
  RECIPIENT_FA=$(printf '%s' "$LINE" | cut -f1)
  BAC_FA=$(printf '%s' "$LINE" | cut -f2)
  [ -z "$BAC_FA" ] && continue
  PAIR_DIR="$REF_DIR/$(basename "$BAC_FA" | sed 's/\.[^.]*$//')"

  if [ -f "$PAIR_DIR/snps.tsv" ] && [ "$FORCE" != "1" ]; then
    echo "[skip] $(basename "$PAIR_DIR") already built"
    continue
  fi
  bash "$HERE/scripts/prepare_reference.sh" \
       "$RECIPIENT_FA" "$BAC_FA" "$PAIR_DIR" "$PLATFORM"
  python "$HERE/scripts/find_synteny.py" \
       "$PAIR_DIR/recipient.fa" "$PAIR_DIR/bac.fa" "$PAIR_DIR" "$MIN_SYNTENY_LEN"
done

# --- 2. per sample ----------------------------------------------------------
echo
echo "--- samples ---"
while IFS= read -r LINE; do
  SAMPLE=$(printf '%s' "$LINE" | cut -f1)
  FASTQ1=$(printf '%s' "$LINE" | cut -f2)
  FASTQ2=$(printf '%s' "$LINE" | cut -f3)
  BAC_FA=$(printf '%s' "$LINE" | cut -f5)
  [ -z "$SAMPLE" ] && continue
  PAIR=$(basename "$BAC_FA" | sed 's/\.[^.]*$//')
  PAIR_DIR="$REF_DIR/$PAIR"

  [ -f "$PAIR_DIR/snps.tsv" ] || { echo "ERROR: no reference for '$PAIR' ($SAMPLE)" >&2; exit 1; }
  [ -f "$FASTQ1" ] || { echo "ERROR: fastq1 not found for $SAMPLE: $FASTQ1" >&2; exit 1; }

  # the junction table is the LAST thing written, so it is what "done" means
  if [ -f "$JUNCTION_DIR/$SAMPLE.junctions.tsv" ] && [ "$FORCE" != "1" ]; then
    echo "[skip] $SAMPLE already done"
    continue
  fi
  echo "[$(date +%H:%M:%S)] $SAMPLE  (BAC: $PAIR)"

  if [ "$TRIM" = "1" ]; then
    if [ "$PLATFORM" = "illumina" ]; then
      fastp -w 4 -i "$FASTQ1" -I "$FASTQ2" \
            -o "$TRIM_DIR/$SAMPLE.1.fq.gz" -O "$TRIM_DIR/$SAMPLE.2.fq.gz" \
            --detect_adapter_for_pe -l 40 \
            -j "$LOG_DIR/$SAMPLE.fastp.json" -h "$LOG_DIR/$SAMPLE.fastp.html" \
            2> "$LOG_DIR/$SAMPLE.fastp.log"
      FASTQ1="$TRIM_DIR/$SAMPLE.1.fq.gz"
      FASTQ2="$TRIM_DIR/$SAMPLE.2.fq.gz"
    else
      fastp -i "$FASTQ1" -o "$TRIM_DIR/$SAMPLE.fq.gz" -l 200 -Q \
            -j "$LOG_DIR/$SAMPLE.fastp.json" -h "$LOG_DIR/$SAMPLE.fastp.html" \
            2> "$LOG_DIR/$SAMPLE.fastp.log"
      FASTQ1="$TRIM_DIR/$SAMPLE.fq.gz"
    fi
  fi

  # combined keeps supplementary alignments: without them a read spanning
  # recipient and BAC is reported as one piece and the junction evidence is
  # gone. recipient and bac drop them - there each read should have one home.
  for REF in combined recipient bac; do
    echo "    align -> $REF"
    if [ "$REF" = "combined" ]; then KEEP=keep; else KEEP=; fi
    bash "$HERE/scripts/align.sh" \
         "$PAIR_DIR/$REF.fa" "$PLATFORM" "$THREADS" \
         "$BAM_DIR/$SAMPLE.$REF.bam" "$FASTQ1" "$FASTQ2" "$KEEP" \
         2> "$LOG_DIR/$SAMPLE.$REF.align.log"
    bash "$HERE/scripts/make_bigwig.sh" \
         "$BAM_DIR/$SAMPLE.$REF.bam" "$BW_DIR/$SAMPLE.$REF.bw" \
         "$BEDGRAPH_TO_BIGWIG"
  done

  # allele fractions, once from each side
  for FRAME in recipient bac; do
    python "$HERE/scripts/allele_freq.py" \
         "$BAM_DIR/$SAMPLE.$FRAME.bam" "$PAIR_DIR/$FRAME.fa" \
         "$PAIR_DIR/snps.tsv" "$FRAME" "$ALLELE_DIR/$SAMPLE.$FRAME" \
         "$MIN_DEPTH" "$MIN_BASEQ"
    python "$HERE/scripts/plot_alleles.py" \
         "$ALLELE_DIR/$SAMPLE.$FRAME.tsv" "$FRAME" \
         "$SAMPLE  -  aligned to $FRAME only" \
         "$PLOT_DIR/$SAMPLE.$FRAME.png" \
         "$SMOOTH_WINDOW" "$RECIPIENT_LABEL" "$BAC_LABEL"
  done

  # junctions. Runs after the allele step because it needs
  # <sample>.recipient.tsv to know which SNPs are still informative here.
  bash "$HERE/scripts/junction_bam.sh" \
       "$BAM_DIR/$SAMPLE.combined.bam" "$JUNCTION_DIR/$SAMPLE.junction.bam"
  python "$HERE/scripts/junction_reads.py" \
       "$JUNCTION_DIR/$SAMPLE.junction.bam" "$PAIR_DIR/homology_all.bed" \
       "$PAIR_DIR/snps.tsv" "$ALLELE_DIR/$SAMPLE.recipient.tsv" \
       "$JUNCTION_DIR/$SAMPLE" "$MIN_ANCHOR"
done < "$SHEET"

# --- 3. summary -------------------------------------------------------------
echo
echo "--- summary ---"
python "$HERE/scripts/summarise.py" "$ALLELE_DIR" "$OUTDIR/summary.tsv"

echo
echo "--- junctions ---"
python "$HERE/scripts/summarise_junctions.py" \
       "$JUNCTION_DIR" "$OUTDIR/summary_junctions.tsv"
