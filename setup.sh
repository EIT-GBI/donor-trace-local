#!/usr/bin/env bash
# One-off setup. Safe to re-run.
#
#   bash setup.sh
set -euo pipefail
cd "$(dirname "$0")"

echo "--- command line tools ---"
brew install minimap2 bwa samtools bcftools seqkit fastp

echo "--- python environment ---"
command -v uv > /dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
uv sync                              # pyproject.toml -> .venv, writes uv.lock

echo "--- bedGraphToBigWig ---"
# UCSC ships this as a binary and Homebrew does not carry it. Note the host is
# hgdownload.SOE.ucsc.edu - the hgdownload.cse name that most instructions
# still quote fails certificate validation.
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64)  BUILD=macOSX.arm64 ;;
  Darwin-x86_64) BUILD=macOSX.x86_64 ;;
  *)             BUILD=linux.x86_64 ;;
esac
mkdir -p tools/ucsc
curl -fL --progress-bar -o tools/ucsc/bedGraphToBigWig \
  "https://hgdownload.soe.ucsc.edu/admin/exe/$BUILD/bedGraphToBigWig" \
  && chmod +x tools/ucsc/bedGraphToBigWig \
  || echo "WARNING: download failed - the pipeline runs, but makes no tracks"

echo
echo "done. Next: edit config.txt and samplesheet_ont.tsv, then"
echo "  bash run.sh config.txt"
