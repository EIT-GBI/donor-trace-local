#!/usr/bin/env python3
"""Count recipient vs BAC alleles at the diagnostic SNPs.

  usage: allele_freq.py <bam> <ref.fa> <snps.tsv> <recipient|bac> <out_prefix>
                        [min_depth] [min_baseq]

Run this on a BAM aligned to recipient.fa alone, or to bac.fa alone - never on
the combined reference.

Why: in the combined reference the shared sequence is present TWICE. A read
from it fits both copies equally, so the aligner gives it MAPQ 0 and picks one
at random. Depth and alleles there are both meaningless.

With only one copy in the reference there is nothing to choose between.
Recipient-derived and BAC-derived reads stack up at the same position, and the
allele fraction is a direct measure of the mixture:

    BAC fraction ~ 0.0   only recipient sequence here
    BAC fraction ~ 0.5   one BAC copy present alongside the recipient locus
    BAC fraction ~ 1.0   the recipient locus itself has been replaced

The 'frame' argument picks which pair of allele columns to read from snps.tsv.
Pass 'recipient' for a recipient.fa alignment and 'bac' for a bac.fa alignment.
The two runs are independent measurements of the same thing.

Counting is allele-SPECIFIC: only the exact recipient base and the exact BAC
base are counted. Any third base is sequencing error and is reported separately
in the 'other' column, which doubles as an error floor measured at the very
positions being scored.
"""
import os
import subprocess
import sys


def read_snps(snps_tsv, frame):
    """-> list of (position, recipient_allele, bac_allele) in the chosen frame."""
    # snps.tsv names its columns with "rec" and "bac", so the frame picks both
    # the position column and which pair of allele columns to read.
    prefix = "rec" if frame == "recipient" else "bac"
    position_column = f"{prefix}_pos"
    recipient_column = f"rec_allele_{prefix}_frame"
    bac_column = f"bac_allele_{prefix}_frame"

    rows = []
    with open(snps_tsv) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            field = line.rstrip("\n").split("\t")
            rows.append((int(field[index[position_column]]),
                         field[index[recipient_column]],
                         field[index[bac_column]]))
    rows.sort()
    return rows


def count_bases(pileup_string):
    """mpileup base column -> {base: count}, plus 'REF' for . and ,

    The base column is not a plain list of letters. It also carries:
       ^X   start of a read, followed by one mapping-quality character
       $    end of a read
       +5ACGTA / -2AC   an indel, written as a length then that many letters
    Those extra letters are NOT aligned bases at this position, so the length
    has to be read and skipped. Counting them is a classic way to get wrong
    allele fractions.
    """
    counts = {}
    i = 0
    while i < len(pileup_string):
        character = pileup_string[i]
        if character == "^":
            i += 2                      # skip the marker and the quality char
        elif character == "$":
            i += 1
        elif character in "+-":
            digits = ""
            j = i + 1
            while j < len(pileup_string) and pileup_string[j].isdigit():
                digits += pileup_string[j]
                j += 1
            i = j + int(digits)         # skip the indel sequence itself
        elif character in "ACGTacgt":
            base = character.upper()
            counts[base] = counts.get(base, 0) + 1
            i += 1
        elif character in ".,":
            counts["REF"] = counts.get("REF", 0) + 1
            i += 1
        else:
            i += 1                      # * deletion placeholder, > < skips
    return counts


def run_mpileup(bam, ref_fa, contig, positions, min_baseq, out_prefix):
    """samtools mpileup restricted to the diagnostic positions."""
    regions_path = out_prefix + ".positions.bed"
    with open(regions_path, "w") as out:
        for position in positions:
            out.write(f"{contig}\t{position - 1}\t{position}\n")

    # -q 0 keeps low-MAPQ reads. On a single-copy reference they are not
    # ambiguous between copies, they are just reads in repetitive sequence,
    # and discarding them would throw away real data.
    result = subprocess.run(
        ["samtools", "mpileup", "-f", ref_fa, "-l", regions_path,
         "-q", "0", "-Q", str(min_baseq), "-d", "100000", "-B", bam],
        capture_output=True, text=True, check=True)
    os.remove(regions_path)
    return result.stdout


def main(bam, ref_fa, snps_tsv, frame, out_prefix, min_depth=10, min_baseq=10):
    if frame not in ("recipient", "bac"):
        sys.exit("frame must be 'recipient' or 'bac'")

    snps = read_snps(snps_tsv, frame)
    alleles_at = {position: (recipient, bac) for position, recipient, bac in snps}
    contig = frame

    pileup = run_mpileup(bam, ref_fa, contig,
                         [position for position, _, _ in snps],
                         min_baseq, out_prefix)

    rows = []
    for line in pileup.splitlines():
        field = line.split("\t")
        position = int(field[1])
        reference_base = field[2].upper()
        bases = field[4]
        if position not in alleles_at:
            continue

        recipient_allele, bac_allele = alleles_at[position]
        counts = count_bases(bases)

        # '.' and ',' mean "same as the reference base", so whichever of the
        # two alleles matches the reference is counted from 'REF'.
        if reference_base == recipient_allele:
            recipient_count = counts.get("REF", 0)
        else:
            recipient_count = counts.get(recipient_allele, 0)
        if reference_base == bac_allele:
            bac_count = counts.get("REF", 0)
        else:
            bac_count = counts.get(bac_allele, 0)

        total = recipient_count + bac_count
        if total < min_depth:
            continue
        other = sum(counts.values()) - total
        rows.append((position, recipient_count, bac_count, other,
                     bac_count / total))

    with open(out_prefix + ".tsv", "w") as out:
        out.write("pos\trecipient_count\tbac_count\tother\tbac_fraction\n")
        for position, recipient_count, bac_count, other, fraction in rows:
            out.write(f"{position}\t{recipient_count}\t{bac_count}"
                      f"\t{other}\t{fraction:.4f}\n")

    with open(out_prefix + ".bedgraph", "w") as out:
        name = os.path.basename(out_prefix)
        out.write(f'track type=bedGraph name="{name}" '
                  'viewLimits=0:1 autoScale=off visibility=full\n')
        for position, _, _, _, fraction in rows:
            out.write(f"{contig}\t{position - 1}\t{position}\t{fraction:.4f}\n")

    total_recipient = sum(row[1] for row in rows)
    total_bac = sum(row[2] for row in rows)
    total_other = sum(row[3] for row in rows)
    scored = total_recipient + total_bac
    overall = total_bac / scored if scored else float("nan")
    error_rate = total_other / (scored + total_other) if scored else 0.0
    print(f"[allele_freq] {os.path.basename(out_prefix)}"
          f"  frame={frame}  sites={len(rows)}/{len(snps)}"
          f"  bac_fraction={overall:.4f}  error_floor={error_rate:.4f}")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1], args[2], args[3], args[4],
         int(args[5]) if len(args) > 5 else 10,
         int(args[6]) if len(args) > 6 else 10)
