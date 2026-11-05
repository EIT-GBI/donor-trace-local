#!/usr/bin/env python3
"""One table across all samples, and a look for conversion tracts.

    usage: summarise.py <allele_dir> <out.tsv>

Reads every *.tsv written by allele_freq.py and reports, per sample and frame:

  sites         diagnostic SNPs with enough depth to score
  bac_fraction  BAC allele fraction over the whole homologous region
  error         third-allele rate - neither expected allele, so sequencing
                error, measured at the same positions. The bac_fraction of a
                sample with no BAC should land near this number, not at zero
  verdict       from bac_fraction (see below)
  tract         longest unbroken run of sites sitting at recipient sequence,
                inside a sample that otherwise carries BAC. That is a block of
                cargo replaced by recipient sequence

A tract is only reported when the sample carries BAC elsewhere, because in a
sample with no BAC at all every site is "recipient" and the whole region would
be reported as one meaningless tract.
"""
import glob
import os
import sys

ABSENT_BELOW = 0.02        # at or under this, call it absent
REPLACED_ABOVE = 0.80      # at or over this, the recipient locus is gone
TRACT_BELOW = 0.05         # a site this low is recipient sequence
MIN_TRACT_SITES = 20       # shorter runs are noise, not a tract
MAX_TRACT_GAP = 3          # isolated high sites inside a tract to step over


def read_table(path):
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            field = line.rstrip("\n").split("\t")
            rows.append((int(field[index["pos"]]),
                         int(field[index["recipient_count"]]),
                         int(field[index["bac_count"]]),
                         int(field[index["other"]]),
                         float(field[index["bac_fraction"]])))
    return rows


def longest_recipient_run(rows):
    """-> (number_of_sites, start_pos, end_pos) for the longest low run.

    Short interruptions are stepped over. A couple of isolated sites can read
    high inside a genuine tract - a mis-called SNP, or a position sitting on an
    indel where the two alignments disagree. Ending the run there would report
    one 38 kb tract as two 19 kb ones. A run only ends after MAX_TRACT_GAP
    consecutive high sites.
    """
    best_count, best_start, best_end = 0, 0, 0
    count, start, end, gap = 0, None, None, 0

    for position, _, _, _, fraction in rows:
        if fraction < TRACT_BELOW:
            if start is None:
                start = position
            count += 1
            end = position
            gap = 0
        elif start is not None:
            gap += 1
            if gap > MAX_TRACT_GAP:
                if count > best_count:
                    best_count, best_start, best_end = count, start, end
                count, start, end, gap = 0, None, None, 0

    if start is not None and count > best_count:
        best_count, best_start, best_end = count, start, end
    return best_count, best_start, best_end


def main(allele_dir, out_tsv):
    lines = []
    for path in sorted(glob.glob(os.path.join(allele_dir, "*.tsv"))):
        name = os.path.basename(path)[:-4]
        sample, _, frame = name.rpartition(".")
        rows = read_table(path)
        if not rows:
            continue

        recipient_total = sum(row[1] for row in rows)
        bac_total = sum(row[2] for row in rows)
        other_total = sum(row[3] for row in rows)
        scored = recipient_total + bac_total
        bac_fraction = bac_total / scored
        error = other_total / (scored + other_total)

        if bac_fraction <= ABSENT_BELOW:
            verdict = "no BAC cargo"
        elif bac_fraction >= REPLACED_ABOVE:
            verdict = "recipient locus replaced"
        else:
            verdict = "BAC cargo present"

        tract = "-"
        if bac_fraction > ABSENT_BELOW:
            count, start, end = longest_recipient_run(rows)
            if count >= MIN_TRACT_SITES:
                tract = f"{start}-{end} ({(end - start) / 1000:.1f} kb, {count} sites)"

        lines.append((sample, frame, len(rows), bac_fraction, error,
                      verdict, tract))

    with open(out_tsv, "w") as out:
        out.write("sample\tframe\tsites\tbac_fraction\terror\tverdict\ttract\n")
        for row in lines:
            out.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]:.4f}"
                      f"\t{row[4]:.4f}\t{row[5]}\t{row[6]}\n")

    width = max([len(row[0]) for row in lines] + [6])
    print(f"\n{'sample':<{width}}  {'frame':<9} {'sites':>6} {'bacFrac':>8}"
          f" {'error':>7}  {'verdict':<24} tract")
    for row in lines:
        print(f"{row[0]:<{width}}  {row[1]:<9} {row[2]:>6} {row[3]:>8.4f}"
              f" {row[4]:>7.4f}  {row[5]:<24} {row[6]}")
    print(f"\nwritten: {out_tsv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
