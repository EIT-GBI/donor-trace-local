#!/usr/bin/env python3
"""One table of the junction results across all samples.

    usage: summarise_junctions.py <junction_dir> <out.tsv>

Reads every *.junctions.tsv and reports, per sample:

  candidates    reads with alignments on BOTH contigs. Mostly meaningless on
                its own - the cargo matches both, so reads from it split
  junctions     candidates that survive: a BAC-only piece AND a piece that is
                really the recipient genome. THIS is the integration count
  why_dropped   what disqualified the rest, so a zero can be read rather than
                taken on trust

A non-zero junction count is not proof on its own. Look at the reads in
junction_bam/, check the positions cluster, and check MAPQ.
"""
import collections
import glob
import os
import sys


def main(junction_dir, out_tsv):
    lines = []
    for path in sorted(glob.glob(os.path.join(junction_dir, "*.junctions.tsv"))):
        sample = os.path.basename(path)[:-len(".junctions.tsv")]
        reads, verdicts = {}, collections.Counter()
        with open(path) as handle:
            index = {n: i for i, n in
                     enumerate(handle.readline().rstrip("\n").split("\t"))}
            for line in handle:
                f = line.rstrip("\n").split("\t")
                reads[f[index["read"]]] = f[index["junction"]]
                verdicts[f[index["verdict"]]] += 1

        junctions = sum(1 for v in reads.values() if v == "yes")
        why = ", ".join(f"{k}={v}" for k, v in verdicts.most_common())
        lines.append((sample, len(reads), junctions, why or "-"))

    with open(out_tsv, "w") as out:
        out.write("sample\tcandidates\tjunction_reads\tpiece_verdicts\n")
        for row in lines:
            out.write("\t".join(str(x) for x in row) + "\n")

    width = max([len(r[0]) for r in lines] + [6])
    print(f"\n{'sample':<{width}}  {'candidates':>10} {'junctions':>9}   verdicts")
    for sample, candidates, junctions, why in lines:
        print(f"{sample:<{width}}  {candidates:>10} {junctions:>9}   {why}")
    total = sum(r[2] for r in lines)
    print(f"\n{'NO INTEGRATION FOUND' if total == 0 else f'{total} junction reads - inspect them'}")
    print(f"written: {out_tsv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
