#!/usr/bin/env python3
"""Does any read run from the recipient genome into the BAC?

  usage: junction_reads.py <junction.bam> <homology_all.bed> <snps.tsv>
                           <allele_freq.tsv> <out_prefix> [min_anchor]

Takes the BAM junction_bam.sh produces - the candidates, already reduced to
reads with alignments on both contigs, with ALL their alignments. Running this
on the full combined BAM gives the same answer but loads every read sequence
in the genome to look at a few hundred: 227,704 records instead of 308 on this
data.

Pass homology_all.bed, NOT synteny.bed. synteny.bed drops blocks below
MIN_SYNTENY_LEN because they are too small to be cargo - but a 731 bp shared
element is still quite capable of making a read align to both sequences, and
every read crossing it then looks like a junction.

Integration means one molecule contains both. A read off such a molecule
cannot align in one piece - it aligns as a primary plus supplementary
alignments, and on a two-contig reference those pieces land on different
contigs.

Most reads with alignments on both contigs mean nothing, though. The cargo is
~99% identical to its target, so a read from it genuinely matches both and
the aligner simply had to choose. Two tests separate those out.

  BAC side    must lie OUTSIDE every syntenic block. Vector has no recipient
              counterpart, so an alignment there is unambiguous on its own.

  recipient   outside the syntenic blocks: unambiguous, accepted.
  side        inside them: ask the SNPs rather than the coordinates -
                  recipient alleles -> really the recipient genome
                  donor alleles     -> cargo, placed here by the aligner

Only sites where BOTH alleles actually occur in this sample are used, which is
what allele_freq.tsv is for. A listed SNP inside a conversion tract is not
informative - both molecules read recipient there - and scoring it anyway
makes BAC reads look like recipient reads. That is not hypothetical; it
invents junctions.

That table must come from the RECIPIENT-ONLY alignment, which is what the
pipeline already produces. Measured on this two-contig BAM the numbers are
wrong: the BAC is in the reference, so cargo reads go to the bac contig and
the recipient side looks almost pure.
"""
import collections
import os
import re
import subprocess
import sys

MIN_ANCHOR = 500          # a shorter piece is not convincing on ONT
MIN_SITES = 5             # fewer diagnostic SNPs than this is not a genotype
CALL = 0.8                # allele fraction needed to call a piece either way
OVERLAP_FRACTION = 0.5    # above this much overlap with synteny: ambiguous
MIN_INFORMATIVE = 0.1     # a site needs both alleles present in the sample


def read_homology(path):
    """-> (bac_blocks, recipient_blocks), each a list of (start, end)."""
    bac_blocks, recipient_blocks = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        recipient_blocks.append((int(f[1]), int(f[2])))
        bac_blocks.append((int(f[3]), int(f[4])))
    return bac_blocks, recipient_blocks


def read_snps(path):
    """recipient position (1-based) -> (recipient_base, bac_base)."""
    sites = {}
    with open(path) as handle:
        index = {n: i for i, n in
                 enumerate(handle.readline().rstrip("\n").split("\t"))}
        for line in handle:
            f = line.rstrip("\n").split("\t")
            sites[int(f[index["rec_pos"]])] = (f[index["rec_allele_rec_frame"]],
                                               f[index["bac_allele_rec_frame"]])
    return sites


def informative_sites(path):
    """Positions where both alleles occur in this sample."""
    sites = set()
    with open(path) as handle:
        index = {n: i for i, n in
                 enumerate(handle.readline().rstrip("\n").split("\t"))}
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if MIN_INFORMATIVE <= float(f[index["bac_fraction"]]) <= 1 - MIN_INFORMATIVE:
                sites.add(int(f[index["pos"]]))
    return sites


def overlap(start, end, blocks):
    return sum(max(0, min(end, b) - max(start, a)) for a, b in blocks)


def bases_at(pos, cigar, seq, wanted):
    """Bases this alignment puts at each reference position in `wanted`."""
    found = {}
    ref, query = pos, 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in "M=X":
            for offset in range(length):
                if ref + offset in wanted:
                    found[ref + offset] = seq[query + offset].upper()
            ref += length
            query += length
        elif op in "DN":
            ref += length
        elif op in "IS":
            query += length
    return found


def main(bam, homology_bed, snps_tsv, allele_tsv, out_prefix,
         min_anchor=MIN_ANCHOR):
    bac_blocks, recipient_blocks = read_homology(homology_bed)
    snps = read_snps(snps_tsv)
    usable = informative_sites(allele_tsv) & set(snps)

    process = subprocess.run(["samtools", "view", "-F", "0x104", bam],
                             capture_output=True, text=True, check=True)
    reads = collections.defaultdict(list)
    for line in process.stdout.splitlines():
        f = line.split("\t")
        if f[5] == "*":
            continue
        span = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5])
                   if op in "MDN=X")
        reads[f[0]].append({"contig": f[2], "start": int(f[3]),
                            "end": int(f[3]) + span, "mapq": int(f[4]),
                            "cigar": f[5], "seq": f[9]})

    rows, junctions = [], []
    for name, alignments in reads.items():
        if not {"recipient", "bac"} <= {a["contig"] for a in alignments}:
            continue

        verdicts = []
        for a in alignments:
            length = a["end"] - a["start"]
            blocks = bac_blocks if a["contig"] == "bac" else recipient_blocks
            shared = overlap(a["start"], a["end"], blocks) / length
            recipient = donor = 0
            if length < min_anchor:
                verdict = "TOO_SHORT"
            elif shared <= OVERLAP_FRACTION:
                verdict = "BAC_ONLY" if a["contig"] == "bac" else "RECIPIENT_ONLY"
            elif a["contig"] == "bac":
                verdict = "SYNTENIC"          # cargo; no way to tell from here
            else:
                # a hard-clipped record can carry SEQ "*", with no bases to read
                wanted = set() if a["seq"] == "*" else {
                    p for p in usable if a["start"] <= p < a["end"]}
                for position, base in bases_at(a["start"], a["cigar"],
                                               a["seq"], wanted).items():
                    recipient_base, bac_base = snps[position]
                    if base == recipient_base:
                        recipient += 1
                    elif base == bac_base:
                        donor += 1
                scored = recipient + donor
                if scored < MIN_SITES:
                    verdict = "TOO_FEW_SNPS"
                elif recipient / scored >= CALL:
                    verdict = "RECIPIENT_BY_SNP"
                elif donor / scored >= CALL:
                    verdict = "CARGO_BY_SNP"
                else:
                    verdict = "MIXED"
            verdicts.append((a, verdict, recipient, donor, shared))

        labels = {v for _, v, _, _, _ in verdicts}
        is_junction = ("BAC_ONLY" in labels
                       and bool(labels & {"RECIPIENT_ONLY", "RECIPIENT_BY_SNP"}))
        if is_junction:
            junctions.append(name)
        for a, verdict, recipient, donor, shared in sorted(
                verdicts, key=lambda v: (v[0]["contig"], v[0]["start"])):
            rows.append((name, a["contig"], a["start"], a["end"],
                         a["end"] - a["start"], a["mapq"], f"{shared:.3f}",
                         recipient, donor, verdict,
                         "yes" if is_junction else "no"))

    with open(out_prefix + ".junctions.tsv", "w") as out:
        out.write("read\tcontig\tstart\tend\tlength\tmapq\tsynteny_overlap"
                  "\trecipient_snps\tdonor_snps\tverdict\tjunction\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    candidates = len({r[0] for r in rows})
    print(f"[junction_reads] {os.path.basename(out_prefix)}"
          f"  candidates={candidates}  junction_reads={len(junctions)}"
          f"  informative_snps={len(usable)}/{len(snps)}")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], a[3], a[4],
         int(a[5]) if len(a) > 5 else MIN_ANCHOR)
