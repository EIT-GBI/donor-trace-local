#!/usr/bin/env python3
"""Align the BAC to the recipient genome and write out what they share.

    usage: find_synteny.py <recipient.fa> <bac.fa> <outdir> [min_block_len]

Two outputs:

  synteny.bed   the blocks where BAC and recipient are homologous, in BOTH
                coordinate systems, keeping only blocks of at least
                min_block_len. This defines the cargo.

  homology_all.bed
                EVERY scrap of shared sequence, found with a deliberately
                over-sensitive search and no length filter. Junction detection
                must use this one. A short shared element - an IS copy, a
                repeat - is far too small to be "the cargo" but is perfectly
                capable of making a read align to both sequences, and a
                junction test that has not excluded it reports that read as a
                junction.

                Erring towards too many blocks is the safe direction here: an
                extra block can only make the junction test more conservative,
                while a missing one invents junctions. A real case: a 731 bp
                element at 81% identity, shared between the miniREXER vector
                and the chromosome, is found by `-x asm20` but dropped once
                `-c` is added. That single missing block turned 103 reads into
                false junctions.

  snps.tsv      every single-base difference inside those blocks. These are the
                diagnostic positions: a read covering one of them can be told
                apart as recipient-derived or BAC-derived.

ABOUT STRAND - read this before using snps.tsv
-----------------------------------------------
A BAC often aligns to the recipient on the minus strand. When it does, a base
written in BAC orientation is the complement of the same base written in
recipient orientation. Comparing them directly silently destroys the data:
A<->T and C<->G differences look like "no difference" at all, and every
transition reads as a transversion.

So this file never asks you to complement anything. It writes all four alleles
out explicitly, two per coordinate frame:

    rec_allele_rec_frame   what the RECIPIENT has, written on the recipient
    bac_allele_rec_frame   what the BAC has,       written on the recipient
    rec_allele_bac_frame   what the RECIPIENT has, written on the BAC
    bac_allele_bac_frame   what the BAC has,       written on the BAC

Scoring a recipient alignment uses the two *_rec_frame columns. Scoring a BAC
alignment uses the two *_bac_frame columns. Neither ever needs the other pair.
"""
import os
import re
import subprocess
import sys

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def run_minimap2(recipient_fa, bac_fa, paf_path, options):
    with open(paf_path, "w") as out:
        subprocess.run(["minimap2"] + options + [recipient_fa, bac_fa],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)


# base-level alignment, for the cs string the SNP caller needs
SYNTENY_OPTIONS = ["-x", "asm20", "-c", "--cs"]
# no -c, and secondary hits kept: finds short, diverged shared elements that
# the call above silently discards
HOMOLOGY_OPTIONS = ["-x", "map-ont", "-N", "50", "-p", "0.1"]


def parse_paf_line(line):
    f = line.rstrip("\n").split("\t")
    cs = None
    for field in f[12:]:
        if field.startswith("cs:Z:"):
            cs = field[5:]
    return {
        "bac_start": int(f[2]), "bac_end": int(f[3]),     # 0-based, forward BAC
        "strand": f[4],
        "rec_start": int(f[7]), "rec_end": int(f[8]),     # 0-based, forward recipient
        "matches": int(f[9]), "block_len": int(f[10]),
        "cs": cs,
    }


# The column order of snps.tsv. snps_in_alignment() yields tuples in exactly
# this order, so the file is written without naming the fields again.
SNP_COLUMNS = ["rec_pos", "bac_pos", "strand",
               "rec_allele_rec_frame", "bac_allele_rec_frame",
               "rec_allele_bac_frame", "bac_allele_bac_frame"]


def snps_in_alignment(aln):
    """Walk the cs string and yield one tuple per single-base difference.

    minimap2 reports target (recipient) coordinates on the forward strand
    always, but walks the QUERY (BAC) in alignment orientation. So on a minus
    alignment the cs query bases come out complemented, and the query position
    counts backwards. Both are undone here, once, in this function.
    """
    minus = aln["strand"] == "-"
    bac_walk = aln["bac_start"]      # position along the ALIGNED bac orientation
    rec_pos = aln["rec_start"]       # position along the forward recipient

    for op, value in re.findall(r"([:*=+-])([0-9]+|[A-Za-z]+)", aln["cs"]):
        if op == ":":                       # identical run, length given
            n = int(value)
            bac_walk += n
            rec_pos += n
        elif op == "*":                     # substitution: value = ref + query
            rec_base = value[0].upper()     # recipient, forward recipient strand
            bac_base = value[1].upper()     # BAC, but in ALIGNMENT orientation

            if minus:
                # count backwards from the far end to get the forward BAC position
                bac_pos = aln["bac_end"] - 1 - (bac_walk - aln["bac_start"])
            else:
                bac_pos = bac_walk

            # In alignment orientation both bases sit on the recipient's strand.
            # Flipping to BAC orientation complements both, or neither.
            if minus:
                rec_in_bac_frame = COMPLEMENT[rec_base]
                bac_in_bac_frame = COMPLEMENT[bac_base]
            else:
                rec_in_bac_frame = rec_base
                bac_in_bac_frame = bac_base

            yield (rec_pos + 1, bac_pos + 1, aln["strand"],   # 1-based
                   rec_base, bac_base,
                   rec_in_bac_frame, bac_in_bac_frame)
            bac_walk += 1
            rec_pos += 1
        elif op == "+":                     # insertion in the BAC
            bac_walk += len(value)
        elif op == "-":                     # deletion from the BAC
            rec_pos += len(value)


def main(recipient_fa, bac_fa, outdir, min_block_len=1000):
    os.makedirs(outdir, exist_ok=True)
    paf_path = os.path.join(outdir, "bac_vs_recipient.paf")
    run_minimap2(recipient_fa, bac_fa, paf_path, SYNTENY_OPTIONS)
    with_cs = [parse_paf_line(line) for line in open(paf_path)]
    alignments = [a for a in with_cs
                  if a["cs"] and a["block_len"] >= min_block_len]

    homology_paf = os.path.join(outdir, "homology_all.paf")
    run_minimap2(recipient_fa, bac_fa, homology_paf, HOMOLOGY_OPTIONS)
    every_block = [parse_paf_line(line) for line in open(homology_paf)]

    def write_blocks(path, blocks):
        with open(path, "w") as out:
            out.write("#recipient\trec_start\trec_end\tbac_start\tbac_end"
                      "\tstrand\tidentity\n")
            for aln in sorted(blocks, key=lambda a: a["rec_start"]):
                identity = aln["matches"] / aln["block_len"]
                out.write(f"recipient\t{aln['rec_start']}\t{aln['rec_end']}"
                          f"\t{aln['bac_start']}\t{aln['bac_end']}"
                          f"\t{aln['strand']}\t{identity:.4f}\n")

    write_blocks(os.path.join(outdir, "synteny.bed"), alignments)
    write_blocks(os.path.join(outdir, "homology_all.bed"), every_block)

    snps = []
    for aln in alignments:
        snps.extend(snps_in_alignment(aln))
    snps.sort()                                  # by rec_pos, the first field

    with open(os.path.join(outdir, "snps.tsv"), "w") as out:
        out.write("\t".join(SNP_COLUMNS) + "\n")
        for snp in snps:
            out.write("\t".join(str(field) for field in snp) + "\n")

    # a track for IGV. Both contigs, so it loads against combined.fa,
    # recipient.fa or bac.fa without conversion.
    with open(os.path.join(outdir, "snps.bed"), "w") as out:
        out.write('track name="diagnostic_snps" '
                  'description="where recipient and BAC differ" '
                  'useScore=1 visibility=2\n')
        rows = []
        for snp in snps:
            rec_pos, bac_pos = snp[0], snp[1]
            rows.append(("recipient", rec_pos - 1, rec_pos,
                         f"{snp[3]}>{snp[4]}"))
            rows.append(("bac", bac_pos - 1, bac_pos, f"{snp[5]}>{snp[6]}"))
        for contig, start, end, label in sorted(rows):
            out.write(f"{contig}\t{start}\t{end}\t{label}\t1000\t.\n")

    syntenic_bp = sum(a["rec_end"] - a["rec_start"] for a in alignments)
    # fields 3 and 4 are the two alleles in recipient orientation
    transitions = sum(1 for snp in snps
                      if snp[3] + snp[4] in ("AG", "GA", "CT", "TC"))
    transversions = len(snps) - transitions

    print(f"[find_synteny] {len(alignments)} syntenic blocks, "
          f"{syntenic_bp:,} bp of the recipient"
          f"  ({len(every_block)} homologous blocks in all)")
    print(f"[find_synteny] {len(snps):,} diagnostic SNPs "
          f"(1 per {syntenic_bp // max(len(snps), 1):,} bp)")
    # A real pair of genomes gives Ti/Tv around 2-4. A value near 0.5 means the
    # strand handling above has gone wrong and the alleles are being compared
    # across orientations. Cheap check, so it runs every time.
    if transversions:
        ratio = transitions / transversions
        flag = "" if ratio > 1.5 else "   <-- SUSPICIOUS, expected > 1.5"
        print(f"[find_synteny] Ti/Tv = {ratio:.2f}{flag}")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1], args[2],
         int(args[3]) if len(args) > 3 else 1000)
