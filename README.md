# donor-trace-local

Is the donor sequence there, across how much of the region, and is it
integrated?

## Why three tiers

The cargo is ~99% identical to the recipient locus it targets — about one
difference every 100–200 bp. A read from that shared sequence fits both
copies, so on a reference holding both it gets MAPQ 0 and is assigned at
random. Everything below works around that.

| tier | reference | measurement | answers |
|---|---|---|---|
| **1** | `combined.fa` — both, two contigs | depth | is a BAC there at all; where the reference is wrong |
| **2** | `recipient.fa` and `bac.fa` separately | allele fraction at diagnostic SNPs | how much cargo, over which interval |
| **3** | tier 1 BAM | junction reads | integrated or episomal |

**Tier 1.** Depth is trustworthy on BAC-only sequence (vector) and on
recipient-only sequence. It is *not* trustworthy on the shared region — MAPQ 0,
split at random between the contigs. Use it for the overview and for spotting
structural surprises, not for cargo.

**Tier 2.** One copy of the shared sequence per reference, so nothing is
assigned: recipient- and BAC-derived reads stack at the same position and the
allele fraction is the mixture.

| BAC allele fraction | meaning |
|---|---|
| ~0.0 | no cargo |
| ~0.5 | one BAC copy alongside the recipient locus |
| ~1.0 | the recipient locus has been replaced |

Run from both sides. They share no alignment step, so agreement is a real
check.

**Tier 3.** An episome at one copy per cell gives the same depth and the same
allele fraction as an integrated copy. Only a molecule containing both
sequences separates them, and a read off one aligns in pieces that land on
different contigs — which is why tier 1 keeps supplementary alignments.

## Running it
Install all the dependencies and set up the environment with:
```bash
bash setup.sh                 # once
```
Create your samplesheet.csv file with long or short reads and pass it to the config.txt, i.e. specify it there.
Edit any other settings you may want in `config.txt` (like output directory) and then run the pipeline with:
```bash
# edit config.txt and samplesheet_ont.csv
bash run.sh config.txt
```

One platform per run; keep two configs. `setup.sh` installs the tools via
Homebrew, installs `uv` and builds `.venv`, and downloads `bedGraphToBigWig`
into `tools/ucsc/`. `run.sh` activates `.venv` itself.

### Sample sheet

CSV. **Columns differ by platform** — ONT has no `fastq2`:

```
samplesheet_ont.csv
sample,fastq1,recipient,bac
SF1-1.2A,/path/SF1.fastq.gz,/path/NZ131.fa,/path/SF1_BAC.fa

samplesheet_illumina.csv
sample,fastq1,fastq2,recipient,bac
RN_1_S17,/path/R1.fq.gz,/path/R2.fq.gz,/path/NZ131.fa,/path/SF1_BAC.fa
```

Samples sharing a BAC share its reference, which is built once.

## Output

```
$OUTDIR/
  reference/<BAC>/  combined.fa recipient.fa bac.fa
                    synteny.bed homology_all.bed snps.bed snps.tsv
  bam/              <sample>.combined.bam          tier 1, tier 3 input
                    <sample>.{recipient,bac}.bam   tier 2
  bigwig/           <sample>.*.bw                  tier 1
  alleles/          <sample>.{recipient,bac}.tsv / .bedgraph   tier 2
  plots/            <sample>.{recipient,bac}.png   tier 2
  junctions/        <sample>.junction.bam          candidates, unfiltered
                    <sample>.junctions.tsv         per-piece verdicts
  summary.tsv            tier 2
  summary_junctions.tsv  tier 3
```

### `summary.tsv`

| column | meaning |
|---|---|
| `frame` | `recipient` or `bac`. **The two rows per sample should agree** |
| `sites` | diagnostic SNPs with enough depth |
| `bac_fraction` | BAC allele fraction over the region |
| `error` | third-allele rate — sequencing error at the same positions. A cargo-free sample lands near this, not at zero |
| `tract` | longest run of recipient sequence inside a sample that carries BAC elsewhere: cargo replaced by recipient sequence |

The two frames are not equally powered. The recipient frame always has reads;
the BAC frame only has reads where cargo exists, so with no cargo it is
computed from almost nothing (0.0023 vs 0.0003 on real data — both "absent",
one meaningless). It also reads a few points high when cargo is present,
because divergent recipient reads fail to map to the BAC. **Take the number
from the recipient frame; use the BAC frame to check the shape.**

### `summary_junctions.tsv`

| column | meaning |
|---|---|
| `candidates` | reads with alignments on both contigs. Largely meaningless alone — cargo reads split and land on each |
| `junction_reads` | candidates that survive. **This is the integration count** |
| `usable` | `ok`, or `UNINFORMATIVE` if no piece ever reached the anchor length — then a zero means "not tested", not "negative" |
| `piece_verdicts` | what disqualified the rest, so a zero can be read rather than trusted |

**Tier 3 needs long reads.** A 150 bp Illumina read cannot carry two
`MIN_ANCHOR`-length pieces, so every piece comes back `TOO_SHORT` and the
sample is reported `UNINFORMATIVE` rather than negative. Lowering the anchor
far enough to admit short reads buys false positives, not sensitivity. Tiers 1
and 2 work normally on Illumina.

A candidate survives only with a piece on the BAC outside every homologous
block *and* a piece that is really the recipient genome. For the second, a
piece inside the homologous blocks is **genotyped** rather than discarded —
filtering on synteny alone throws away the target window, which is exactly
where the recipient half of a real junction read would be.

Only SNPs still informative in that sample are used, which is why this runs
after tier 2 and reads `<sample>.recipient.tsv`. A SNP inside a conversion
tract tells you nothing — both molecules read recipient — and scoring it
anyway invents junctions.

It also reads `homology_all.bed`, not `synteny.bed`: every scrap of shared
sequence, unfiltered. A 731 bp element shared by vector and chromosome once
turned 103 reads into false junctions.

### Plots

Two complementary lines, so it is one measurement drawn twice; the crossovers
are easier to see that way. A flat run at 1.0 recipient inside an otherwise
50/50 region is cargo **replaced by recipient sequence** — check the bigWig,
since converted sequence keeps its depth and deleted sequence does not.

### IGV

Load `reference/<BAC>/combined.fa`, then `bigwig/<sample>.combined.bw`,
`alleles/<sample>.recipient.bedgraph` and `junctions/<sample>.junction.bam`.
The reference folder's BEDs load against the same coordinates — no conversion,
which is the point of the two-contig layout. Keep `homology_all.bed` on
screen: it marks where tier 1 depth means nothing. Fix the allele track to
0–1, or autoscale makes a 0.0005 noise floor look like signal.

### Comparing tiers

| observation | meaning |
|---|---|
| full vector depth (T1), ~0 allele fraction (T2) | vector retained, cargo lost |
| flat 0.0 run inside a 0.5 region (T2) | cargo replaced by recipient sequence |
| ~0.5 (T2), zero junctions (T3) | episomal |
| ~1.0 (T2), junctions (T3) | integrated |
| breakpoint in T1, nothing in T2 or T3 | the recipient reference is wrong |

## The strand trap

A BAC often aligns to the recipient on the **minus strand**. A base in BAC
orientation is then the complement of the same base in recipient orientation.
Compare them directly and A↔T and C↔G differences vanish, every transition
reads as a transversion, and the allele calls are wrong at every site —
quietly, with both copies looking pure.

`snps.tsv` stores all four alleles explicitly so this cannot happen:

```
rec_pos  bac_pos  strand  rec_allele_rec_frame  bac_allele_rec_frame
                          rec_allele_bac_frame  bac_allele_bac_frame
```

Score a recipient alignment with the `*_rec_frame` pair, a BAC alignment with
the `*_bac_frame` pair. Nothing needs complementing. `find_synteny.py` prints
Ti/Tv every run — real genomes give 2–4, anything near 0.5 means orientations
are mixed.

Coordinates run backwards too, so the same block is "the start of the cargo"
in BAC coordinates and "the far end of the locus" in recipient coordinates.
The two plots look mirrored. They agree.

## Scripts

| script | does |
|---|---|
| `setup.sh` | brew packages, uv + `.venv`, `bedGraphToBigWig` |
| `run.sh` | the driver |
| `scripts/prepare_reference.sh` | three FASTAs per BAC, indexed |
| `scripts/find_synteny.py` | BAC vs recipient → synteny, homology, SNPs |
| `scripts/align.sh` | one sample, one reference → sorted BAM |
| `scripts/make_bigwig.sh` | BAM → bigWig |
| `scripts/allele_freq.py` | pileup at the SNPs → allele fractions |
| `scripts/plot_alleles.py` | allele table → PNG |
| `scripts/summarise.py` | tier 2 table, finds tracts |
| `scripts/junction_bam.sh` | reads on both contigs → BAM |
| `scripts/junction_reads.py` | which of those are real junctions |
| `scripts/summarise_junctions.py` | tier 3 table |

## Requirements

`bash setup.sh` installs everything: Homebrew for `minimap2`, `bwa`,
`samtools`, `bcftools`, `seqkit`, `fastp`; `uv` for the Python environment;
and the UCSC `bedGraphToBigWig` binary in `tools/ucsc/`. Without the last one
the pipeline runs but makes no tracks.
