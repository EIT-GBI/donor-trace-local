# donor-trace-local

Trace donor sequence through a recipient genome, allele by allele.

For each sequenced clone: **is the donor sequence (here, a BAC cargo) there,
and across how much of the region?**

## Approach: three tiers, three different questions

No single alignment can answer "is the donor sequence there, and is it
integrated". The pipeline runs three, each against a reference chosen to make
one question answerable, and each tier is independent of the others — they can
agree, and when they disagree that is information.

The constraint behind all of it: in the initial experiments the BAC cargo is
~99% identical to the recipient locus it targets, about one difference every
100–200 bp. If a reference holds **both** copies, a read from the shared
sequence fits each equally well, gets MAPQ 0, and is assigned to one of them at
random. Everything below follows from working around that.

### Tier 1 — combined reference: look at the alignment

Both sequences, as two contigs. Align, make depth tracks, look at it.

| reads here are | so the depth is |
|---|---|
| BAC-only sequence (vector) | **trustworthy** — nothing competes for these reads, so this is a clean "is a BAC present at all, and how much of the backbone survives" |
| recipient-only sequence | **trustworthy** — chromosome coverage, and where it drops or breaks |
| the shared region | **not trustworthy per contig** — MAPQ 0, split at random between the two. The *sum* over both contigs is still meaningful; neither contig alone is |



### Tier 2 — single references: look at the allele frequency

Align separately to `recipient.fa` and to `bac.fa`. Each holds only **one**
copy of the shared sequence, so there is nothing to assign: recipient-derived
and BAC-derived reads stack at the same position, and the allele fraction at
the diagnostic SNPs is a direct measure of the mixture.

| BAC allele fraction | meaning |
|---|---|
| ~0.0 | only recipient sequence here — no cargo |
| ~0.5 | one BAC copy present alongside the recipient locus |
| ~1.0 | the recipient locus itself has been replaced |

This is the quantitative tier: how much cargo, and across exactly which
interval. Running it from both sides is deliberate — the two share no alignment
step, so agreement is a real check rather than the same number restated.

### Tier 3 — junction reads: look for integration

Tiers 1 and 2 both measure how much cargo is present. Neither can say whether
it is **integrated**, because an episome at one copy per cell and an integrated
copy give the same depth and the same allele fraction. The only thing that
distinguishes them is a single molecule containing both sequences — and a read
off such a molecule cannot align in one piece. It aligns as a primary plus
supplementary alignments, which on the two-contig reference land on different
contigs.

So this tier reuses the Tier 1 BAM (supplementary alignments are kept there for
exactly this reason) and asks which of those split reads are real. Most are
not: the cargo matches both contigs, so reads from it split for alignment
reasons. Separating them needs the synteny map *and* the Tier 2 allele calls —
which is why it runs last.

## What it does

1. **Build the references.** Three per BAC, in `$OUTDIR/reference/<BAC>/`:
   - `combined.fa` — recipient and BAC as **two contigs**, `recipient` and
     `bac`. No spacer, no offsets: a BAC position is just a position on `bac`,
     so every BED loads against any of the three references unchanged.
   - `recipient.fa` — the recipient alone
   - `bac.fa` — the BAC alone

   The last two are not a spare set. Tier 2 needs a reference holding **one**
   copy of the shared sequence; `combined.fa` holds both.
2. **Find synteny.** Align BAC to recipient, record the homologous blocks
   (`synteny.bed`, `homology_all.bed`) and every single-base difference inside
   them (`snps.tsv`, `snps.bed`).
3. **Tier 1 — align three times** per sample, to `combined`, `recipient` and
   `bac`, with a bigWig from each. Secondary alignments are always dropped.
   Supplementary alignments are **kept for `combined`** (Tier 3 needs them) and
   dropped for the other two, where each read should have exactly one home.
4. **Tier 2 — count alleles twice**, once on each single-reference alignment.
5. **Plot** both, and write `summary.tsv`.
6. **Tier 3 — look for junctions**, off the `combined` BAM from step 3. The
   candidate reads are pulled into `junctions/<sample>.junction.bam` first and
   the classifier reads that — 308 alignment records instead of 227,704 on this
   data, which is the difference between 27 MB and 4.2 GB of memory. Writes
   `summary_junctions.tsv`.

## Running it

```bash
bash setup.sh                 # once
# edit config.txt and samplesheet_ont.tsv
bash run.sh config.txt
```

`setup.sh` is safe to re-run and does three things:

1. `brew install minimap2 bwa samtools bcftools seqkit fastp`
2. installs [uv](https://docs.astral.sh/uv/) if missing, then `uv sync` to
   build `.venv` from `pyproject.toml` (matplotlib, for the plots)
3. downloads `bedGraphToBigWig` into `tools/ucsc/` — UCSC ships it as a binary
   and Homebrew does not carry it

`run.sh` activates `.venv` itself, so there is nothing to activate by hand.
Note the download host is `hgdownload.**soe**.ucsc.edu`; the `cse` name that
most instructions still quote now fails certificate validation.

One platform per run. For both, keep two configs:

```bash
bash run.sh config_ont.txt
bash run.sh config_illumina.txt
```

### Sample sheet

Tab-separated, with this header. **The columns differ by platform** — ONT has
no `fastq2`:

`samplesheet_ont.tsv`

```
sample      fastq1               recipient       bac
SF1-1.2A    /path/SF1.fastq.gz   /path/NZ131.fa  /path/SF1_BAC.fa
```

`samplesheet_illumina.tsv`

```
sample      fastq1               fastq2               recipient       bac
RN_1_S17    /path/RN_1_R1.fq.gz  /path/RN_1_R2.fq.gz  /path/NZ131.fa  /path/SF1_BAC.fa
```

Do not leave an empty `fastq2` column in an ONT sheet. Bash treats a tab as
whitespace, so two consecutive tabs collapse into one and every later field
shifts left — which is why ONT gets its own four-column layout rather than a
blank third column.

`recipient` and `bac` are paths to FASTA files. Samples sharing a BAC share its
reference, which is built once.

## Output

```
$OUTDIR/
  reference/<BAC>/   combined.fa, recipient.fa, bac.fa,
                     synteny.bed, homology_all.bed, snps.bed, snps.tsv
  bam/               <sample>.combined.bam        tier 1 (+ tier 3 input)
                     <sample>.{recipient,bac}.bam tier 2
  bigwig/            <sample>.{combined,recipient,bac}.bw   tier 1
  alleles/           <sample>.{recipient,bac}.tsv and .bedgraph   tier 2
  plots/             <sample>.{recipient,bac}.png                 tier 2
  junctions/         <sample>.junction.bam   reads on BOTH contigs, unfiltered
                     <sample>.junctions.tsv  per-piece verdicts      tier 3
  summary.tsv            tier 2, one row per sample per frame
  summary_junctions.tsv  tier 3, one row per sample
```

### Reading `summary.tsv` (Tier 2)

| column | meaning |
|---|---|
| `frame` | which alignment it came from — `recipient` or `bac`. **The two rows for one sample should agree.** If they don't, something is wrong upstream, not biologically interesting |
| `sites` | diagnostic SNPs with enough depth to score |
| `bac_fraction` | BAC allele fraction over the whole homologous region |
| `error` | third-allele rate: reads matching *neither* expected allele, so sequencing error, measured at the same positions. **A sample with no BAC should land near this number, not at zero.** If `bac_fraction` is only 2–3× `error`, that is noise |
| `tract` | longest unbroken run of recipient sequence inside a sample that carries BAC elsewhere — a block of cargo replaced by recipient sequence |

### Tier 2: the two frames are not equally powered

They measure the same thing, but not with the same confidence, and the
difference is structural rather than noise.

The **recipient frame** always has reads: the recipient chromosome is present
in every sample, so every diagnostic site is covered whether or not a BAC is
there. This is the primary readout.

The **BAC frame** only has reads where BAC sequence exists. In a sample with no
cargo, the few reads that land on `bac.fa` are mismapped ones, so the fraction
is computed from almost nothing and comes out noisy — on the test data,
`bac_fraction` 0.0023 in the BAC frame against 0.0003 in the recipient frame
for the same cargo-free sample. Both say "absent"; only one of them is
well-powered.

There is also a mild bias in the BAC frame even when cargo IS present.
Recipient-derived reads have to map to the BAC reference, and the most
divergent ones fail, which pushes the BAC fraction up — 0.63 in the BAC frame
versus 0.55 in the recipient frame for the same clone. Expect the BAC frame to
read a few points high.

So: **take the number from the recipient frame, and use the BAC frame to check
that the shape agrees.** Tract boundaries should match to within a few hundred
bp. If the two frames disagree about whether cargo is present at all, that is a
problem with the reference or the sample sheet, not a biological result.

### Reading the plots (Tier 2)

The two lines are complementary, so it is one measurement drawn twice; the
crossovers are just easier to see that way. A long flat stretch at 1.0
recipient inside an otherwise 50/50 region is cargo **replaced by recipient
sequence**. Before calling that: converted sequence keeps its depth, deleted
sequence does not — so check the bigWig.

### In IGV (Tier 1)

Load `reference/<BAC>/combined.fa`, then `bigwig/<sample>.combined.bw` for
depth (Tier 1), `alleles/<sample>.recipient.bedgraph` for the allele track
(Tier 2), and `junctions/<sample>.junction.bam` for the candidate junction
reads (Tier 3). The reference folder's `synteny.bed`, `homology_all.bed` and
`snps.bed` all load against the same coordinates — no conversion, which is the
point of the two-contig layout.

Keep `homology_all.bed` on screen. It marks where the two contigs share
sequence, and therefore where Tier 1 depth means nothing. **Fix the
allele track to 0–1**; on autoscale a noise floor of 0.0005 fills the panel and
looks like signal.

## Reading `summary_junctions.tsv` (Tier 3)

`summary_junctions.tsv` has two numbers per sample and they mean different
things:

| column | meaning |
|---|---|
| `candidates` | reads with alignments on **both** contigs. Largely meaningless on its own — the cargo matches both, so reads from it split and land on each |
| `junction_reads` | candidates that survive both tests. **This is the integration count** |
| `piece_verdicts` | what disqualified the rest, so a zero can be read rather than taken on trust |

A candidate survives only if it has a piece on the BAC lying **outside** every
syntenic block — vector, which has no recipient counterpart — *and* a piece
that is really the recipient genome. For the second, a piece outside the
syntenic blocks is accepted as-is; a piece inside them is **genotyped**:
recipient alleles mean the recipient genome, donor alleles mean cargo the
aligner put there.

Genotyping rather than discarding matters. Filtering on synteny alone throws
away every recipient-side alignment in the target window — and if a BAC had
integrated, that is exactly where the recipient half of the junction read
would be. The test would be unable to see the thing it is testing for.

**Only SNPs that are still informative in that sample are used**, which is why
this step runs after the allele step and reads `<sample>.recipient.tsv`. A
listed SNP inside a conversion tract tells you nothing — both molecules read
recipient there — and scoring it anyway makes BAC reads look like recipient
reads and invents junctions. That table has to come from the recipient-only
alignment: measured on the two-contig BAM, cargo reads go to the `bac` contig
and the recipient side looks almost pure.

`junctions/<sample>.junction.bam` holds every candidate with **no** filtering,
all alignments including supplementary, for looking at in IGV against
`reference/<BAC>/combined.fa` — with `synteny.bed`, `homology_all.bed` and
`snps.bed` from the same folder as tracks.

A non-zero count is not proof by itself. Check that the positions cluster,
check MAPQ, and check the read actually spans the join rather than ending at
it.

**The homology map has to be complete.** This step reads
`reference/<BAC>/homology_all.bed`, not `synteny.bed`: every scrap of shared
sequence, found with a deliberately over-sensitive search and no length
filter. A short shared element is far too small to be cargo but is perfectly
capable of making a read align to both sequences. Erring towards too many
blocks is the safe direction — an extra one only makes the test more
conservative, a missing one invents junctions. A real case from this data: a
731 bp element at 81% identity, shared between the miniREXER vector and the
chromosome, is found by `minimap2 -x asm20` but dropped the moment `-c` is
added. That one missing block turned 103 reads into false junctions.

## Comparing the tiers

They are independent, so disagreement is worth reading rather than averaging
away:

| | what it says | typical disagreement, and what it means |
|---|---|---|
| Tier 1 depth on BAC-only sequence | a BAC is present | full vector depth with Tier 2 at ~0 → **vector retained, cargo lost** |
| Tier 2 allele fraction | how much cargo, and where | a flat run at 0.0 inside an otherwise 0.5 region → a block of cargo **replaced by recipient sequence**, not deleted. Check Tier 1 depth: converted sequence keeps its depth, deleted sequence does not |
| Tier 3 junction reads | integrated or episomal | Tier 2 at ~0.5 with zero junctions → **episomal**. Tier 2 at ~1.0 with junctions → **integrated** |

A Tier 1 breakpoint with nothing in Tiers 2 or 3 usually means the recipient
reference is wrong, not that anything happened to the BAC.

## The strand trap

A BAC often aligns to the recipient on the **minus strand** — all three SF BACs
here do. A base written in BAC orientation is then the complement of the same
base written in recipient orientation. Compare them directly and A↔T and C↔G
differences vanish into "no difference", while every transition reads as a
transversion. The resulting allele calls are wrong at every site, and they fail
quietly: both copies look perfectly pure.

`snps.tsv` is written so this cannot happen. All four alleles are stored
explicitly, two per coordinate frame:

```
rec_pos  bac_pos  strand  rec_allele_rec_frame  bac_allele_rec_frame
                          rec_allele_bac_frame  bac_allele_bac_frame
```

Score a recipient alignment with the `*_rec_frame` pair, a BAC alignment with
the `*_bac_frame` pair. Nothing ever needs complementing. `find_synteny.py`
prints Ti/Tv on every run as a check — real genomes give 2–4, and anything near
0.5 means orientations have been mixed.

Note this also means **coordinates run backwards** on a minus-strand BAC. The
same block is "the start of the cargo" in BAC coordinates and "the far end of
the locus" in recipient coordinates. The two plots will look mirrored. They
agree.

## Scripts

| script | what it does |
|---|---|
| `setup.sh` | brew packages, uv + `.venv`, `bedGraphToBigWig` |
| `run.sh` | the driver — activates `.venv`, loops over samples |
| `scripts/prepare_reference.sh` | three FASTAs per BAC, indexed |
| `scripts/find_synteny.py` | minimap2 BAC vs recipient → `synteny.bed`, `homology_all.bed`, `snps.bed`, `snps.tsv` |
| `scripts/align.sh` | one sample, one reference → sorted indexed BAM |
| `scripts/make_bigwig.sh` | BAM → bigWig depth track |
| `scripts/allele_freq.py` | pileup at the diagnostic SNPs → allele fractions |
| `scripts/plot_alleles.py` | allele table → PNG |
| `scripts/summarise.py` | one table across all samples, finds tracts |
| `scripts/junction_bam.sh` | reads with alignments on both contigs → BAM |
| `scripts/junction_reads.py` | which of those are real junctions |
| `scripts/summarise_junctions.py` | one junction table across all samples |

## Requirements

`bash setup.sh` installs everything. By hand it is: Homebrew for `minimap2`,
`bwa`, `samtools`, `bcftools`, `seqkit`, `fastp`; `uv` for the Python
environment; and the UCSC `bedGraphToBigWig` binary in `tools/ucsc/`.

Without `bedGraphToBigWig` the pipeline still runs — it just makes no tracks.
