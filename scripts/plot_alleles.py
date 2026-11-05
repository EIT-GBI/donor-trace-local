#!/usr/bin/env python3
"""Plot recipient vs BAC allele frequency across the homologous region.

  usage: plot_alleles.py <allele_freq.tsv> <recipient|bac> <title> <out.png>
                         [smooth_window] [recipient_label] [bac_label]

One line per allele. They are complementary - they add to 1 at every point -
so the plot is really one measurement drawn twice. Drawn that way because the
crossover points are much easier to see.

How to read it:
  both lines near 0.5    BAC and recipient sequence both present, one copy each
  recipient line at 1.0  only recipient sequence here; no BAC cargo
  BAC line at 1.0        the recipient locus has been replaced by BAC sequence

A long flat stretch at 1.0 inside an otherwise 50/50 region is a block of the
cargo that is present as recipient sequence - converted, not deleted. Check the
depth track before calling it: converted sequence keeps its depth, deleted
sequence does not.

smooth_window is a rolling mean in NUMBER OF SNPS, not base pairs, because the
SNPs are not evenly spaced. 1 disables smoothing.
"""
import sys

import matplotlib
matplotlib.use("Agg")                   # no display needed
import matplotlib.pyplot as plt

# A stretch with no diagnostic SNPs is not a measurement, but a line plot will
# happily draw straight through it and make it look like one. Break the line
# instead wherever consecutive SNPs are further apart than this.
MAX_GAP_BP = 5000

RECIPIENT_COLOUR = "#1b7f4b"
BAC_COLOUR = "#e8468c"


def read_table(path):
    positions, bac_fractions = [], []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            field = line.rstrip("\n").split("\t")
            positions.append(int(field[index["pos"]]))
            bac_fractions.append(float(field[index["bac_fraction"]]))
    return positions, bac_fractions


def rolling_mean(values, window):
    """Mean of up to `window` neighbouring values, centred on each point."""
    if window <= 1:
        return list(values)
    half = window // 2
    smoothed = []
    for i in range(len(values)):
        start = max(0, i - half)
        end = min(len(values), i + half + 1)
        smoothed.append(sum(values[start:end]) / (end - start))
    return smoothed


def break_at_gaps(positions, values):
    """Insert a gap in the line wherever the SNPs stop.

    matplotlib skips any point that is NaN, so putting one between two distant
    positions leaves a visible hole rather than a straight line across sequence
    that was never scored.
    """
    x_out, y_out = [], []
    for i, (position, value) in enumerate(zip(positions, values)):
        if i > 0 and position - positions[i - 1] > MAX_GAP_BP:
            x_out.append(float("nan"))
            y_out.append(float("nan"))
        x_out.append(position / 1000.0)
        y_out.append(value)
    return x_out, y_out


def main(table, frame, title, out_png, smooth_window=11,
         recipient_label="recipient", bac_label="BAC"):
    positions, bac_fractions = read_table(table)
    if not positions:
        print(f"[plot_alleles] no scored sites, skipping {out_png}")
        return

    bac_smoothed = rolling_mean(bac_fractions, smooth_window)
    recipient_smoothed = [1.0 - value for value in bac_smoothed]

    x_kb, bac_line = break_at_gaps(positions, bac_smoothed)
    _, recipient_line = break_at_gaps(positions, recipient_smoothed)

    figure, axes = plt.subplots(figsize=(16, 5))
    axes.plot(x_kb, recipient_line, color=RECIPIENT_COLOUR, linewidth=1.2,
              label=recipient_label)
    axes.plot(x_kb, bac_line, color=BAC_COLOUR, linewidth=1.2, label=bac_label)

    # 0.5 is where a single episomal copy sits; worth a guide line.
    axes.axhline(0.5, color="#999999", linewidth=0.6, linestyle=":", zorder=0)

    axes.set_ylim(-0.03, 1.03)
    axes.set_xlim(min(positions) / 1000.0, max(positions) / 1000.0)
    label = "BAC" if frame == "bac" else "recipient"
    axes.set_xlabel(f"{label} position (kb)")
    axes.set_ylabel("relative allele frequency")
    # title left, legend above right - neither sits on top of the data
    axes.set_title(title, loc="left")
    axes.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0),
                frameon=False, ncol=2)
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)

    figure.tight_layout()
    figure.savefig(out_png, dpi=150)
    plt.close(figure)
    print(f"[plot_alleles] {out_png}  ({len(positions)} sites,"
          f" smoothed over {smooth_window})")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args[0], args[1], args[2], args[3],
         int(args[4]) if len(args) > 4 else 11,
         args[5] if len(args) > 5 else "recipient",
         args[6] if len(args) > 6 else "BAC")
