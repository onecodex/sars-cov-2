#!/usr/bin/env python

"""Generate summary statistics of read coverage by region defined in a bedfile.
Use to evaluate per-amplicon performance of primers for a sample."""

import argparse
from io import StringIO
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bedfile", help="Amplicon insert regions bedfile path")
    parser.add_argument("alignments", help="Alignments .bam path")
    parser.add_argument(
        "--prop_cov",
        help="Proportion of a read that must be mapped to a region "
        "for it to count towards coverage. Default: 0.90",
        default="0.90",
    )
    args = parser.parse_args()
    return args


def run(args):
    # Run bedtools to get per-base coverage for reads that map to each insert
    cmd = [
        "bedtools",
        "coverage",
        "-F",
        args.prop_cov,
        "-d",
        "-a",
        args.bedfile,
        "-b",
        args.alignments,
    ]
    a = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    b = StringIO(a.communicate()[0].decode("utf-8"))
    df = pd.read_csv(
        b,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "chrom_start",
            "chrom_end",
            "region_name",
            "primer_pool",
            "strand",
            "position",
            "depth",
        ],
    )

    # Generate summary statistics of coverage by insert and save to file
    insert_stats_data = []
    for insert in df["region_name"].unique():
        insert_depth_stats = df[df["region_name"] == insert]["depth"].describe()
        insert_depth_stats["insert_name"] = str(insert)
        insert_stats_data.append(insert_depth_stats)
    insert_stats_df = pd.DataFrame(insert_stats_data).reset_index(drop=True)
    insert_stats_df.rename(columns={"count": "insert_length"}, inplace=True)
    insert_stats_df = insert_stats_df[
        ["insert_name", "insert_length", "mean", "std", "min", "25%", "50%", "75%", "max"]
    ]
    insert_stats_df.to_csv("depth_by_insert_stats.tsv", index=False, sep="\t")

    # Save a boxplot of coverage by insert
    ax, fig = plt.subplots(figsize=[18, 5])
    sns.boxplot(x="region_name", y="depth", data=df, color="lightcyan", linewidth=1)
    plt.xticks(rotation=90, size=8)
    plt.savefig("depth_by_insert_boxplot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    args = parse_args()
    run(args)
