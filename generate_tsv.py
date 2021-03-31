import pandas as pd

colnames = [
    "region",
    "position",
    "ref allele",
    "alt allele",
    "frequencies",
    "protein sequence variant",
]

df_variants = pd.read_csv("snpeff.tsv", skiprows=[0], names=colnames, sep="\t| ", engine="python")

for i in df_variants.index:
    ref_reads = int(df_variants.loc[i, "frequencies"].rsplit(",", 3)[0]) + int(
        df_variants.loc[i, "frequencies"].rsplit(",", 3)[1]
    )
    alt_reads = int(df_variants.loc[i, "frequencies"].rsplit(",", 3)[2]) + int(
        df_variants.loc[i, "frequencies"].rsplit(",", 3)[3]
    )
    all_reads = ref_reads + alt_reads
    alt_freq = alt_reads / (all_reads) * 100
    df_variants.loc[i, "depth"] = "{:0.0f}".format(all_reads)
    if alt_freq == 100:
        df_variants.loc[i, "alt frequency"] = "100%"
    else:
        df_variants.loc[i, "alt frequency"] = "{:0.2f}%".format(alt_freq)

df_variants = df_variants.drop(columns=["frequencies"])


df_orfs = pd.read_csv(
    "/annot_table.orfs.txt",
    sep="\t",
    header=None,
    usecols=[0, 1, 2],
    names=["orf", "start", "stop"],
)
for i in df_variants.index:
    for j in df_orfs.index:
        if df_orfs.loc[j, "start"] <= df_variants.loc[i, "position"] <= df_orfs.loc[j, "stop"]:
            df_variants.loc[i, "orf"] = df_orfs.loc[j, "orf"]

# os.path.join(workdir, "variants.tsv")
df_variants.to_csv("variants.tsv", sep="\t")
