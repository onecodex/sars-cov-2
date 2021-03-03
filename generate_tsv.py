import sys
import pandas as pd

samplename = sys.argv[1]
wd = "/" + samplename + "/"
print(wd)

df_variants = pd.DataFrame(
    columns=[
        "region",
        "position",
        "ref",
        "alt",
        "depth",
        "alt frequency",
        "S/N",
        "orf",
        "AA mutation",
    ]
)

# Import the Medaka-generated vcf file. These are our variants
df_variants = pd.read_csv(
    wd + samplename + ".pass.snpeff.vcf.noheaders",
    sep="\t",
    header=None,
    usecols=[0, 1, 3, 4, 7],
    names=["region", "position", "ref", "alt", "info_vcf"],
)

# Import the bcf file. The info column will give us the allele frequency
df_bcf = pd.read_csv(
    wd + samplename + ".sorted.bcf.noheaders",
    sep="\t",
    header=None,
    usecols=[0, 1, 3, 4, 7],
    names=["region", "position", "ref", "alt", "info_bcf"],
)
df_bcf = df_bcf[df_bcf.alt != "."]

# Import the variant headers from the snpEff-generated faa file. This will give us AA mutations
aa = open(wd + samplename + ".pass.snpeff.vcf.faa.headers", "r")
for line in aa.readlines():
    pos = int(line.split(":")[1].split("-")[0])
    aa_mut = line.split("p:p.")[1].rstrip("\n")
    ind = df_variants[df_variants["position"] == pos].index
    df_variants.loc[ind, "AA mutation"] = aa_mut

# Add in depth, allele frequency, S/N info
for i in df_variants.index:
    info_bcf = str(df_bcf[df_bcf.position == df_variants.loc[i, "position"]]["info_bcf"])
    depth = info_bcf.split()[1].split(";")[0].split("=")[1]
    df_variants.loc[i, "depth"] = depth
    alt_freq = info_bcf.split()[1].split(";")[1].split("=")[1].split(",")[1]
    df_variants.loc[i, "alt frequency"] = alt_freq
    var_effect = df_variants.loc[i, "info_vcf"].split("|")[1]
    if var_effect == "synonymous_variant":
        var_type = "S"
    elif var_effect == "missense_variant":
        var_type = "N"
    else:
        var_type = "non-coding"
    df_variants.loc[i, "S/N"] = var_type

df_variants = df_variants.drop(columns=["info_vcf"])

# Add in orf info
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

# Change AA mutations from 3-letter to 1-letter codes
aa_code = pd.read_csv("/aa_codes.txt", sep="\t", index_col="Three-Letter Code")

for i in df_variants.index:
    aa_mut = str(df_variants.loc[i, "AA mutation"])
    if aa_mut != "nan":
        aa_ref = aa_mut[0:3]
        aa_ref = aa_code.loc[aa_ref, "One-Letter Code"]
        aa_alt = aa_mut[-3:]
        aa_alt = aa_code.loc[aa_alt, "One-Letter Code"]
        aa_num = ""
        for s in aa_mut:
            if s.isdigit():
                aa_num = aa_num + s
        df_variants.loc[i, "AA mutation"] = aa_ref + aa_num + aa_alt

df_variants.to_csv("/variants.tsv", sep="\t")
