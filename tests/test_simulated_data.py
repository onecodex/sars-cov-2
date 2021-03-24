from Bio import SeqIO
import pandas as pd
import pytest


def read_vcf_as_dataframe(path):
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        usecols=[0, 1, 2, 3, 4, 5, 6, 7],
    )


@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only_fastq(tmp_path, n, run_art, run_covid_pipeline, run_snp_mutator):
    """Tests insert of N snps
    """
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=n)

    # Run ART
    run_art()

    # Run pipeline on simulated data
    run_covid_pipeline(input_filename="simulated_reads.fq")

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = read_vcf_as_dataframe(tmp_path / "variants.vcf")

    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])

    # We add these tests to ensure we have a high percent of reads aligning
    # We simulate at 50x, so low end variants with coverage variability should
    # be ~25-30x, and then another ~33-50% due to Q scores <20

    # parse alt depth from VCF info column
    alt_dp = []
    for row in called["INFO"].values:
        items = row.split(";")
        for i in items:
            if i.startswith("DP="):
                alt_dp.append(int(i.split("=")[1]))

    assert len(alt_dp) == called.shape[0]
    assert all([i > 10 for i in alt_dp])
    assert sum(alt_dp) / len(alt_dp) > 15

    # Finally, test that the FASTAs match
    # Note we ignore the first 50bp which may have low coverage and N masking
    # plus the final 120bps due to a polyA tail
    reference = list(SeqIO.parse(f"{tmp_path}/nCoV-2019.reference_mutated_1.fasta", "fasta"))[0]
    consensus = list(SeqIO.parse(f"{tmp_path}/consensus.fa", "fasta"))[0]
    assert consensus.seq[50:-120] in reference.seq
