import pandas as pd
import pytest


@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only_fasta(tmp_path, n, run_covid_pipeline, run_snp_mutator):
    """Tests insert of N snps
    """
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=n)

    # Run pipeline on simulated data
    run_covid_pipeline()

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])


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
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])

    # Note: `ivar trim -e` on WGS removes ~50% of the data
    # We add these tests to ensure we have a high percent of reads aligning
    # We simulate at 50x, so low end variants with coverage variability should
    # be ~25-30x, and then another ~33-50% due to Q scores <20
    assert (called["ALT_DP"] > 10).all()
    assert called["ALT_DP"].mean() > 15
