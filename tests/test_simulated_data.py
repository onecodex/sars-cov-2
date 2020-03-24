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
    run_covid_pipeline(input_filename="simulated_reads_1.fq")

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])

    # We simulate at 50x, so all calls should be >= 10x
    assert (called["ALT_DP"] > 10).all()
