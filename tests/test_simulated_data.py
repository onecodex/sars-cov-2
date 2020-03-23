from conftest import run_art
from conftest import run_covid_pipeline
from conftest import run_snp_mutator
import pandas as pd
import pytest


@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only_fasta(tmp_path, n):
    """Tests insert of N snps
    """
    run_snp_mutator("reference/nCoV-2019.reference.fasta", tmp_path, num_subs=n)

    # Run pipeline on simulated data
    run_covid_pipeline(tmp_path)

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])


# @pytest.mark.xfail(strict=False)  # some of these fail at the moment, but *should not*
@pytest.mark.parametrize("n", [x for x in range(1, 10)])
# @pytest.mark.parametrize("n", [x for x in [10, 20, 30]])
# @pytest.mark.parametrize("n", [x for x in [30]])
def test_snps_only_fastq(tmp_path, n):
    """Tests insert of N snps
    """
    run_snp_mutator("reference/nCoV-2019.reference.fasta", tmp_path, num_subs=n)

    # Run ART
    run_art(tmp_path)

    # Run pipeline on simulated data
    run_covid_pipeline(tmp_path, input_filename="simulated_reads_1.fq")

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])
