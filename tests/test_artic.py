import pytest

import pandas as pd
import os


def test_sra_illumina_artic(tmp_path, run_covid_pipeline):
    """
    Test that the pipeline generates the same SNPs as detected by NextStrain:
    https://nextstrain.org/ncov?s=USA/CA-PC101P/2020
    """

    run_covid_pipeline(input_filename="/repo/data/ARTIC/SRR11314339.ARTICv1.100k.fastq.gz")

    assert os.path.exists(tmp_path / "consensus.fa")
    assert os.path.exists(tmp_path / "variants.vcf")
    assert os.path.exists(tmp_path / "covid19.bam")
    assert os.path.exists(tmp_path / "covid19.bam.fa")

    # for backwards compatibility
    # need to rewrite notebook to not depend on this TSV (use nextclade output instead)
    called = pd.read_csv(tmp_path / "variants.tsv", sep="\t")

    truth = pd.read_csv("data/ARTIC/SRR11314339.ARTICv1.100k.truth.tsv", sep="\t")

    called = called[called["ALT_DP"] >= 10].reset_index()
    # Subset to depth >= 10 and then compare
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])


def test_ont_artic(tmp_path, run_artic_covid_pipeline):
    run_artic_covid_pipeline(input_filename="/repo/data/ARTIC/RT-1212.barcode04.21_02_25.fastq.gz")

    assert os.path.exists(tmp_path / "consensus.fa")
    assert os.path.exists(tmp_path / "variants.vcf")
    assert os.path.exists(tmp_path / "covid19.bam")
    assert os.path.exists(tmp_path / "covid19.bam.fa")


def test_post_process_variants(tmp_path, run_post_process_variants):
    run_post_process_variants("/repo/data/test-consensus.fa")
    assert os.path.exists(tmp_path / "nextclade.tsv")
    assert os.path.exists(tmp_path / "nextclade.json")
    assert os.path.exists(tmp_path / "pangolin.csv")
