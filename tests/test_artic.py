import pandas as pd


def test_sra_illumina_artic(tmp_path, run_covid_pipeline):
    """
    Test that the pipeline generates the same SNPs as detected by NextStrain:
    https://nextstrain.org/ncov?s=USA/CA-PC101P/2020
    """
    run_covid_pipeline(input_filename="/repo/data/ARTIC/SRR11314339.ARTICv1.100k.fastq.gz")
    called = pd.read_csv(tmp_path / "variants.tsv", sep="\t")

    truth = pd.read_csv("data/ARTIC/SRR11314339.ARTICv1.100k.truth.tsv", sep="\t")

    called = called[called["ALT_DP"] >= 10].reset_index()
    # Subset to depth >= 10 and then compare
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])
