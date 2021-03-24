import pandas as pd


def test_twist_truth_data(tmp_path, run_covid_pipeline, read_vcf_as_dataframe):
    """Tests insert of N snps"""
    # Run pipeline on simulated data
    run_covid_pipeline(
        input_filename=(
            "/repo/data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz"
        )
    )
    truth = pd.read_csv(open("data/twist-target-capture/truth.tsv"), sep="\t")
    called = read_vcf_as_dataframe(tmp_path / "variants.vcf")

    # Subset to depth >= 10 and then compare
    # lost last variant because depth went down. differences in filtering/mapping?
    called = called[called["ALT_DP"] >= 10].reset_index()
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])
