import pandas as pd


def test_twist_truth_data(tmp_path, run_call_variants_illumina, read_vcf_as_dataframe):
    """Tests insert of N snps"""
    # Run pipeline on simulated data
    run_call_variants_illumina(
        input_filename=(
            "/repo/data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz"
        )
    )
    truth = pd.read_csv(
        open("data/twist-target-capture/truth.tsv"),
        sep="\t",
        dtype={
            "Replicate": "str",
            "Position": "str",
            "OriginalBase": "str",
            "NewBase": "str",
        },
    )
    called = read_vcf_as_dataframe(tmp_path / "variants.vcf")
    # Subset to depth >= 10 and then compare
    # lost last variant because depth went down. differences in filtering/mapping?
    called = called[called["ALT_DP"] >= 10]
    assert all(truth["Position"].values == called["POS"].values)
    assert all(truth["OriginalBase"].values == called["REF"].values)
    assert all(truth["NewBase"].values == called["ALT"].values)
