import os

import pandas as pd


def test_sra_illumina_artic(tmp_path, run_call_variants_illumina, read_vcf_as_dataframe):
    """
    Test that the pipeline generates the same SNPs as detected by NextStrain:
    https://nextstrain.org/ncov?s=USA/CA-PC101P/2020
    """

    run_call_variants_illumina(input_filename="/repo/data/ARTIC/SRR11314339.ARTICv1.100k.fastq.gz")

    assert os.path.exists(tmp_path / "consensus.fa")
    assert os.path.exists(tmp_path / "variants.vcf")
    assert os.path.exists(tmp_path / "covid19.bam")
    assert os.path.exists(tmp_path / "covid19.bam.bai")

    # for backwards compatibility
    # need to rewrite notebook to not depend on this TSV (use nextclade output instead)
    called = read_vcf_as_dataframe(tmp_path / "variants.vcf")

    truth = pd.read_csv(
        "data/ARTIC/SRR11314339.ARTICv1.100k.truth.tsv",
        sep="\t",
        dtype={
            "Replicate": "str",
            "Position": "str",
            "OriginalBase": "str",
            "NewBase": "str",
        },
    )

    # Subset to depth >= 10 and then compare
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])


def test_ont_artic(tmp_path, run_call_variants_ont, read_vcf_as_dataframe):
    run_call_variants_ont(input_filename="/repo/data/ARTIC/ERR5284916.ONT.ARTICv3.40k.fastq.gz")
    assert os.path.exists(tmp_path / "consensus.fa")
    assert os.path.exists(tmp_path / "variants.vcf")
    assert os.path.exists(tmp_path / "covid19.bam")
    assert os.path.exists(tmp_path / "covid19.bam.bai")
