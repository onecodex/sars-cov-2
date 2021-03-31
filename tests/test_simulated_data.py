from Bio import SeqIO
import pandas as pd
import pytest


@pytest.mark.parametrize("n", [x for x in range(1, 2)])
def test_snps_only_fastq(
    tmp_path, n, run_art, run_call_variants_illumina, run_snp_mutator, read_vcf_as_dataframe
):
    """Tests insert of N snps"""
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=n)

    # Run ART
    run_art()

    # Run pipeline on simulated data
    run_call_variants_illumina(input_filename="simulated_reads.fastq.gz")

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = read_vcf_as_dataframe(tmp_path / "variants.vcf")

    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])

    # We add these tests to ensure we have a high percent of reads aligning
    # We simulate at 50x, so low end variants with coverage variability should
    # be ~25-30x, and then another ~33-50% due to Q scores <20
    assert (called["ALT_DP"] > 10).all()
    assert called["ALT_DP"].mean() > 15

    # Finally, test that the FASTAs match
    # Note we ignore the first 50bp which may have low coverage and N masking
    # plus the final 120bps due to a polyA tail
    reference = list(SeqIO.parse(f"{tmp_path}/nCoV-2019.reference_mutated_1.fasta", "fasta"))[0]
    consensus = list(SeqIO.parse(f"{tmp_path}/consensus.fa", "fasta"))[0]
    assert consensus.seq[50:-120] in reference.seq
