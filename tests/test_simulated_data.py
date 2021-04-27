from Bio import SeqIO
import pandas as pd
import pytest


@pytest.mark.parametrize("n", [x for x in range(1, 2)])
def test_snps_only_fastq(
    tmp_path, n, run_art, run_jobscript, run_snp_mutator, read_vcf_as_dataframe
):
    """Tests insert of N snps"""
    run_snp_mutator(
        input_fasta_file="reference/nCoV-2019.reference.fasta",
        num_subs=n,
        num_insertions=2,
        num_deletions=2,
    )

    # Run ART
    run_art()

    # Run pipeline on simulated data
    run_jobscript(input_filename="simulated_reads.fastq.gz")

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants_table.tsv"), sep="\t")

    # The position, ref, and alt can differ between snpmutator and the vcf if
    # an indel occurs at an ambiguous position (i.e. AGGG -> AGG can give three different
    # variant calls/positions that are equivalent).
    assert len(list(truth["Position"])) == len(list(called["Position (snpmutator)"]))

    # We add these tests to ensure we have a high percent of reads aligning
    # We simulate at 50x, so low end variants with coverage variability should
    # be ~25-30x, and then another ~33-50% due to Q scores <20
    assert (called["Alt depth"] > 10).all()
    assert called["Alt depth"].mean() > 15

    # Finally, test that the FASTAs match
    # Note we ignore the first 50bp which may have low coverage and N masking
    # plus the final 120bps due to a polyA tail
    reference = list(SeqIO.parse(f"{tmp_path}/nCoV-2019.reference_mutated_1.fasta", "fasta"))[0]
    consensus = list(SeqIO.parse(f"{tmp_path}/consensus.fa", "fasta"))[0]
    assert consensus.seq[50:-120] in reference.seq
