import os

import pytest


def test_jobscript_ont(tmp_path, run_jobscript, n_variants):
    run_jobscript(
        input_filename="/repo/data/ARTIC/ERR5284916.ONT.ARTICv3.40k.fastq.gz",
        instrument_vendor="Oxford Nanopore",
    )

    assert os.path.exists(tmp_path / "report.pdf")


def test_jobscript_ont_no_variants(tmp_path, run_jobscript, n_variants):
    run_jobscript(
        input_filename="/repo/data/ARTIC/ERR5284916.ONT.ARTICv3.40k.fastq.gz",
        instrument_vendor="Oxford Nanopore",
    )

    # simulate having no variants detected
    with open(tmp_path / "variants.vcf", "wt") as handle:
        handle.write("# no variants\n")

    assert os.path.exists(tmp_path / "report.pdf")


@pytest.mark.parametrize("n_variants", [0, 3])
def test_jobscript_illumina(tmp_path, run_snp_mutator, run_art, run_jobscript, n_variants):
    # generate some data
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=n_variants)
    run_art(coverage=15)
    run_jobscript(input_filename="simulated_reads.fastq.gz")

    assert os.path.exists(tmp_path / "report.pdf")
