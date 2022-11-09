# these tests will copy the reports to the current working directory
# there is a GitHub action setup to upload these reports as artifacts

import os
import shutil

import pytest



def test_jobscript_ont(tmp_path, run_jobscript):
    run_jobscript(
        input_filename="/repo/data/ARTIC/ERR5284916.ONT.ARTICv3.40k.fastq.gz",
        instrument_vendor="Oxford Nanopore",
    )

    assert os.path.exists(tmp_path / "report.pdf")

    shutil.copy(tmp_path / "report.pdf", f"report-ont.pdf-{num_subs}_subs.pdf")

@pytest.mark.parametrize("num_subs", [0, 3])
def test_jobscript_illumina(tmp_path, run_snp_mutator, run_art, run_jobscript, num_subs):
    # generate some data
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=num_subs)
    run_art(coverage=15)
    run_jobscript(input_filename="simulated_reads.fastq.gz")

    assert os.path.exists(tmp_path / "report.pdf")

    shutil.copy(tmp_path / "report.pdf", f"report-illumina-{num_subs}_subs.pdf")
