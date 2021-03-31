import os


def test_jobscript(tmp_path, run_snp_mutator, run_art, run_jobscript):
    # generate some data
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=3)
    run_art(coverage=15)
    run_jobscript(input_filename="simulated_reads.fastq.gz")

    assert os.path.exists(tmp_path / "report.pdf")
