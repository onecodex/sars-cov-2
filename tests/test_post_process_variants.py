import os


def test_post_process_variants(tmp_path, run_post_process_variants, run_snp_mutator):
    run_snp_mutator(input_fasta_file="reference/nCoV-2019.reference.fasta", num_subs=1)
    run_post_process_variants("nCoV-2019.reference_mutated_1.fasta")
    assert os.path.exists(tmp_path / "nextclade.tsv")
    assert os.path.exists(tmp_path / "nextclade.json")
    assert os.path.exists(tmp_path / "pangolin.csv")
