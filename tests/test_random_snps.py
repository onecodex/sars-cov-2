from argparse import Namespace
import os
import subprocess

import pandas as pd
import pytest
from snpmutator.script import run_from_args


def _generate_snp_mutator_args(
    input_fasta_file, tmp_path, num_subs=10, num_insertions=0, num_deletions=0, random_seed=None
):
    """Wrapper to generate args for SNP mutator
    """
    kwargs = locals()
    kwargs.pop("tmp_path")
    args = Namespace(
        # Unused but needed defaults
        group_size=None,
        seq_id=None,
        subset_len=0,
        mono=False,
        metrics_file=None,
        concat_ref_file=None,
        # Hard-coded
        fasta_output_dir=tmp_path,
        num_sims=1,
        vcf_file=tmp_path / "variants.vcf",
        summary_file=tmp_path / "summary.tsv",
        # Passed
        **kwargs,
    )
    return args


@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only(tmp_path, n):
    """Tests insert of N snps
    """
    args = _generate_snp_mutator_args("reference/nCoV-2019.reference.fasta", tmp_path, num_subs=n)
    run_from_args(args)

    # Check FASTA files are equal

    # Run pipeline on simulated data
    p = subprocess.call(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{os.getcwd()}:/repo",
            "-v",
            f"{tmp_path}:/pytest",
            "-w",
            "/pytest",
            "covid19",
            "/bin/bash",
            "/repo/covid19_call_variants.sh",
            "/repo/reference/nCoV-2019.reference.fasta",
            "/pytest/nCoV-2019.reference_mutated_1.fasta",
            "/repo/reference/artic-v1/ARTIC-V1.bed",
        ]
    )
    assert p == 0
    assert len(open(tmp_path / "variants.tsv").readlines()) == n + 1

    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])
