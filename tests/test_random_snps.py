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


def run_art(
    tmp_path,
    input_reference="nCoV-2019.reference_mutated_1.fasta",
    system="MSv3",
    read_length=150,
    coverage=50,
):
    return subprocess.call(
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
            "art_illumina",
            "--seqSys",
            system,
            "--len",
            f"{read_length}",
            "--mflen",
            "200",
            "--sdev",
            "10",
            "--paired",
            "--rndSeed",
            "1",
            "--fcov",
            f"{coverage:.2f}",
            "--in",
            f"/pytest/{input_reference}",
            "--out",
            "simulated_reads_",
        ]
    )


def run_covid_pipeline(tmp_path, input_filename="nCoV-2019.reference_mutated_1.fasta"):
    return subprocess.call(
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
            f"/pytest/{input_filename}",
            "/repo/reference/artic-v1/ARTIC-V1.bed",
        ]
    )


@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only_fasta(tmp_path, n):
    """Tests insert of N snps
    """
    args = _generate_snp_mutator_args("reference/nCoV-2019.reference.fasta", tmp_path, num_subs=n)
    run_from_args(args)

    # Run pipeline on simulated data
    p = run_covid_pipeline(tmp_path)
    assert p == 0

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])


@pytest.mark.xfail(strict=False)  # some of these fail at the moment, but *should not*
@pytest.mark.parametrize("n", [x for x in range(1, 10)])
def test_snps_only_fastq(tmp_path, n):
    """Tests insert of N snps
    """
    args = _generate_snp_mutator_args("reference/nCoV-2019.reference.fasta", tmp_path, num_subs=n)
    run_from_args(args)

    # Run ART
    p = run_art(tmp_path)
    assert p == 0

    # Run pipeline on simulated data
    p = run_covid_pipeline(tmp_path, input_filename="simulated_reads_1.fq")
    assert p == 0

    # Check that all variants are detected and there are no extras
    truth = pd.read_csv(open(tmp_path / "summary.tsv"), sep="\t")
    called = pd.read_csv(open(tmp_path / "variants.tsv"), sep="\t")
    assert all(truth["Position"] == called["POS"])
    assert all(truth["OriginalBase"] == called["REF"])
    assert all(truth["NewBase"] == called["ALT"])
