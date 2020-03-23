# Helper functions for running the Docker pipeline and snp-mutator
# TODO: Refactor into fixtures?
from argparse import Namespace
import os
import subprocess

from snpmutator.script import run_from_args


def run_snp_mutator(*args, **kwargs):
    return run_from_args(_generate_snp_mutator_args(*args, **kwargs))


def _generate_snp_mutator_args(
    input_fasta_file, tmp_path, num_subs=10, num_insertions=0, num_deletions=0, random_seed=42
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
    subprocess.check_output(
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
    subprocess.call(
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
            input_filename,
            "/repo/reference/artic-v1/ARTIC-V1.bed",
        ]
    )
