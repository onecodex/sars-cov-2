# Helper functions for running the Docker pipeline and snp-mutator
from argparse import Namespace
import os
import subprocess

import pytest
from snpmutator.script import run_from_args


@pytest.fixture
def run_snp_mutator(tmp_path):
    def _run_snp_mutator(**kwargs):
        return run_from_args(_generate_snp_mutator_args(tmp_path, **kwargs))

    return _run_snp_mutator


def _generate_snp_mutator_args(
    tmp_path, input_fasta_file=None, num_subs=10, num_insertions=0, num_deletions=0, random_seed=42
):
    """Wrapper to generate args for SNP mutator"""
    if input_fasta_file is None:
        raise Exception("Requires an input reference FASTA")
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


def _interleave_fastqs(r1, r2, output_filename):
    """Helper function to interleave ART outputs for minimap2"""
    with open(r1) as f1, open(r2) as f2, open(output_filename, "w") as fout:
        while True:
            line = f1.readline()
            if line.strip() == "":
                break
            fout.write(line)

            for _ in range(3):
                fout.write(f1.readline())

            for _ in range(4):
                fout.write(f2.readline())


@pytest.fixture
def run_art(tmp_path):
    def _run_art(
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
                *("--volume", f"{os.getcwd()}:/repo"),
                *("--volume", f"{tmp_path}:/pytest"),
                "--workdir",
                "/pytest",
                "covid19",
                *(
                    *("conda", "run", "-n", "report"),
                    "art_illumina",
                    "--paired",
                    *("--seqSys", system),
                    *("--len", f"{read_length}"),
                    *("--mflen", "200"),
                    *("--sdev", "10"),
                    *("--rndSeed", "1"),
                    *("--fcov", f"{coverage:.2f}"),
                    *("--in", f"/pytest/{input_reference}"),
                    *("--out", "simulated_reads_"),
                ),
            ]
        )
        _interleave_fastqs(
            f"{tmp_path}/simulated_reads_1.fq",
            f"{tmp_path}/simulated_reads_2.fq",
            f"{tmp_path}/simulated_reads.fq",
        )

    return _run_art


def run_docker_container(tmp_path, container_command):
    command = [
        "docker",
        "run",
        "--rm",
        *("--volume", f"{os.getcwd()}:/repo"),
        *("--volume", f"{tmp_path}:/pytest"),
        *("--volume", f"{os.getcwd()}/reference/:/share"),
        *("--workdir", "/pytest"),
        "covid19",
        *container_command,
    ]

    result = subprocess.run(
        command,
        capture_output=True,
    )

    if result.returncode != 0:
        raise Exception(
            "\n".join(
                [
                    "Command Failed!",
                    "command:",
                    " ".join(command),
                    "stdout:",
                    result.stdout.decode("utf-8"),
                    "stderr:",
                    result.stderr.decode("utf-8"),
                ]
            )
        )

    return result


@pytest.fixture
def run_covid_pipeline(tmp_path):
    def _run_covid_pipeline(
        input_filename="nCoV-2019.reference_mutated_1.fasta",
    ):
        container_command = [
            "/bin/bash",
            "/repo/covid19_call_variants.sh",
            "/share/nCoV-2019.reference.fasta",
            input_filename,
            "/share/ARTIC-V1.bed",
        ]

        run_docker_container(tmp_path, container_command)

    return _run_covid_pipeline


@pytest.fixture
def run_artic_covid_pipeline(tmp_path):
    def _run_artic_covid_pipeline(input_filename):
        container_command = [
            "/bin/bash",
            "/repo/covid19_call_variants.artic.sh",
            input_filename,
        ]

        run_docker_container(tmp_path, container_command)

    return _run_artic_covid_pipeline
