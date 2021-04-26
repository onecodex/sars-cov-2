# Helper functions for running the Docker pipeline and snp-mutator
from argparse import Namespace
import gzip
import os
import shutil
import subprocess

import pandas as pd
import pytest
from snpmutator.script import run_from_args


@pytest.fixture
def read_vcf_as_dataframe():
    def _read_vcf_as_dataframe(path):
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
            usecols=[0, 1, 2, 3, 4, 5, 6, 7],
        )

        # parse info add key/value items to data-frame
        info_data = []
        for row in df["INFO"].values:
            items = [i.split("=") for i in row.split(";")]
            parsed = {}
            for item in items:
                # boolean items. skip for now
                if len(item) == 1:
                    pass
                elif len(item) == 2:
                    parsed[item[0]] = item[1]
                else:
                    raise Exception(f"Item of unexpected length in vcf info column: {item}")
            # parse depth column into ALT and REF depths
            depths = [int(x) for x in parsed["DP4"].split(",")]
            parsed["ALT_DP"] = sum(depths[2:])
            parsed["REF_DP"] = sum(depths[:2])
            info_data.append(parsed)

        df = pd.DataFrame([{**r, **i} for r, i in zip(df.to_dict(orient="records"), info_data)])

        return df

    return _read_vcf_as_dataframe


@pytest.fixture
def run_snp_mutator(tmp_path):
    def _run_snp_mutator(**kwargs):
        return run_from_args(_generate_snp_mutator_args(tmp_path, **kwargs))

    return _run_snp_mutator


def _generate_snp_mutator_args(
    tmp_path,
    input_fasta_file=None,
    num_subs=10,
    num_insertions=0,
    num_deletions=0,
    random_seed=42,
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
    with open(r1) as f1, open(r2) as f2, gzip.open(output_filename, "wt") as fout:
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
            f"{tmp_path}/simulated_reads.fastq.gz",
        )

    return _run_art


def run_docker_container(tmp_path, container_command, env=None):
    if env is None:
        env = {}

    command = [
        "docker",
        "run",
        "--rm",
        *("--volume", f"{os.getcwd()}:/repo"),
        *("--volume", f"{tmp_path}:/pytest"),
        *("--volume", f"{os.getcwd()}/reference/:/share"),
        *("--workdir", "/pytest"),
        *[i for r in [("--env", f"{k}={v}") for k, v in env.items()] for i in r],
        "covid19",
        *container_command,
    ]

    result = subprocess.run(command, capture_output=True)

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
def run_post_process_variants(tmp_path):
    def _run_post_process_variants(input_filename=None):
        container_command = [
            "/bin/bash",
            "/repo/post_process_variants.sh",
            input_filename,
        ]

        run_docker_container(tmp_path, container_command)

    return _run_post_process_variants


@pytest.fixture
def run_jobscript(tmp_path):
    def _run_jobscript(input_filename, instrument_vendor="Illumina"):
        container_command = ["/bin/bash", "/repo/jobscript.sh", input_filename]

        run_docker_container(
            tmp_path, container_command, env={"INSTRUMENT_VENDOR": instrument_vendor}
        )

    return _run_jobscript


@pytest.fixture
def run_call_variants_illumina(tmp_path):
    def _run_call_variants_illumina(
        input_filename="nCoV-2019.reference_mutated_1.fasta",
    ):
        container_command = [
            "/bin/bash",
            "/repo/covid19_call_variants.sh",
            "/share/nCoV-2019.reference.fasta",
            input_filename,
            "/share/ARTIC-V3.bed",
        ]

        run_docker_container(tmp_path, container_command)

    return _run_call_variants_illumina


@pytest.fixture
def run_call_variants_ont(tmp_path):
    def _run_call_variants_ont(input_filename):
        container_command = [
            "/bin/bash",
            "/repo/covid19_call_variants.ont.sh",
            input_filename,
        ]

        run_docker_container(tmp_path, container_command)

    return _run_call_variants_ont


@pytest.fixture
def run_render_notebook(tmp_path):
    def _run_render_notebook():
        # notebook must exist in tmp_path because file paths are relative to
        # the notebook
        shutil.copy("report.ipynb", tmp_path / "report.ipynb")
        container_command = [
            *("conda", "run", "-n", "report"),
            "jupyter",
            "nbconvert",
            "--execute",
            "--to=onecodex_pdf",
            "--ExecutePreprocessor.timeout=-1",
            "--output='test.pdf'",
            "--output-dir='.'",
            "report.ipynb",
        ]
        run_docker_container(tmp_path, container_command)

    return _run_render_notebook
