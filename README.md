# SARS-CoV-2 variant calling

[![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/test/badge.svg)](https://github.com/onecodex/sars-cov-2/actions) [![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/pre-commit/badge.svg)](https://github.com/onecodex/sars-cov-2/actions) [![Docker Repository on Quay](https://quay.io/repository/refgenomics/covid19/status "Docker Repository on Quay")](https://quay.io/repository/refgenomics/covid19)

This pipeline performs consensus assembly and variant calling for amplicon sequencing data (Illumina or Oxford Nanopore) generated using the [`ARTIC protocol`](https://artic.network/ncov-2019). The user can specify a primer set to be used for trimming alignments; this is assumed to be **ARTIC V4.1** if not specified.

# Pipeline overview

The pipeline takes in a single FASTQ file (interleaved if Illumina) and processes it as follows:

1. Map reads to the Wuhan-Hu-1 reference and trim ARTIC primer sequences
2. Generate a consensus sequence (bcftools for Illumina; medaka for Oxford Nanopore)
3. Call variants with a limit of 200x coverage, as recommended by the ARTIC network. While indels and SNVs are reported for Illumina data, only SNVs are reported for Oxford Nanopore based on benchmarking studies that indicate [`small indel detection is unreliable.`](https://doi.org/10.1038/s41467-020-20075-6)
4. Assign Pangolin and Nextclade lineages
5. Predict amino acid mutations
6. Predict consequences of compound variants (ex: adjacent SNVs on the same codon; frame-shifting indels followed by frame-restoring indels)

Rigorous quality checks are implemented throughout the pipeline, including flagging of variants in low complexity regions for error-prone Oxford Nanopore data and conservative lineage calls (no lineage assignments will be reported if the consensus sequence has too many N’s or is too fragmented).

In addition to a results JSON, a PDF report is generated for each sample that tells you at a glance whether primer dropout has occurred, which amino acid mutations are present, and whether the sample contains a variant of concern. An example report is shown below.

![Example report](https://www.onecodex.com/uploads/sars-cov-2-report-2021-example.png)

# Quick start

```sh
docker build -t covid19 .
```

Run the pipeline in the Docker image (note that fastq files are stored in git lfs so you may need to `git lfs pull` before executing):

```sh
docker \
  run \
  --rm \
  --workdir /data \
  --volume `pwd`:/data \
  --entrypoint /bin/bash \
  --env ONE_CODEX_REPORT_FILENAME=report.pdf \
  --env INSTRUMENT_VENDOR=Illumina \
  --env ARTIC_PRIMER_VERSION=4.1 \
  covid19 \
  jobscript.sh \
  data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz
```

For Oxford Nanopore:

```sh
docker \
  run \
  --rm \
  --workdir /data \
  --volume `pwd`:/data \
  --entrypoint /bin/bash \
  --env ONE_CODEX_REPORT_FILENAME=report.pdf \
  --env INSTRUMENT_VENDOR="Oxford Nanopore" \
  --env ARTIC_PRIMER_VERSION=4.1 \
  covid19 \
  jobscript.sh \
  data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz
```

# Development & Testing

To run tests, run `pytest`.

The `requirements.txt` file lists dependencies for quickly running some golden output tests across a variety of datasets. This repository is set up to use Github Actions to automatically build the Docker image and run these tests, to ensure that parameter and pipeline changes don't affect variant calls or consensus sequence generation.

Currently, integration tests are run on:
* Simulated Illumina data from the SARS-CoV-2 reference including _simulated variants across the genome_
* Example Twist hybrid capture data (Illumina)
* Example ARTIC v1 amplicon sequencing data (Illumina)

It also uses [`pre-commit`](https://pre-commit.com/) to keep things clean and orderly. To get started, first install the requirements (Python 3 required): `pip install -r requirements.txt`. Then install the `pre-commit` hooks: `pre-commit install --install-hooks`.

# Acknowledgments

Many thanks are due across the community, including _but not limited_ to:
- [@tseemann](https://github.com/tseemann), [@gkarthik](https://github.com/gkarthik), [@nickloman](https://github.com/nickloman), and many others for quick discussions on optimal SNP calling for both amplicon (ARTIC primers) and non-amplicon sequencing approaches
- [@nickloman](https://github.com/nickloman), [@joshquick](https://github.com/joshquick), [@rambaut](https://github.com/rambaut), [@k-florek](https://github.com/k-florek) and others working on the [ARTIC protocol for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019)
- [@pangolin](https://github.com/cov-lineages/pangolin) and [@nextclade](https://github.com/nextstrain/nextclade) for surveillance tools
- [Voigt lab](http://web.mit.edu/voigtlab/) for [dnaplotlib](https://github.com/VoigtLab/dnaplotlib)
