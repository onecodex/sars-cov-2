# SARS-CoV-2

Version v0.4.10

[![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/test/badge.svg)](https://github.com/onecodex/sars-cov-2/actions) [![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/pre-commit/badge.svg)](https://github.com/onecodex/sars-cov-2/actions) [![Docker Repository on Quay](https://quay.io/repository/refgenomics/covid19/status "Docker Repository on Quay")](https://quay.io/repository/refgenomics/covid19)

SARS-CoV-2 variant calling and consensus assembly pipeline for ARTIC v3 amplicons sequenced on Illumina or Oxford Nanopore platforms

# Quick start

```sh
docker build -t covid19 .
```

Run the pipeline in the Docker image:

```sh
docker \
  run \
  --rm \
  --workdir /data \
  --volume `pwd`:/data \
  --entrypoint /bin/bash \
  --env prefix=test-covid19 \
  --env reference=reference/nCoV-2019.reference.fasta \
  --env input_fastq=data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz \
  --env primer_bed_file=reference/artic-v1/ARTIC-V3.bed \
  covid19 \
  jobscript.sh
```

For Oxford Nanopore:

```sh
docker \
  run \
  --rm \
  --workdir /data \
  --volume `pwd`:/data \
  --entrypoint /bin/bash \
  --env prefix=test-covid19 \
  --env INSTRUMENT_VENDOR="Oxford Nanopore" \
  --env reference=reference/nCoV-2019.reference.fasta \
  --env input_fastq=data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz \
  --env primer_bed_file=reference/artic-v1/ARTIC-V3.bed \
  covid19 \
  jobscript.sh
```

This currently produces a `consensus.fa` file, a `variants.vcf`, a BAM file (`covid19.bam`), nextstrain results (`nextstrain.json`) and pangolin results (`pangolin.csv`).

# Development & Testing

To run tests, run `pytest`.

This repository includes a local `requirements.txt` file for quickly running some golden output tests across a variety of datasets. This repository is set up to use Github Actions to automatically build the Docker image and run those tests to ensure there are no regressions. These ensure that parameter and pipeline changes don't affect variant calls or consensus sequence generation.

Currently, the following integration tests are run:
* Simulated Illumina data from the SARS-CoV-2 reference including _simulated variants across the genome_
* Example Twist hybrid capture data (Illumina)
* Example ARTIC v1 amplicon sequencing data (Illumina)

It also uses [`pre-commit`](https://pre-commit.com/) to keep things clean and orderly. To get started, first install the requirements (Python 3 required): `pip install -r requirements.txt`. Then install the `pre-commit` hooks: `pre-commit install --install-hooks`. Note that you'll also need [`shellcheck`](https://www.shellcheck.net/) installed on your system (`brew install shellcheck` on a Mac).

# Acknowledgments

Many thanks are due across the community, including _but not limited_ to:
- [@tseemann](https://github.com/tseemann), [@gkarthik](https://github.com/gkarthik), [@nickloman](https://github.com/nickloman), and many others for quick discussions on optimal SNP calling for both amplicon (ARTIC primers) and non-amplicon sequencing approaches
- [@nickloman](https://github.com/nickloman), [@joshquick](https://github.com/joshquick), [@rambaut](https://github.com/rambaut), [@k-florek](https://github.com/k-florek) and others working on the [ARTIC protocol for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019)
- [@pangolin](https://github.com/cov-lineages/pangolin) and [@nextclade](https://github.com/nextstrain/nextclade) for surveillance tools
- [Voigt lab](http://web.mit.edu/voigtlab/) for [dnaplotlib](https://github.com/VoigtLab/dnaplotlib)
