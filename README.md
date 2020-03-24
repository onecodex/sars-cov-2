# SARS-CoV-2

[![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/test/badge.svg)](https://github.com/onecodex/sars-cov-2/actions) [![Actions Status](https://github.com/onecodex/sars-cov-2/workflows/pre-commit/badge.svg)](https://github.com/onecodex/sars-cov-2/actions)

SARS-CoV-2 variant calling and consensus assembly pipeline using [ivar](https://github.com/andersen-lab/ivar) and [minimap2](https://github.com/lh3/minimap2).

⚠️ **As of March 24th, this pipeline support ARTIC v1 primers, capture-based approaches, and RNA-seq _on Illumina platforms_. We expect to add ONT support in the coming days, but recommend looking at the [ARTIC Network protocol](https://artic.network/ncov-2019) in the meantime.**

# Quick start

This repo includes a number of compressed FASTQ and other files for use as "golden output" test data. In order to access this data, you'll need to have `git lfs` installed (`git-lfs` can be installed via `brew install git-lfs` on Mac OS):

```sh
git lfs install
```

You can then build the Docker image as follows:

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
  --env primer_bed_file=reference/artic-v1/ARTIC-V1.bed \
  covid19 \
  covid19_call_variants.sh
```

This currently produces a `consensus.fa` file and a `variants.tsv`.

# Development & Testing

This repository includes a local `requirements.txt` file for quickly running some golden output tests across a variety of datasets. This repository is set up to use Github Actions to automatically build the Docker image and run those tests to ensure there are no regressions. These ensure that parameter and pipeline changes don't affect variant calls or consensus sequence generation.

It also uses [`pre-commit`](https://pre-commit.com/) to keep things clean and orderly. To get started, first install the requirements (Python 3 required): `pip install -r requirements.txt`. Then install the `pre-commit` hooks: `pre-commit install --install-hooks`. Note that you'll also need [`shellcheck`](https://www.shellcheck.net/) installed on your system (`brew install shellcheck` on a Mac).

# Acknowledgments

Many thanks are due across the community, including _but not limited_ to:
- [@tseemann](https://github.com/tseemann), [@gkarthik](https://github.com/gkarthik), [@nickloman](https://github.com/nickloman), and many others for quick discussions on optimal SNP calling for both amplicon (ARTIC primers) and non-amplicon sequencing approaches
- [@nickloman](https://github.com/nickloman), [@joshquick](https://github.com/joshquick), [@rambaut](https://github.com/rambaut), [@k-florek](https://github.com/k-florek) and others working on the [ARTIC protocol for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019)
