# SARS-CoV-2
SARS-CoV-2 variant calling and consensus assembly pipeline

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
docker run --rm -w /data -v `pwd`:/data covid19 bash covid19_call_variants.sh \
    reference/nCoV-2019.reference.fasta \
    data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz \
    reference/artic-v1/ARTIC-V1.bed
```

This currently produces a `consensus.fa` file and a `variants.tsv`.

# Testing
This repository includes a local `requirements.txt` file for quickly running some golden output tests across a variety of datasets. This repository is set up to use Github Actions to automatically build the Docker image and run those tests to ensure there are no regressions.

This is to ensure that parameter and pipeline changes don't affect variant calls or consensus sequence generation.

# Acknowledgments
Many thanks are due across the community, including _but not limited_ to:
- X, Y, and Z for sharing data
- @torstenseeman, @gkarthik, @nickloman, and many others for quick discussions on optimal SNP calling for both amplicon (ARTIC primers) and non-amplicon sequencing approaches
- @nickloman, @joshquick, @rambaut, @k-florek and others working on the [ARTIC protocol for SARS-CoV-2](https://github.com/artic-network/artic-ncov2019)