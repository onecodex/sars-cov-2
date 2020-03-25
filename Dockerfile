FROM continuumio/miniconda3:latest

LABEL authors="nick@onecodex.com"
LABEL description="Docker image for SARS-CoV-2 analysis"
LABEL software.version="0.1.0"

RUN apt-get update && apt-get install -y g++ git make procps && apt-get clean -y

COPY environment.yml /

RUN /opt/conda/bin/conda env create -f /environment.yml && /opt/conda/bin/conda clean -a

# Needed for report generation
RUN apt-get install -yq \
    fonts-dejavu \
    fonts-texgyre \
    texlive-fonts-recommended \
    texlive-generic-recommended \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD covid19_call_variants.sh /usr/local/bin/
ADD report.ipynb .

ENV PATH /opt/conda/envs/covid-19/bin:$PATH
