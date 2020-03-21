FROM continuumio/miniconda3:latest

LABEL authors="nick@onecodex.com"
LABEL description="Docker image for SARS-CoV-2 analysis"
LABEL software.version="0.1.0"

COPY environment.yml /
RUN apt-get update && apt-get install -y g++ git make procps && apt-get clean -y
RUN /opt/conda/bin/conda env create -f /environment.yml && /opt/conda/bin/conda clean -a
ENV PATH /opt/conda/envs/covid-19/bin:$PATH
