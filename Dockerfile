FROM continuumio/miniconda3:latest

LABEL authors="nick@onecodex.com"
LABEL description="Docker image for SARS-CoV-2 analysis"
LABEL software.version="0.1.0"

RUN apt-get update && apt-get install -y g++ git make procps && apt-get clean -y

# Needed for report generation
RUN apt-get install -yq \
    fonts-dejavu \
    fonts-texgyre \
    texlive-fonts-recommended \
    texlive-generic-recommended \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /

RUN /opt/conda/bin/conda env create -f /environment.yml && /opt/conda/bin/conda clean -a
ENV PATH /opt/conda/envs/covid-19/bin:$PATH

# Altair rendering requirements
USER root
RUN apt-get update \
    && apt-get install -y gnupg \
    && curl -sL https://deb.nodesource.com/setup_13.x  | bash - \
    && apt-get install -y nodejs \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN npm install -g --unsafe-perm vega-lite vega-cli canvas

RUN pip install git+https://github.com/onecodex/altair_saver@a4be53dad5df23d68f37f2bb9f0782ab79ef7c96#egg=altair_saver git+https://github.com/onecodex/onecodex@3adfb95bc639435d41c9e5a337a9887e82024cb0#egg=onecodex[all,reports]

ADD covid19_call_variants.sh /usr/local/bin/
ADD report.ipynb .
