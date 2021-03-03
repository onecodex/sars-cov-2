FROM continuumio/miniconda3:latest

LABEL authors="nick@onecodex.com"
LABEL description="Docker image for SARS-CoV-2 analysis"
LABEL software.version="0.1.0"

RUN apt-get update && apt-get install -y g++ \
	git \
	make \
	procps \
	curl \
	build-essential \
	unzip \
	&& apt-get clean -y

# Needed for report generation
RUN apt-get install -yq \
	fonts-dejavu \
	fonts-texgyre \
	texlive-fonts-recommended \
	texlive-generic-recommended \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Altair rendering requirements
USER root
RUN apt-get update \
	&& apt-get install -y gnupg \
	&& curl -sL https://deb.nodesource.com/setup_13.x  | bash - \
	&& apt-get install -y nodejs \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Create ARTIC environment
RUN git clone https://github.com/artic-network/artic-ncov2019.git \
        && cd artic-ncov2019 \
        && /opt/conda/bin/conda env create -f environment.yml \
        && /opt/conda/bin/conda clean -a

# Install snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
        && unzip snpEff_latest_core.zip \
        && rm snpEff_latest_core.zip

ENV PATH /opt/conda/envs/covid-19/bin:$PATH

# Install nextclade CLI
RUN npm install --global @neherlab/nextclade

# Create post-ARTIC environment (Pangolin plus additional dependencies; still called pangolin)
COPY post-artic.yml .
RUN git clone https://github.com/cov-lineages/pangolin.git \
        && cd pangolin \
        && mv ../post-artic.yml environment.yml \
        && /opt/conda/bin/conda env create -f environment.yml \
        && /opt/conda/bin/conda run -n pangolin python setup.py install \
        && /opt/conda/bin/conda clean -a

# Report generation
RUN npm install -g --unsafe-perm vega-lite vega-cli canvas

RUN pip install \
  git+https://github.com/onecodex/altair_saver@a4be53dad5df23d68f37f2bb9f0782ab79ef7c96#egg=altair_saver \
  git+https://github.com/onecodex/onecodex@74c502956a033c335363ccbf5d461791814adad8#egg=onecodex[all,reports]

ADD covid19_call_variants.sh /usr/local/bin/
ADD covid19_call_variants.artic.sh .
ADD report.ipynb .
ADD nCoV-2019.reference.fasta .
ADD nCoV-2019.reference.gtf .
ADD aa_codes.txt .
ADD annot_table.orfs.txt .
ADD generate_tsv.py .
