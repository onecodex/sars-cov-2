FROM continuumio/miniconda3:latest

LABEL authors="nick.greenfield@invitae.com,austin.richardson@invitae.com,christine.he@invitae.com"
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
  texlive-xetex \
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
  && curl -sL https://deb.nodesource.com/setup_13.x | bash - \
  && apt-get install -y nodejs \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# install artic into conda environment "artic"
RUN git clone https://github.com/artic-network/fieldbioinformatics.git \
        && cd fieldbioinformatics \
        && /opt/conda/bin/conda env create -f environment.yml \
        && /opt/conda/bin/conda run -n artic python setup.py install \
        && /opt/conda/bin/conda clean -a

## install pangolin into conda environment "pangolin"
RUN git clone https://github.com/cov-lineages/pangolin.git \
        && cd pangolin \
        && /opt/conda/bin/conda env create -f environment.yml \
        && /opt/conda/bin/conda run -n pangolin python setup.py install \
        && /opt/conda/bin/conda clean -a


# Install nextclade CLI
RUN npm install --global @neherlab/nextclade

# install main report environment
COPY environment.yml /
RUN conda env create -f environment.yml

# Install snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
        && unzip snpEff_latest_core.zip \
        && rm snpEff_latest_core.zip

# Report generation
RUN npm install -g --unsafe-perm vega-lite vega-cli canvas

ADD covid19_call_variants.sh /usr/local/bin/
ADD covid19_call_variants.artic.sh /usr/local/bin/
ADD post_process_variants.sh /usr/local/bin/
ADD generate_tsv.py /usr/local/bin

ADD report.ipynb /
ADD nCoV-2019.reference.fasta /
ADD nCoV-2019.reference.gtf /
ADD annot_table.orfs.txt /
ADD ivar_variants_to_vcf.py /
