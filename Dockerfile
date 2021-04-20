FROM python:3.8

# dependencies for generating report
ENV LD_LIBRARY_PATH=/usr/local/lib

RUN pip install numpy

RUN pip install pysam==0.16 biopython==1.78 PyVCF


USER root
RUN apt-get update \
    && apt-get autoclean \
    && apt-get install -y gnupg \
    && curl -sL https://deb.nodesource.com/setup_14.x  | bash - \
    && apt-get install -y nodejs \
    unzip \
    default-jre \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN npm install -g --unsafe-perm vega vega-lite vega-cli canvas

RUN apt-get update && apt-get install -y curl

# install Conda
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  > Miniconda3-latest-Linux-x86_64.sh \
  && yes \
  | bash Miniconda3-latest-Linux-x86_64.sh -b

# put system path first so that conda doesn't override python
ENV PATH=$PATH:/root/miniconda3/bin/

# install "report" environment's dependencies
COPY environment.yml /
RUN conda env create -f environment.yml

# install artic into conda environment "artic"
RUN git clone https://github.com/artic-network/fieldbioinformatics.git \
        && cd fieldbioinformatics \
        && conda env create -f environment.yml \
        && conda run -n artic python setup.py install \
        && conda clean -a

# install pangolin into conda environment "pangolin"
RUN git clone https://github.com/cov-lineages/pangolin.git \
        && cd pangolin \
        && conda env create -f environment.yml \
        && conda run -n pangolin python setup.py install \
        && conda clean -a

# install nextclade
RUN npm install --global @neherlab/nextclade

# update pangolin database 2021-03-18
RUN conda run -n pangolin pangolin --update

# install dnaplotlib for creating the genome diagram
RUN pip install dnaplotlib

RUN pip install nbconvert==6.0.7

# Install onecodex_pdf export option
RUN pip install onecodex[all,reports]==v0.9.4
RUN mkdir -p /usr/local/share/fonts \
    && cp /usr/local/lib/python3.8/site-packages/onecodex/assets/fonts/*.otf /usr/local/share/fonts \
    && fc-cache

# install snpeff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
        && unzip snpEff_latest_core.zip \
	&& mv snpEff /usr/local/bin \
        && rm snpEff_latest_core.zip

ADD covid19_call_variants.sh /usr/local/bin/
ADD covid19_call_variants.artic.sh /usr/local/bin/
ADD post_process_variants.sh /usr/local/bin/
ADD jobscript.sh /usr/local/bin/
ADD generate_tsv.py /usr/local/bin

ADD report.ipynb /
ADD nCoV-2019.reference.fasta /
ADD nCoV-2019.reference.gbk /
ADD annot_table.orfs.txt /
