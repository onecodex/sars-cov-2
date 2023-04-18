FROM python:3.8

# dependencies for generating report
ENV LD_LIBRARY_PATH=/usr/local/lib

# install dnaplotlib for creating the genome diagram
# hard-pin some dependencies for onecodex 0.9.6
RUN pip install numpy pysam==0.16 biopython==1.78 PyVCF dnaplotlib onecodex[all,reports]==v0.9.6 \
    nbconvert==5.6.1 click==8.0.4 Jinja2==3.0.3

USER root
RUN apt-get update \
    && apt-get autoclean \
    && apt-get install -y gnupg curl \
    && curl -sL https://deb.nodesource.com/setup_14.x  | bash - \
    && apt-get install -y nodejs \
    unzip \
    default-jre \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN npm install -g --unsafe-perm vega vega-lite vega-cli canvas

WORKDIR /opt

# install Conda
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  > Miniconda3-latest-Linux-x86_64.sh \
  && yes \
  | bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3

# put system path first so that conda doesn't override python
ENV PATH=$PATH:/opt/miniconda3/bin/

# install environment's dependencies
COPY environment.yml /opt/
RUN conda env create -f /opt/environment.yml

# install artic into conda environment "artic"
RUN git clone https://github.com/artic-network/fieldbioinformatics.git \
        && cd fieldbioinformatics \
        && git checkout 1.2.1 \
      	&& conda env create -f environment.yml \
        && conda run -n artic python setup.py install \
        && conda clean -a

# install snpeff
RUN curl -k -L https://sourceforge.net/projects/snpeff/files/snpEff_v4_5covid19_core.zip/download --output snpEff_v4_5covid19_core.zip\
        && unzip snpEff_v4_5covid19_core.zip \
        && mv snpEff /usr/local/bin \
        && rm snpEff_v4_5covid19_core.zip

COPY reference /reference

# ARTIC's vcf_filter.py breaks when a variant's call score is "."
# I fixed this in reference/vcf_filter_edited.py
COPY /reference/vcf_filter.edited.py /root/miniconda3/envs/artic/lib/python3.6/site-packages/artic-1.2.1-py3.6.egg/artic/vcf_filter.py

# Add ARTIC 4.1 patch as a primer scheme
RUN mkdir /primer_schemes
COPY /reference/primer_schemes/ /primer_schemes/

# install pangolin into conda environment "pangolin"
RUN git clone https://github.com/cov-lineages/pangolin.git \
        && cd pangolin \
        && git checkout v4.2 \
        && conda env create -f environment.yml \
        && conda run -n pangolin python setup.py install \
        && conda clean -a

# install nextclade & download sars-cov-2 dataset
RUN curl -fsSL 'https://github.com/nextstrain/nextclade/releases/download/2.9.1/nextclade-x86_64-unknown-linux-gnu' -o '/usr/local/bin/nextclade' && chmod +x /usr/local/bin/nextclade
RUN /usr/local/bin/nextclade dataset get --name 'sars-cov-2' --output-dir '/usr/local/bin/data/sars-cov-2'

# Setup onecodex_pdf export option
RUN mkdir -p /usr/local/share/fonts \
    && cp /usr/local/lib/python3.8/site-packages/onecodex/assets/fonts/*.otf /usr/local/share/fonts \
    && fc-cache

ADD jobscript.sh /usr/local/bin/
ADD covid19_call_variants.sh /usr/local/bin/
ADD covid19_call_variants.ont.sh /usr/local/bin/
ADD post_process_variants.sh /usr/local/bin/
ADD generate_tsv.py /usr/local/bin
ADD insert_coverage_stats.py /usr/local/bin
ADD report.ipynb /

ENV MPLCONFIGDIR /tmp
RUN chmod -R a+rwx /reference /primer_schemes /report.ipynb

# so we can include git hash in report for tracking
COPY .git /.git
