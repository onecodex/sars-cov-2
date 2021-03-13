FROM python:3.8

# dependencies for generating report
ENV LD_LIBRARY_PATH=/usr/local/lib

RUN pip install numpy

RUN pip install pysam==0.16 biopython==1.78 PyVCF

RUN pip install onecodex[all,reports]==v0.9.4

RUN mkdir -p /usr/local/share/fonts \
    && cp /usr/local/lib/python3.8/site-packages/onecodex/assets/fonts/*.otf /usr/local/share/fonts \
    && fc-cache

USER root
RUN apt-get update \
    && apt-get install -y gnupg \
    && curl -sL https://deb.nodesource.com/setup_14.x  | bash - \
    && apt-get install -y nodejs \
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

ADD covid19_call_variants.sh /usr/local/bin/
ADD covid19_call_variants.artic.sh /usr/local/bin/
ADD post_process_variants.sh /usr/local/bin/
ADD generate_tsv.py /usr/local/bin

ADD report.ipynb /
ADD nCoV-2019.reference.fasta /
ADD nCoV-2019.reference.gtf /
ADD annot_table.orfs.txt /
ADD ivar_variants_to_vcf.py /
