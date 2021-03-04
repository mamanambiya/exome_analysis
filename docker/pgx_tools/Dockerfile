############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
# MAINTAINER Mamana Mbiyavanga "mamana.mbiyavanga@uct.ac.za"


################## BEGIN INSTALLATION ######################
# Install Basic tools

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y build-essential autoconf gcc libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && \
    tar -xvf htslib-1.10.2.tar.bz2 && \
    cd htslib-1.10.2 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.10.2*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
    tar -xvf bcftools-1.10.2.tar.bz2 && \
    cd bcftools-1.10.2 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf bcftools-1.10.2*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    echo "conda activate base" >> ~/.bashrc

RUN conda update conda && \
    conda clean --all --yes && \
    conda install -y -c bioconda samtools vcftools bedtools plink plink2 snpeff=4.3 snpsift numpy pandas python=3.7 eigensoft

# Install R and R packages
RUN conda clean --all --yes && \
    conda install -y -c bioconda r-ggplot2 r-UpSetR r-dplyr r-tidyr r-data.table r-circlize && \
    conda install -c r r-hmisc
# conflict: r-optparse, r-sm, r-tidyverse

# Install python3.7 packages
RUN conda clean --all --yes && \
    conda install -y -c bioconda numpy pandas 

# Install minimac3
RUN wget ftp://share.sph.umich.edu/minimac3/Minimac3.v2.0.1.tar.gz  && \
    tar -xzvf Minimac3.v2.0.1.tar.gz  && \
    cd Minimac3/  && \
    make  && \
    mv ./bin/Minimac3 /usr/local/bin/minimac3

# Install minimac4
RUN wget http://debian.mirror.ac.za/debian/pool/main/libs/libstatgen/libstatgen0_1.0.14-5_amd64.deb && \
    dpkg -i libstatgen0_1.0.14-5_amd64.deb && \
    wget http://debian.mirror.ac.za/debian/pool/main/m/minimac4/minimac4_1.0.0-2_amd64.deb && \
    dpkg -i minimac4_1.0.0-2_amd64.deb

# Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    gunzip Eagle_v2.4.1.tar.gz && \
    tar xvf Eagle_v2.4.1.tar && \
    mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
    rm -rf Eagle_v2.4.1

## Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

ENV PATH=/usr/local/bin/snpEff/scripts:$PATH

CMD ["/bin/bash", "-i"]
