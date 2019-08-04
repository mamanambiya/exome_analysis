############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Mamana Mbiyavanga "mamana.mbiyavanga@uct.ac.za"


################## BEGIN INSTALLATION ######################
# Install Basic tools

# Install wget
RUN apt-get update && \
    apt-get install -y \
    autoconf \
    build-essential \
    git \
    libncurses5-dev \
    pkg-config \
    unzip \
    wget curl \
    python python-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    gfortran \
    libgsl0-dev \
    libopenblas-base \
    libopenblas-dev \
    liblapacke-dev \
    libfontconfig1 \
    libxi6 \
    libxrender1 \
    libxext6 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/miniconda2 && \
    rm ~/miniconda.sh

# use conda to install some bioinfo tools
RUN export PATH=/opt/miniconda2/bin:$PATH && \
    conda config --add channels bioconda && \
    conda config --add channels r && \
    conda install -y vcftools bcftools htslib samtools bedtools plink2 snpeff snpsift && \
    conda install r-base=3.4.1 r-ggplot2 r-dplyr r-reshape2 r-data.table r-cairo r-sm r-optparse r-ggsci r-tidyr

# Compile EIG
RUN git clone --depth 1 --branch master https://github.com/DReichLab/EIG /tmp/EI && \
    cd /tmp/EIG/sr && \
    make clobber OPENBLAS=/opt/openblas LDLIBS="-llapacke" && \
    make install OPENBLAS=/opt/openblas LDLIBS="-llapacke"  && \
    mv ../bin/* /usr/local/bin

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

# Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]
