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
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN conda install -y -c bioconda samtools vcftools bedtools plink2
RUN conda install -y -c bioconda snpeff snpsift
RUN conda install -y -c bioconda numpy pandas


# Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

ENV PATH=/usr/local/bin/snpEff/scripts:$PATH

CMD ["/bin/bash", "-i"]