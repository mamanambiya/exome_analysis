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
    libxext6 \
    openjdk-8-jdk \
    r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# r-base r-ggplot2 r-dplyr r-reshape2
#r-cairo r-optparse r-tidyr && \
#conda install -c conda-forge r-sm r-ggsci

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xvf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.9*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.9*

# Install VCFTools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xvf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
  tar -xvf bcftools-1.9.tar.bz2 && \
  cd bcftools-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf bcftools-1.9*

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mv bin/bedtools /usr/local/bin/ && \
    cd .. && rm -r bedtools2

# Install PLINK2
# there is an undocumented stable url (without the date)
RUN wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O plink.zip && \
  unzip plink.zip -d /usr/local/bin/ && \
  rm -f plink.zip


# Compile EIG
RUN git clone --depth 1 --branch master https://github.com/DReichLab/EIG /tmp/EIG && \
    cd /tmp/EIG/src && \
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

# Install snpEff
# Setup ENV variables
ENV \
    SNPEFF_VERSION=latest\
    SNPEFF_HOME=/usr/local/bin/snpEff
RUN \
    wget --quiet -O snpEff_${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_${SNPEFF_VERSION}_core.zip && \
    unzip snpEff_${SNPEFF_VERSION}_core.zip -d /usr/local/bin/ && \
    bash -c 'echo -e "#!/bin/bash\njava -jar /usr/local/bin/snpEff/SnpSift.jar \$@" > /usr/local/bin/snpEff/scripts/snpsift' && \
    chmod +x /usr/local/bin/snpEff/scripts/snpsift && \
    rm snpEff_${SNPEFF_VERSION}_core.zip


## Install python3.7
RUN apt-get -y update && apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev
RUN cd /tmp && \
        wget https://www.python.org/ftp/python/3.7.2/Python-3.7.2.tar.xz && \
        tar -xf Python-3.7.2.tar.xz && \
        cd Python-3.7.2 && \
        ./configure --enable-optimizations && \
        make -j 4 && \
        make altinstall

## install R packages
RUN R --slave -e 'install.packages("dplyr", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("data.table", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("sm", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("reshape2", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("cairo", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("optparse", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("tidyr", repos="https://cloud.r-project.org/")'

## Install pandas for python3.7 
RUN apt-get -y update && apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install numpy 
RUN pip3 install pandas 

# Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

ENV PATH=/usr/local/bin/snpEff/scripts:$PATH

CMD ["/bin/bash","-i"]
