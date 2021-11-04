############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:18.04

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
    openjdk-8-jdk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R from source
ENV R_VERSION=3.6.3
RUN curl -O https://cran.rstudio.com/src/base/R-3/R-${R_VERSION}.tar.gz && \
    tar -xzvf R-${R_VERSION}.tar.gz && \
    cd R-${R_VERSION} && \
    ./configure \
    --prefix=/opt/R/${R_VERSION} \
    --enable-memory-profiling \
    --enable-R-shlib \
    --with-blas \
    --with-lapack && \
    make && \
    sudo make install

## install R packages
RUN R --slave -e 'install.packages("dplyr", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("data.table", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("sm", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("reshape2", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("cairo", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("optparse", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("tidyr", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("rJava", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("UpSetR", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("tidyverse", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("venneuler", repos="https://cloud.r-project.org/")'
RUN R --slave -e 'install.packages("grid", repos="https://cloud.r-project.org/")'


# Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]