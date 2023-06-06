############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM r-base

# File Author / Maintainer
# MAINTAINER Mamana Mbiyavanga "mamana.mbiyavanga@uct.ac.za"


################## BEGIN INSTALLATION ######################
# Install Basic tools

## install R packages
RUN R --slave -e 'install.packages(c("venneuler", "tidyr",, "optparse", "data.table", "ggplot2", "dplyr", "reshape2", "grid", "cairo", "sm"), repos="https://cloud.r-project.org/", dependencies=TRUE)'

# Tidy up
RUN rm -rf /tmp/*

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]