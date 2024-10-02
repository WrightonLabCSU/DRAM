# Use miniconda3 as the base image
FROM continuumio/miniconda3

# Set environment variables
ENV KMP_DUPLICATE_LIB_OK=TRUE
ENV TZ='America/Denver'
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_CTYPE=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# Label the image
LABEL maintainer="<mescyphers@gmail.com>"

# Copy necessary files into the image
#COPY ../assets/bbmap /opt
COPY environment.yml /opt/environment.yml

# Make the /opt directory and its contents writable
RUN chmod 777 /opt && chmod 777 /opt/*

# Install required dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    git \
    pigz \
    openjdk-11-jre \
    r-base-core \
    build-essential

# Install required R packages
RUN R --slave -e 'install.packages("tidyverse", repos="https://cloud.r-project.org/")'

# Create conda environment
RUN conda update conda \
    && conda env create -f /opt/environment.yml

# Add conda environment to PATH
ENV PATH /opt/conda/envs/dram2-env/bin:$PATH
ENV CONDA_DEFAULT_ENV dram2-env

# Activate conda environment
RUN echo "conda activate dram2-env" >> ~/.bashrc

# Clean up
RUN apt-get clean
