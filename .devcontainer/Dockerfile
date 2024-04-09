FROM mcr.microsoft.com/vscode/devcontainers/base:debian

RUN sudo apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && sudo apt-get install -y \
        cmake \
        git \
        libeigen3-dev \
        libopenmpi-dev \
        r-base-dev 

RUN git clone https://github.com/doctest/doctest.git /doctest \
    && cd /doctest \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make install \
    && cd / \
    && rm -rf /doctest

RUN  /usr/bin/R -q -e "install.packages(c('Rcpp', 'RInside'))"