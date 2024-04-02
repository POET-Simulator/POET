FROM mcr.microsoft.com/vscode/devcontainers/base:debian

RUN sudo apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && sudo apt-get install -y \
        cmake-curses-gui \
        git \
        libeigen3-dev \
        libopenmpi-dev \
        r-cran-rcpp \
        r-cran-rinside \
        gdb

RUN git clone https://github.com/doctest/doctest.git /doctest \
    && cd /doctest \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make install \
    && cd / \
    && rm -rf /doctest
