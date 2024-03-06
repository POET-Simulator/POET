FROM mcr.microsoft.com/vscode/devcontainers/base:debian

RUN sudo apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && sudo apt-get install -y \
        cmake-curses-gui \
        clangd \
        git \
        libeigen3-dev \
        libopenmpi-dev \
        r-cran-rcpp \
        r-cran-rinside 

