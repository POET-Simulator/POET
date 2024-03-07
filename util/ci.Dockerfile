FROM debian

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y \
        cmake-curses-gui \
        clangd \
        git \
        libeigen3-dev \
        libopenmpi-dev \
        r-cran-rcpp \
        r-cran-rinside \
        libssl-dev