FROM gcc:11.2.0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
    git \
    sudo \
    ninja-build \
    fzf \
    libblas-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

RUN curl -Ls https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.5.tar.gz | tar xz && \
    cd openmpi-5.0.5 && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    rm -rf /tmp/openmpi-5.0.5

RUN curl -Lo cmake.sh https://github.com/Kitware/CMake/releases/download/v3.30.5/cmake-3.30.5-linux-x86_64.sh \
    && bash ./cmake.sh --skip-license --prefix=/usr/local \
    && rm cmake.sh

RUN curl -Ls https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.0.tar.gz | tar xz && \
    cd lapack-3.12.0 && \
    mkdir build && \
    cd build && \
    cmake .. -G Ninja -DBUILD_SHARED_LIBS=ON && \
    ninja install && \
    rm -rf /tmp/lapack-3.12.0

RUN curl -Ls https://cran.r-project.org/src/base/R-4/R-4.4.2.tar.gz | tar xz && \
    cd R-4.4.2 && \
    ./configure --prefix=/usr/local --enable-R-shlib --with-blas --with-lapack && \
    make -j $(nproc) && \
    make install && \
    rm -rf /tmp/R-4.4.2

RUN  /usr/local/bin/R -q -e "install.packages(c('Rcpp', 'RInside', 'qs'), repos='https://cran.rstudio.com/')"

RUN curl -Ls https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2 | tar xj \
    && mkdir /tmp/eigen-3.4.0/build \
    && cd /tmp/eigen-3.4.0/build \
    && cmake .. -G Ninja \
    && ninja install \
    && rm -rf /tmp/eigen-3.4.0

RUN git clone https://github.com/doctest/doctest.git /doctest \
    && cd /doctest \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make install \
    && rm -rf /doctest

RUN useradd -m -s /bin/bash -G sudo vscode \
    && echo "vscode ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

USER vscode

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.2.1/zsh-in-docker.sh)" -- \ 
    -t agnoster 

WORKDIR /home/vscode
