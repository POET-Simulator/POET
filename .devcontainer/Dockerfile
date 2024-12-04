FROM gcc:11.2.0 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
    sudo \
    git \
    ninja-build \
    libmpfr-dev \
    python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

ARG OPENMPI_VERSION=4.1.1
ADD https://download.open-mpi.org/release/open-mpi/v${OPENMPI_VERSION%.*}/openmpi-${OPENMPI_VERSION}.tar.gz /tmp/openmpi.tar.gz

RUN mkdir openmpi && \ 
    tar xf openmpi.tar.gz -C openmpi --strip-components 1 && \
    cd openmpi && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    rm -rf /tmp/openmpi tmp/openmpi.tar.gz

ARG CMAKE_VERSION=3.30.5
ADD https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh /tmp/cmake.sh

RUN bash ./cmake.sh --skip-license --prefix=/usr/local \
    && rm cmake.sh

ARG LAPACK_VERSION=3.12.0
ADD https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${LAPACK_VERSION}.tar.gz /tmp/lapack.tar.gz

RUN mkdir lapack && \
    tar xf lapack.tar.gz -C lapack --strip-components 1 && \
    cd lapack && \
    mkdir build && \
    cd build && \
    cmake .. -G Ninja -DBUILD_SHARED_LIBS=ON && \
    ninja install && \
    rm -rf /tmp/lapack tmp/lapack.tar.gz

ARG R_VERSION=4.4.2
ADD https://cran.r-project.org/src/base/R-${R_VERSION%%.*}/R-${R_VERSION}.tar.gz /tmp/R.tar.gz

RUN mkdir R && \
    tar xf R.tar.gz -C R --strip-components 1 && \
    cd R && \
    ./configure --prefix=/usr/local --enable-R-shlib --with-blas --with-lapack && \
    make -j $(nproc) && \
    make install && \
    rm -rf /tmp/R tmp/R.tar.gz

RUN  /usr/local/bin/R -q -e "install.packages(c('Rcpp', 'RInside', 'qs'), repos='https://cran.rstudio.com/')"

ARG EIGEN_VERSION=3.4.0
ADD https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.bz2 /tmp/eigen.tar.bz2

RUN mkdir eigen && \
    tar xf eigen.tar.bz2 -C eigen --strip-components 1 && \
    cd eigen && \
    mkdir build && \
    cd build && \
    cmake .. -G Ninja && \
    ninja install && \
    rm -rf /tmp/eigen tmp/eigen.tar.bz2

ARG GDB_VERSION=15.2
ADD https://ftp.gnu.org/gnu/gdb/gdb-${GDB_VERSION}.tar.xz /tmp/gdb.tar.xz

RUN mkdir gdb && \
    tar xf gdb.tar.xz -C gdb --strip-components 1 && \
    cd gdb && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    rm -rf /tmp/gdb tmp/gdb.tar.xz

RUN useradd -m -s /bin/bash -G sudo vscode \
    && echo "vscode ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

USER vscode

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN sudo apt-get update && \
    sudo apt-get install -y zsh && \
    sudo apt-get clean && \
    sudo rm -rf /var/lib/apt/lists/*

RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.2.1/zsh-in-docker.sh)" -- \ 
    -t agnoster \
    -p zsh-syntax-highlighting

RUN zsh -c "git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting"

RUN zsh -c "git clone --depth 1 https://github.com/junegunn/fzf.git ~/.fzf && ~/.fzf/install"

RUN mkdir -p /home/vscode/.config/gdb \
    && echo "set auto-load safe-path /" > /home/vscode/.config/gdb/gdbinit

ENV CMAKE_GENERATOR=Ninja
ENV CMAKE_EXPORT_COMPILE_COMMANDS=ON

WORKDIR /home/vscode
