FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV GIT_COMMIT=2467e31596ca4b2c8a0925d595b3968e541c8150

# update and install dependencies
RUN apt-get update && \
    apt-get -y install git make wget autoconf gcc g++ bzip2 libcurl4-openssl-dev libpthread-stubs0-dev \
        libbz2-dev liblzma-dev libssl-dev time && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install cmake
WORKDIR /tmp
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.sh && \
    sh /tmp/cmake-3.11.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# get marginPolish
WORKDIR /opt
RUN git clone https://github.com/benedictpaten/marginPhase && \
    cd /opt/marginPhase && \
    git fetch && \
    git checkout $GIT_COMMIT && \
    git submodule update --init && \
    rm -rf /opt/marginPhase/.git

WORKDIR /opt/marginPhase/build
RUN cmake .. && make && ./marginPolish

# setup entrypoint
COPY wrapper.sh /opt/
ENTRYPOINT ["sh", "/opt/wrapper.sh"]
