FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV GIT_COMMIT=a5eafb75f91b627cc8a5112c5dda1da7106250f0

# update and install dependencies
RUN apt-get update && apt-get -y install git make gcc g++ time zlib1g-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# get napper
WORKDIR /opt
RUN git clone https://github.com/lh3/minimap2.git && \
    cd /opt/minimap2 && \
    git fetch && \
    git checkout $GIT_COMMIT && \
    rm -rf /opt/minimap2/.git

# build
WORKDIR /opt/minimap2
RUN make

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
