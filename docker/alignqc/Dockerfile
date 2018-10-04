FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV GIT_COMMIT=2b471c2bd76f10383ea4f1d951486c83e6816e15

# update and install dependencies
RUN apt-get update && apt-get -y install git r-base-core python2.7 python-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# get napper
WORKDIR /opt
RUN git clone https://github.com/jason-weirather/AlignQC.git && \
    cd /opt/AlignQC && \
    git fetch && \
    git checkout $GIT_COMMIT && \
    rm -rf /opt/AlignQC/.git

# build
WORKDIR /opt/AlignQC
RUN pip install .

RUN apt-get update && apt-get -y install time

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
