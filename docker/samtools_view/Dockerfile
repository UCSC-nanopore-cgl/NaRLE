FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV ZIP=samtools-1.9.tar.bz2
ENV URL=https://github.com/samtools/samtools/releases/download/1.9/
ENV FOLDER=samtools-1.9
ENV DST=/tmp

# update and install dependencies
RUN apt-get update && \
    apt-get install -y libncurses5-dev libncursesw5-dev build-essential zlib1g-dev libbz2-dev liblzma-dev wget time && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# get files and install
RUN wget $URL/$ZIP -O $DST/$ZIP && \
    tar xvf $DST/$ZIP -C $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    make && \
    make install && \
    cd / && \
    rm -rf $DST/$FOLDER

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
