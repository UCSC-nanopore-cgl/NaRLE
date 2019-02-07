FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV QUAST_DIR=quast-5.0.2
ENV QUAST_TAR=$QUAST_DIR.tar.gz
ENV QUAST_URL=https://downloads.sourceforge.net/project/quast/$QUAST_TAR

# update and install dependencies
RUN apt-get update && \
    apt-get install -y git time wget pkg-config libfreetype6-dev libpng-dev python-matplotlib  \
        gcc g++ zlib1g-dev make && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# get app
WORKDIR /opt
RUN wget --quiet $QUAST_URL && \
    tar -xvf $QUAST_TAR && \
    rm $QUAST_TAR && \
    cd /opt/$QUAST_DIR/ && \
    /opt/conda/bin/python setup.py install_full

RUN /opt/conda/bin/python3 /opt/$QUAST_DIR/quast.py --help

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
