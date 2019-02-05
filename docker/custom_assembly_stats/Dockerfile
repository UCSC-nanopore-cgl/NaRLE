FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV MM2_GIT_COMMIT=a5eafb75f91b627cc8a5112c5dda1da7106250f0
ENV CAS_GIT_COMMIT=01d769ef68af1a3cace1e0e33a0e64916ed9e088

# update and install dependencies
RUN apt-get update && \
    apt-get install -y git time libpng12-dev libfreetype6-dev pkg-config libbz2-dev liblzma-dev python3.5 python3-pip \
        gcc g++ zlib1g-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN pip3 install numpy pysam matplotlib==2.1.0
ENV LANG C.UTF-8
RUN pip3 install subprocess.run


# get minimap
WORKDIR /opt
RUN git clone https://github.com/lh3/minimap2.git && \
    cd /opt/minimap2 && \
    git fetch && \
    git checkout $MM2_GIT_COMMIT && \
    rm -rf /opt/minimap2/.git
WORKDIR /opt/minimap2
RUN make
RUN cp /opt/minimap2/minimap2 /usr/local/bin

# get app
WORKDIR /opt
RUN git clone https://github.com/rlorigro/nanopore_assembly_and_polishing_assessment.git && \
    cd /opt/nanopore_assembly_and_polishing_assessment && \
    git fetch && \
    git checkout $CAS_GIT_COMMIT && \
    rm -rf /opt/nanopore_assembly_and_polishing_assessment/.git
RUN python3 /opt/nanopore_assembly_and_polishing_assessment/align_and_summarize_contigs.py --help

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
