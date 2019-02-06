FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV MM2_GIT_COMMIT=a5eafb75f91b627cc8a5112c5dda1da7106250f0
ENV CAS_GIT_COMMIT=36ff7a7d04f0d9fe6786a7957754c01da2dc0453
ENV ST_ZIP=samtools-1.9.tar.bz2
ENV ST_URL=https://github.com/samtools/samtools/releases/download/1.9/
ENV ST_FOLDER=samtools-1.9

# update and install dependencies
RUN apt-get update && \
    apt-get install -y git time wget bzip2 ca-certificates \
        libglib2.0-0 libxext6 libsm6 libxrender1 python3-pyqt5 \
        libncurses5-dev libncursesw5-dev build-essential libbz2-dev liblzma-dev \
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

RUN /opt/conda/bin/pip install pysam

# get samtools
RUN wget $ST_URL/$ST_ZIP -O /tmp/$ST_ZIP && \
    tar xvf /tmp/$ST_ZIP -C /tmp && \
    rm /tmp/$ST_ZIP && \
    cd /tmp/$ST_FOLDER && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/$ST_FOLDER

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
RUN /opt/conda/bin/python3 /opt/nanopore_assembly_and_polishing_assessment/align_and_summarize_contigs.py --help

# setup entrypoint
RUN mkdir /data
WORKDIR /data
COPY wrapper.sh /opt/
ENTRYPOINT ["bash", "/opt/wrapper.sh"]
