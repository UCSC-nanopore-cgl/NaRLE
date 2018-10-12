FROM java:8-jdk
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# install deps
RUN apt-get update && apt-get -y install time tabix && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Copy over jar and place in /opt/cgl-docker-lib
RUN mkdir /opt/rtg_tools
COPY rtg-tools.jar  /opt/rtg_tools/
COPY LICENSE.txt /opt/rtg_tools/
COPY wrapper.sh /opt/rtg_tools/

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/rtg_tools/wrapper.sh"]
