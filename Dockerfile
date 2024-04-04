FROM debian:bookworm-slim AS runenv

LABEL Maintainer "Rob Egan"

# install base requirements

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update  && \
    apt-get install --no-install-recommends -y python3 libcurl4 binutils build-essential perl cmake git zlib1g zlib1g-dev && \
    apt-get install -y wget curl && \
    apt-get autoremove -y && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* 

RUN git clone https://bitbucket.org/berkeleylab/mhm2.git /var/tmp/mhm2-build
RUN cd /var/tmp/mhm2-build && \
    git submodule init && \
    git submodule update

ENV UPCXXVER=2023.9.0
ENV GASNET_PHYSMEM_MAX=1/2
ENV GASNET_PHYSMEM_PROBE=NO

RUN cd /var/tmp/mhm2-build && \
    TMPDIR=/var/tmp/mhm2-build ./upcxx-utils/contrib/install_upcxx.sh /usr/local --with-default-network=smp --disable-ibv && \
    cd /var/tmp/mhm2-build && \
    MHM2_INSTALL_PATH=/usr/local ./build.sh Release && \
    find / -name '*.pyc' -exec rm {} \;

ENV PATH=/usr/local/bin:${PATH}

CMD /bin/bash
