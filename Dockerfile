# Use base image
FROM ubuntu:14.04

# Define directories
ENV OUTPUT_DIR /data
ENV WORK_DIR /code
ENV OPT_DIR /opt

# Mount the output volume as persistant
VOLUME ${OUTPUT_DIR}

RUN \
    # Install Packages Dependencies
    apt-get update -yqq && \
    apt-get install -yqq \
        zlib1g-dev \
        pkg-config \
        curl \
        git \
        locales \
        python-pip \
        wget && \
    apt-get clean && \
    \
    ## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Install click_mergevcfs
COPY . ${WORK_DIR}
WORKDIR ${WORK_DIR}
RUN \
    pip install --upgrade setuptools && \
    pip install --editable .

# Install HTSlib
RUN \
    cd /tmp && \
    wget "https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2" && \
    tar -vxjf htslib-1.8.tar.bz2 && \
    cd htslib-1.8 && \
    ./configure --disable-bz2 --disable-lzma && \
    make && \
    make install

# Install vcftools
RUN \
    cd /tmp && \
    wget "https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz" && \
    tar xvf vcftools-0.1.15.tar.gz && \
    cd vcftools-0.1.15 && \
    ./configure && \
    make && \
    make install

# Run command
ENTRYPOINT ["click_mergevcfs"]
