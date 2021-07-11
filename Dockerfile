FROM bitnami/minideb:stretch

MAINTAINER Arnav Gupta

RUN apt update \
  && apt install -y \
  openjdk-8-jre-headless \
  libncursesw5-dev \
  build-essential \
  liblzma-dev \
  zlib1g-dev \
  libbz2-dev \
  bedtools \
  unzip \
  wget \
  make

ENV SAMTOOLS_VERSION="1.13"
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd samtools-${SAMTOOLS_VERSION} \
  && ./configure \
  && make \
  && make install

RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
  && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd htslib-${SAMTOOLS_VERSION} \
  && ./configure \
  && make \
  && make install

ENV GATK_VERSION="4.2.0.0"
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
  && unzip gatk-${GATK_VERSION}.zip \
  && rm gatk-${GATK_VERSION}.zip \
  && mv gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-local.jar .. \
  && rm -r gatk-${GATK_VERSION}

ENV BWA_VERSION="0.7.17"
RUN wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
  && tar xjf bwa-${BWA_VERSION}.tar.bz2 \
  && rm bwa-${BWA_VERSION}.tar.bz2 \
  && cd bwa-${BWA_VERSION} \
  && make
ENV PATH=/bwa-${BWA_VERSION}:$PATH