# hash:sha256:333e8199ccf9d100d33ce8b0ad6640e5e8878dd6dde7deb6a4158698239bc303
FROM registry.codeocean.com/codeocean/r-base:3.6.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/shiny-server-trusty/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        bedtools=2.26.0+dfsg-5 \
        gsl-bin=2.4+dfsg-6 \
        libgsl-dev=2.4+dfsg-6 \
        libmagick++-dev=8:6.9.7.4+dfsg-16ubuntu6.9 \
        pandoc=1.19.2.4~dfsg-1build4 \
        perl=5.26.1-6ubuntu0.5 \
        samtools=1.7-1 \
        shiny-server=1.5.12.933 \
        wget=1.19.4-1ubuntu2.2 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("Matrix", "1.3-2")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-2")' \
    && Rscript -e 'remotes::install_version("UpSetR", "1.4.0")' \
    && Rscript -e 'remotes::install_version("broom", "0.7.5")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("data.table", "1.14.0")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.0.4")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.3.3")' \
    && Rscript -e 'remotes::install_version("ggpubr", "0.4.0")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("magick", "2.6.0")' \
    && Rscript -e 'remotes::install_version("pheatmap", "1.0.12")' \
    && Rscript -e 'remotes::install_version("reshape", "0.8.8")' \
    && Rscript -e 'remotes::install_version("shiny", "1.6.0")' \
    && Rscript -e 'remotes::install_version("simpleCache", "0.4.1")' \
    && Rscript -e 'remotes::install_version("stringr", "1.4.0")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-7")' \
    && Rscript -e 'remotes::install_version("tibble", "3.0.6")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BSgenome.Hsapiens.UCSC.hg38", \
        "DESeq2", \
        "GenomicAlignments", \
        "GenomicRanges", \
        "LOLA", \
        "SummarizedExperiment", \
        "chromVAR", \
        "edgeR", \
        "genomation", \
        "qvalue" \
    ))' # Original versions: 1.4.1 1.26.0 1.22.1 1.38.0 1.16.0 1.16.1 1.8.0 3.28.1 1.18.0 2.18.0

COPY postInstall /
RUN /postInstall
