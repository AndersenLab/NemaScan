FROM mambaorg/micromamba:1.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml .

RUN micromamba install -n base -f conda.yml -y \
	&& micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN Rscript -e "install.packages('valr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('cairo', dependencies=TRUE, repos='http://cran.us.r-project.org')"

RUN micromamba install -n base -c conda-forge r-fuzzyjoin
RUN micromamba install -n base -c bioconda bioconductor-iranges
RUN micromamba install -n base -c conda-forge mscorefonts
RUN micromamba install -n base -c bioconda tabix

USER root

RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y procps make wget gcc zlib1g-dev libgsl-dev libperl-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libxrender1 libcairo2-dev && \
    apt-get clean && \
	rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools*

RUN apt-get purge -y make gcc

USER $MAMBA_USER
ENV PATH "/opt/conda/bin:$PATH"
ENTRYPOINT source /usr/local/bin/_entrypoint.sh