FROM r-base:4.5.1

RUN apt update \
    && apt install -y r-cran-tidyverse \
                      r-cran-ggbeeswarm \
                      r-cran-genetics \
                      r-cran-ggrepel \
                      r-cran-markdown \
                      r-cran-dt \
                      r-cran-rspectra \
                      r-cran-cowplot \
                      r-cran-plotly \
                      r-cran-data.table \
                      r-cran-knitr \
                      r-cran-biocmanager \
                      build-essential \
                      libcurl4-gnutls-dev \
                      libxml2-dev libssl-dev \
                      libgit2-dev \
                      libc6-dev \
                      pandoc \
                      bedtools \
                      tabix \
    && apt clean

RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('ggnewscale')"
RUN Rscript -e "install.packages('coop')"
RUN Rscript -e "install.packages('fuzzyjoin')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('IRanges')"
RUN Rscript -e "install.packages('sommer')"
RUN Rscript -e "install.packages('valr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('cairo', dependencies=TRUE, repos='http://cran.us.r-project.org')"

    # RUN grep -p "^deb" /etc/apt/sources.list | awk '{printf "%s contrib non-free\n", $0}' > sources.list \
    #     && grep -v -p "^deb" /etc/apt/sources.list >> sources.list \
    #     && mv sources.list /etc/apt/sources.list \
    #     && apt update \
    #     && apt install ttf-mscorefounts-installer \


RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250615.zip -O plink.zip \
    && unzip plink.zip \
    && mv plink /usr/bin/ \
    && rm plink.zip

RUN apt install -y procps make wget gcc zlib1g-dev libgsl-dev libperl-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libxrender1 libcairo2-dev && \
    apt clean && \
	rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 && \
    tar xjf bcftools-1.20.tar.bz2 && \
    cd bcftools-1.20 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools*

RUN Rscript -e "install.packages('R.utils')"
RUN Rscript -e "BiocManager::install('MultiMed')"
RUN apt update && \
    apt purge -y make gcc

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-4-x86_64.zip \
    && unzip gcta-1.94.1-linux-kernel-4-x86_64.zip \
    && mv gcta-1.94.1-linux-kernel-4-x86_64/gcta64 /usr/bin/ \
    && rm -rf gcta*
RUN Rscript -e "install.packages('mediation')"
