FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN Rscript -e "install.packages('valr', dependencies=TRUE, repos='http://cran.us.r-project.org')"

# install other tools not avalible on conda cloud
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

RUN conda install -c conda-forge r-fuzzyjoin
RUN conda install -c bioconda bioconductor-iranges

# RUN Rscript -e "install.packages('roperators',dependencies=TRUE, repos='http://cran.us.r-project.org')"
# RUN Rscript -e "install.packages('tidyverse', dependencies=TRUE, repos='http://cran.us.r-project.org')"
# RUN conda install -c conda-forge r-tidyverse
# RUN conda install r-valr
# RUN Rscript -e "devtools::install_version('tidyverse',version='1.3.1', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

RUN conda install -c conda-forge mscorefonts
