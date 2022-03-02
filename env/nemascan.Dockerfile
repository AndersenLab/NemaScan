FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

# install other tools not avalible on conda cloud
RUN apt-get update && apt-get install -y procps  
RUN R -e "install.packages('roperators',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_version('tidyverse', version = '1.3.0', repos = 'http://cran.us.r-project.org')"
