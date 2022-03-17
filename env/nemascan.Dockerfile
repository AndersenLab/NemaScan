FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

# install other tools not avalible on conda cloud
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
# RUN Rscript -e "install.packages('roperators',dependencies=TRUE, repos='http://cran.us.r-project.org')"
# RUN Rscript -e "install.packages('tidyverse', dependencies=TRUE, repos='http://cran.us.r-project.org')"
# RUN conda install -c conda-forge r-tidyverse
# RUN conda install r-valr
RUN Rscript -e "install.packages('valr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('tidyverse', dependencies=TRUE, repos='http://cran.us.r-project.org')"
