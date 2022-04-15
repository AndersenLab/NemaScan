FROM continuumio/miniconda
MAINTAINER Katie Evans <kathrynevans2015@u.northwestern.edu>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
RUN conda install -c conda-forge r-lme4
RUN R -e "install.packages('mediation',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('MultiMed')"
RUN conda install -c conda-forge mscorefonts
