FROM mambaorg/micromamba:1.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER gcta_conda.yml .

RUN micromamba install -n base -f gcta_conda.yml -y \
	&& micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps && \
    apt-get clean && \
	rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

ENTRYPOINT source /usr/local/bin/_entrypoint.sh
