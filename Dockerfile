##############################################################################################################################
#
# This container includes all necessary components for initializing the NemaScan Nextflow pipeline in Google Cloud
# Additional configuration options can be passed in via environment variables
#
##############################################################################################################################

# Base image includes Google Cloud SDK tools
FROM google/cloud-sdk:slim

# Install OpenJDK JRE for Nextflow
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends openjdk-11-jre wget procps

LABEL Name="NemaScan-NXF" Author="Sam Wachspress"

# Specify Nextflow version and mode 
# (21.05.0-edge is the first version to support configuring which service account acts as pipeline-runner)
ENV NXF_VER=21.05.0-edge \
  NXF_MODE=google \
  NXF_EDGE=1

WORKDIR /nemascan

# Run the Nextflow install script (version and mode must be piped in to bash during install 
# or nextflow will initially download the latest version and only download and switch to NXF_VER when the container runs)
RUN NXF_VER=21.05.0-edge NXF_MODE=google NXF_EDGE=1 \
    wget -qO- https://get.nextflow.io | bash

COPY nemascan-nxf.sh /nemascan/nemascan-nxf.sh
COPY nextflow.config /nemascan/nextflow.config
COPY main.nf /nemascan/main.nf
COPY conf/* /nemascan/conf/
COPY bin/* /nemascan/bin/


# add nextflow and nemarun directory to te system path and make them executable
ENV PATH="/nemascan:${PATH}"
RUN chmod +x /nemascan/nemascan-nxf.sh /nemascan/nextflow

WORKDIR /nemascan
