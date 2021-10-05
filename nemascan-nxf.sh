#!/bin/bash
#
# This script acts as a wrapper around the execution of Nextflow passing environment variables as arguments
#
###################################################################################################################

DEFAULT_DATA_DIR="gs://nf-pipelines/NemaScan/input_data"
DEFAULT_VCF_VERSION="20210121"

# Environment variables with default values:

if [[ -z "${VCF_VERSION}" ]]; then
  VCF_VERSION=${DEFAULT_VCF_VERSION}
  echo "VCF_VERSION environment variable is not set - defaulting to ${VCF_VERSION}"
fi

if [[ -z "${DATA_DIR}" ]]; then
  DATA_DIR=${DEFAULT_DATA_DIR}
  echo "DATA_DIR environment variable is not set - defaulting to ${DATA_DIR}"
fi


# Environment variables that MUST be set

if [[ -z "${TRAIT_FILE}" ]]; then
  echo "TRAIT_FILE environment variable must be set to the Google Storage path of the data"
  exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  echo "OUTPUT_DIR environment variable must be set to the Google Storage path of the output directory"
  exit 1
fi

if [[ -z "${WORK_DIR}" ]]; then
  echo "WORK_DIR environment variable must be set to the Google Storage path of the working directory"
  exit 1
fi


nextflow run main.nf \
  -profile gcp \
  --traitfile "${TRAIT_FILE}" \
  --vcf "${VCF_VERSION}" \
  --work_dir "${WORK_DIR}" \
  --out "${OUTPUT_DIR}" \
  --data_dir "${DATA_DIR}"
