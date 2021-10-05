#!/bin/bash
#
# This script acts as a wrapper around the execution of Nextflow passing environment variables as arguments
#
###################################################################################################################


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

if [[ -z "${VCF_VERSION}" ]]; then
  VCF_VERSION="20210121"
  echo "VCF_VERSION environment variable is not set - defaulting to ${VCF_VERSION}"
fi

nextflow run main.nf \
  -profile gcp \
  --traitfile "${TRAIT_FILE}" \
  --vcf "${VCF_VERSION}" \
  --work_dir "${WORK_DIR}" \
  --out "${OUTPUT_DIR}" \
  --data_dir "gs://nf-pipelines/NemaScan/new_input_data/input_data/"
