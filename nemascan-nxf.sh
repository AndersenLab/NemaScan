#!/bin/bash
#
# This script acts as a wrapper around the execution of Nextflow passing environment variables as arguments
#
###################################################################################################################

DEFAULT_DATA_DIR="gs://nf-pipelines/NemaScan/input_data"
DEFAULT_VCF_VERSION="20220216"
DEFAULT_GOOGLE_PROJECT="andersen-lab"
DEFAULT_QUEUE_REGION="us-east1"
DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL="nscalc-201573431837@andersen-lab.iam.gserviceaccount.com"
DEFAULT_SPECIES="c_elegans"


# Environment variables with default values:

if [[ -z "${VCF_VERSION}" ]]; then
  VCF_VERSION=${DEFAULT_VCF_VERSION}
  echo "VCF_VERSION environment variable is not set - defaulting to ${VCF_VERSION}"
fi

if [[ -z "${DATA_DIR}" ]]; then
  DATA_DIR=${DEFAULT_DATA_DIR}
  echo "DATA_DIR environment variable is not set - defaulting to ${DATA_DIR}"
fi

if [[ -z "${GOOGLE_PROJECT}" ]]; then
  GOOGLE_PROJECT=${DEFAULT_GOOGLE_PROJECT}
  echo "GOOGLE_PROJECT environment variable is not set - defaulting to ${GOOGLE_PROJECT}"
fi

if [[ -z "${QUEUE_REGION}" ]]; then
  QUEUE_REGION=${DEFAULT_QUEUE_REGION}
  echo "QUEUE_REGION environment variable is not set - defaulting to ${QUEUE_REGION}"
fi

if [[ -z "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" ]]; then
  GOOGLE_SERVICE_ACCOUNT_EMAIL=${DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL}
  echo "GOOGLE_SERVICE_ACCOUNT_EMAIL environment variable is not set - defaulting to ${GOOGLE_SERVICE_ACCOUNT_EMAIL}"
fi

if [[ -z "${SPECIES}" ]]; then
  SPECIES=${DEFAULT_SPECIES}
  echo "SPECIES environment variable is not set - defaulting to ${SPECIES}"
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
  --google_project "${GOOGLE_PROJECT}" \
  --google_zone "${QUEUE_REGION}" \
  --google_service_account_email "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" \
  --traitfile "${TRAIT_FILE}" \
  --species "${SPECIES}" \
  --vcf "${VCF_VERSION}" \
  --work_dir "${WORK_DIR}" \
  --out "${OUTPUT_DIR}" \
  --data_dir "${DATA_DIR}" \
  --work_dir "${WORK_DIR}"
  