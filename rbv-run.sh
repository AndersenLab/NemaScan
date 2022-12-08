#!/bin/bash
docker run -i -t \
  -v "$(pwd)"/data:/data \
  -v "/Users/rbv218/.gcp/andersen-lab-616df5b9ec72.json:/secret" \
  --env-file rbv.env \
  gcr.io/caendr/nemascan-nxf:v0.03-debug \
  bash


#  nemascan-nxf.sh
  
