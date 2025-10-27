#!/bin/bash

docker stop mini-viewer-image_v2
docker rm mini-viewer-image_v2

# Change to the root directory (where R/ lives)
cd "$(dirname "$0")/../.." || exit 1

# Build using the root directory as context
docker build --build-arg LOCAL_DATA=/data/dev/miniViewer_3.0 -t mini-viewer-image_v2 -f kalynn_R/latest_app_kalynn/Dockerfile .

# Run the container (as root by default to access mounted directories)
docker run -m 30g -d -p 51175:3838 \
  -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0 \
  -v /data/dev/DO_mapping_files:/data/dev/DO_mapping_files:ro \
  -v /home/khwillis@ad.wisc.edu/qtlApp/data/correlations:/data/correlations:ro \
  --name mini-viewer-image_v2 mini-viewer-image_v2