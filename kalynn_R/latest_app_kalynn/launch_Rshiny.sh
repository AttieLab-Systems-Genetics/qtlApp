#!/bin/bash

docker stop mini-viewer-image_v2

docker rm mini-viewer-image_v2

# Change to the directory containing the Dockerfile
cd "$(dirname "$0")" || exit 1

# Remove --no-cache to allow using cached layers for faster builds
# Add --no-cache to force a rebuild without using cache
docker build --no-cache --build-arg LOCAL_DATA=/data/dev/miniViewer_3.0 -t mini-viewer-image_v2 .

# Change back to the original directory (optional, but good practice)
cd -

#RUN IT
docker run -m 30g -d -p 51175:3838 -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0 --name mini-viewer-image_v2 mini-viewer-image_v2
