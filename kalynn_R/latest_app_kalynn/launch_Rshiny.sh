#!/bin/bash

docker stop mini-viewer-image_v2

docker rm mini-viewer-image_v2

# Change to the GRANDPARENT directory (qtlApp root) containing the Dockerfile's context
# cd "$(dirname "$0")" || exit 1
# cd "$(dirname "$0")/.." || exit 1
cd "$(dirname "$0")/../.." || exit 1

# DEBUG: List directory contents before building to confirm context
pwd
echo "Listing contents of build context directory:"
ls -al
echo "---------------------------------------------"

# Build using the parent directory as context and specify Dockerfile path
docker build --build-arg LOCAL_DATA=/data/dev/miniViewer_3.0 -t mini-viewer-image_v2 -f kalynn_R/latest_app_kalynn/Dockerfile .

# Change back to the original directory (optional, but good practice)
cd -

#RUN IT
docker run -m 30g -d -p 51175:3838 -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0 --name mini-viewer-image_v2 mini-viewer-image_v2
