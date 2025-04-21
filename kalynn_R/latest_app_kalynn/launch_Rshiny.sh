
docker stop mini-viewer-image_v2

docker rm mini-viewer-image_v2

docker build --build-arg LOCAL_DATA=/data/dev/miniViewer_3.0 -t mini-viewer-image_v2 .

#RUN IT
docker run -m 30g -d -p 51175:3838 -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0 --name mini-viewer-image_v2 mini-viewer-image_v2
