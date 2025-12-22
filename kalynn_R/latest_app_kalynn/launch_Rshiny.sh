#!/bin/bash
set -euo pipefail

# Usage:
#   ./kalynn_R/latest_app_kalynn/launch_Rshiny.sh prod
#   ./kalynn_R/latest_app_kalynn/launch_Rshiny.sh dev
#
# Defaults:
# - prod: http://attie.diabetes.wisc.edu:51175/
# - dev:  http://attie.diabetes.wisc.edu:51176/

target="${1:-prod}"

current_branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"

if [[ "${target}" == "prod" ]]; then
  if [[ "${current_branch}" != "main" ]]; then
    echo "Refusing to deploy 'prod' from branch '${current_branch}'."
    echo "Use a checkout/worktree on 'main' to deploy prod (port 51175)."
    exit 3
  fi
  container_name="mini-viewer-prod"
  image_name="mini-viewer:prod"
  host_port="51175"
  data_root="/data/prod/miniViewer_3.0"
elif [[ "${target}" == "dev" ]]; then
  if [[ "${current_branch}" != "develop" && "${current_branch}" != "dev" ]]; then
    echo "Refusing to deploy 'dev' from branch '${current_branch}'."
    echo "Use a checkout/worktree on 'develop' (or 'dev') to deploy dev (port 51176)."
    exit 3
  fi
  container_name="mini-viewer-dev"
  image_name="mini-viewer:dev"
  host_port="51176"
  data_root="/data/dev/miniViewer_3.0"
else
  echo "Unknown target '${target}'. Use 'prod' or 'dev'."
  exit 2
fi

# Change to the repo root directory (where app.R and R/ live)
cd "$(dirname "$0")/../.." || exit 1
repo_root="$(pwd)"

# Stop only the target container (so prod/dev can run simultaneously)
docker rm -f "${container_name}" >/dev/null 2>&1 || true

# Build using the repo root as context
docker build -t "${image_name}" -f kalynn_R/latest_app_kalynn/Dockerfile .

# Run the container
docker run -m 30g -d -p "${host_port}:3838" \
  -e MINIVIEWER_DATA_ROOT="${data_root}" \
  -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0:ro \
  -v /data/prod/miniViewer_3.0:/data/prod/miniViewer_3.0:ro \
  -v /data/dev/DO_mapping_files:/data/dev/DO_mapping_files:ro \
  -v "${repo_root}/data/correlations:/data/correlations:ro" \
  --name "${container_name}" "${image_name}"
