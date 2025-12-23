#!/bin/bash
set -euo pipefail

# Usage:
#   ./kalynn_R/latest_app_kalynn/launch_Rshiny.sh prod [host_port]
#   ./kalynn_R/latest_app_kalynn/launch_Rshiny.sh dev  [host_port]
#
# Defaults:
# - prod: http://attie.diabetes.wisc.edu:51175/
# - dev:  http://attie.diabetes.wisc.edu:51173/

target="${1:-prod}"
override_port="${2:-}"

current_branch="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"

if [[ "${target}" == "prod" ]]; then
  if [[ "${current_branch}" != "main" ]]; then
    echo "Refusing to deploy 'prod' from branch '${current_branch}'."
    echo "Use a checkout/worktree on 'main' to deploy prod (port 51175)."
    exit 3
  fi
  container_name="mini-viewer-prod"
  image_name="mini-viewer:prod"
  host_port="${override_port:-51175}"
  data_root="/data/prod/miniViewer_3.0"
elif [[ "${target}" == "dev" ]]; then
  if [[ "${current_branch}" != "develop" && "${current_branch}" != "dev" ]]; then
    echo "Refusing to deploy 'dev' from branch '${current_branch}'."
    echo "Use a checkout/worktree on 'develop' (or 'dev') to deploy dev (port 51173)."
    exit 3
  fi
  container_name="mini-viewer-dev"
  image_name="mini-viewer:dev"
  host_port="${override_port:-51173}"
  data_root="/data/dev/miniViewer_3.0"
else
  echo "Unknown target '${target}'. Use 'prod' or 'dev'."
  exit 2
fi

# Change to the repo root directory (where app.R and R/ live)
cd "$(dirname "$0")/../.." || exit 1
repo_root="$(pwd)"

# Determine where correlation CSVs live on the host.
# Prefer an explicit override, otherwise try the current worktree, then fall back
# to a sibling main worktree (common on the server), and finally warn.
correlations_src="${MINIVIEWER_CORRELATIONS_SRC:-}"
if [[ -z "${correlations_src}" ]]; then
  if compgen -G "${repo_root}/data/correlations/*_corr.csv" > /dev/null; then
    correlations_src="${repo_root}/data/correlations"
  elif compgen -G "${repo_root}/../qtlApp-main/data/correlations/*_corr.csv" > /dev/null; then
    correlations_src="$(readlink -f "${repo_root}/../qtlApp-main/data/correlations" 2>/dev/null || echo "${repo_root}/../qtlApp-main/data/correlations")"
  elif compgen -G "/home/*/qtlApp-main/data/correlations/*_corr.csv" > /dev/null; then
    # Support deployments where dev/prod are launched from different users' worktrees.
    # Pick the first matching qtlApp-main correlations directory under /home.
    candidate_dir="$(ls -d /home/*/qtlApp-main/data/correlations 2>/dev/null | head -n 1)"
    correlations_src="$(readlink -f "${candidate_dir}" 2>/dev/null || echo "${candidate_dir}")"
  else
    correlations_src="${repo_root}/data/correlations"
    echo "WARNING: No '*_corr.csv' files found under '${repo_root}/data/correlations' (or sibling 'qtlApp-main')." >&2
    echo "WARNING: Correlation tab may appear empty unless you set MINIVIEWER_CORRELATIONS_SRC to a directory containing the CSVs." >&2
  fi
fi

# Stop only the target container (so prod/dev can run simultaneously)
docker rm -f "${container_name}" >/dev/null 2>&1 || true

# Build using the repo root as context
docker build -t "${image_name}" -f kalynn_R/latest_app_kalynn/Dockerfile .

# Run the container (cap RAM at 35GB to simulate heavy caching use; set swap equal to RAM to avoid extra swap headroom)
docker run --memory=35g --memory-swap=35g -d -p "${host_port}:3838" \
  -e MINIVIEWER_DATA_ROOT="${data_root}" \
  -v /data/dev/miniViewer_3.0:/data/dev/miniViewer_3.0:ro \
  -v /data/prod/miniViewer_3.0:/data/prod/miniViewer_3.0:ro \
  -v /data/dev/DO_mapping_files:/data/dev/DO_mapping_files:ro \
  -v "${correlations_src}:/data/correlations:ro" \
  --name "${container_name}" "${image_name}"

echo "Container: ${container_name}"
echo "Port: ${host_port}"
echo "Correlations mount: ${correlations_src} -> /data/correlations (ro)"
echo "Observe memory: docker stats --no-stream ${container_name}"
