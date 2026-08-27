#!/usr/bin/env bash
set -euo pipefail

experiment_root="${G12SS4E_ROOT:-/root/v0lv-g12-focused-reachability-stage4e-20260827-v1}"
parent_root="${G12SS4E_PARENT_STAGE4C_ROOT:-/root/v0lv-g12-stopping-semantics-stage4c-20260826-v1}"
cd "$experiment_root"

test -d "$parent_root"
test ! -e DATA_REUSE_COMPLETE.txt
test ! -d data
mkdir -p data

for data_id in g12ss4c_base_01 g12ss4c_sparse_01 g12ss4c_sparse_02; do
  source_path="$parent_root/data/$data_id/observation_bundle.rds"
  target_dir="data/$data_id"
  test -f "$source_path"
  mkdir -p "$target_dir"
  cp --preserve=mode,timestamps "$source_path" "$target_dir/observation_bundle.rds"
done

printf '%s  %s\n' \
  f1a183dd568ad0838f8c037c5de526304550f9c9feef72d0f35f918670ca8257 \
  data/g12ss4c_base_01/observation_bundle.rds \
  bcfe2f856d9e8c41c6b6810dd52abbc76bb1e3664d9701178138beac588f131b \
  data/g12ss4c_sparse_01/observation_bundle.rds \
  895afc709321aad5f8916316c0569c69d86970ca0013f31009d791c41e858c81 \
  data/g12ss4c_sparse_02/observation_bundle.rds | sha256sum -c --status

if find data -type f | grep -Eq 'sealed_truth|truth_bundle'; then
  echo "Truth material entered the Stage-4E fit root." >&2
  exit 2
fi

printf 'status=PASS\nprepared_utc=%s\nobservation_only_bundles=3\ntruth_bundles=0\nsource_stage=4C\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > DATA_REUSE_COMPLETE.txt
printf 'STAGE4E_DATA_REUSE_PASS observations=3 truth=0\n'
