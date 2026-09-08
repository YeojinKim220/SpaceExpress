#!/bin/bash
set -euo pipefail
PYTHON=${1:?Python executable}
PYTHON=$(cd "$(dirname "$PYTHON")" && pwd)/$(basename "$PYTHON")
CONFIG=$(realpath "${2:?Config JSON}")
OUTPUT=$(realpath -m "${3:?New output directory}")
HERE=$(cd "$(dirname "$0")" && pwd)
if [[ -e "$OUTPUT/job_ids.tsv" || -e "$OUTPUT/config_resolved.json" ]]; then
  echo "Output already contains a run; choose a new output directory: $OUTPUT" >&2
  exit 1
fi
mkdir -p "$OUTPUT/logs"
COMMON=(--parsable --account="${SE_ACCOUNT:-gts-ssinha338-paid}" --qos="${SE_QOS:-inferno}"
        --nodes=1 --ntasks=1 --time=04:00:00 --kill-on-invalid-dep=yes
        --output="$OUTPUT/logs/%x-%j.out" --error="$OUTPUT/logs/%x-%j.err")
CPU_PARTITION=${SE_CPU_PARTITION:-inferno}
STAGE="$HERE/run_public_pair_stage.sh"
prep=$(sbatch "${COMMON[@]}" --job-name=se015-prep --partition="$CPU_PARTITION" \
       --cpus-per-task=12 --mem=128G "$STAGE" "$PYTHON" "$CONFIG" "$OUTPUT" prepare)
prep=${prep%%;*}
printf 'stage\tjob_id\nprepare\t%s\n' "$prep" > "$OUTPUT/job_ids.tsv"
gpu=$(sbatch "${COMMON[@]}" --job-name=se015-gpu --partition="${SE_GPU_PARTITION:-gpu-h100}" \
      --cpus-per-task=12 --mem=256G --gres=gpu:1 --dependency="afterok:$prep" \
      "$STAGE" "$PYTHON" "$CONFIG" "$OUTPUT" gpu)
gpu=${gpu%%;*}
printf 'probe_and_embedding\t%s\n' "$gpu" >> "$OUTPUT/job_ids.tsv"
deps=""
for k in 30 50 100; do
  job=$(sbatch "${COMMON[@]}" --job-name="se015-dse$k" --partition="$CPU_PARTITION" \
        --cpus-per-task=12 --mem=128G --dependency="afterok:$gpu" \
        "$STAGE" "$PYTHON" "$CONFIG" "$OUTPUT" dse "$k")
  job=${job%%;*}
  printf 'dse_k%s\t%s\n' "$k" "$job" >> "$OUTPUT/job_ids.tsv"
  deps="$deps:$job"
done
report=$(sbatch "${COMMON[@]}" --job-name=se015-report --partition="$CPU_PARTITION" \
         --cpus-per-task=4 --mem=64G --dependency="afterok$deps" \
         "$STAGE" "$PYTHON" "$CONFIG" "$OUTPUT" report)
report=${report%%;*}
printf 'report\t%s\n' "$report" >> "$OUTPUT/job_ids.tsv"
cat "$OUTPUT/job_ids.tsv"
