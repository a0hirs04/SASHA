#!/bin/bash
# Live monitor for a single RC2 job using the SLURM log as the progress source.
# Usage: bash watch_rc2_eta.sh [--once] <SLURM_JOB_ID> [RUN_DIR]

set -euo pipefail

ONCE=0
if [[ "${1:-}" == "--once" ]]; then
    ONCE=1
    shift
fi

JOB_ID="${1:?Usage: bash watch_rc2_eta.sh [--once] <SLURM_JOB_ID> [RUN_DIR]}"
RUN_DIR="${2:-/home/a0hirs04/PROJECT-NORTHSTAR/build/rc2_full_seed42/replicate_01_seed42}"
OUT_DIR="${RUN_DIR}/output"
OUT_LOG="${RUN_DIR}/slurm_${JOB_ID}.out"
ERR_LOG="${RUN_DIR}/slurm_${JOB_ID}.err"
POLL_SECONDS=15

job_state() {
    local state
    state=$(squeue -j "$JOB_ID" --format="%T" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ -n "$state" ]]; then
        printf "%s" "$state"
        return
    fi

    state=$(sacct -j "$JOB_ID" --format=State --noheader --parsable2 2>/dev/null | \
        head -1 | cut -d'|' -f1 | tr -d ' ')
    printf "%s" "${state:-UNKNOWN}"
}

job_elapsed() {
    local elapsed
    elapsed=$(squeue -j "$JOB_ID" --format="%M" --noheader 2>/dev/null | head -1 | tr -d ' ')
    if [[ -n "$elapsed" ]]; then
        printf "%s" "$elapsed"
        return
    fi

    elapsed=$(sacct -j "$JOB_ID" --format=Elapsed --noheader --parsable2 2>/dev/null | \
        head -1 | cut -d'|' -f1 | tr -d ' ')
    printf "%s" "${elapsed:-unknown}"
}

format_seconds() {
    python3 - "$1" <<'PY'
import sys

secs = float(sys.argv[1])
if secs < 0:
    secs = 0.0

days = int(secs // 86400)
secs -= days * 86400
hours = int(secs // 3600)
secs -= hours * 3600
minutes = int(secs // 60)
secs -= minutes * 60
seconds = int(round(secs))

if seconds == 60:
    seconds = 0
    minutes += 1
if minutes == 60:
    minutes = 0
    hours += 1
if hours == 24:
    hours = 0
    days += 1

parts = []
if days:
    parts.append(f"{days}d")
if hours:
    parts.append(f"{hours}h")
if minutes:
    parts.append(f"{minutes}m")
if not parts or seconds:
    parts.append(f"{seconds}s")
print(" ".join(parts))
PY
}

parse_progress() {
    python3 - "$OUT_LOG" <<'PY'
import os
import re
import sys

path = sys.argv[1]
if not os.path.exists(path):
    print("0\t60480\t0\t0\tSTARTUP")
    raise SystemExit(0)

sim_re = re.compile(r"current simulated time:\s*([0-9.eE+-]+)\s*min\s*\(max:\s*([0-9.eE+-]+)\s*min\)")
agents_re = re.compile(r"total agents:\s*(\d+)")
wall_re = re.compile(
    r"total wall time:\s*(\d+)\s*days,\s*(\d+)\s*hours,\s*(\d+)\s*minutes,\s*and\s*([0-9.eE+-]+)\s*seconds"
)

sim_time = 0.0
max_time = 60480.0
agents = 0
wall_seconds = 0.0

with open(path, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        m = sim_re.search(line)
        if m:
            sim_time = float(m.group(1))
            max_time = float(m.group(2))
            continue

        m = agents_re.search(line)
        if m:
            agents = int(m.group(1))
            continue

        m = wall_re.search(line)
        if m:
            days = int(m.group(1))
            hours = int(m.group(2))
            minutes = int(m.group(3))
            seconds = float(m.group(4))
            wall_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds

if sim_time < 20160.0:
    phase = "BARRIER"
elif sim_time < 40320.0:
    phase = "DRUG-ON"
elif sim_time < max_time:
    phase = "REGROWTH"
else:
    phase = "DONE"

print(f"{sim_time}\t{max_time}\t{wall_seconds}\t{agents}\t{phase}")
PY
}

render_bar() {
    python3 - "$1" "$2" <<'PY'
import sys

progress = float(sys.argv[1])
width = int(sys.argv[2])
progress = 0.0 if progress < 0 else progress
progress = 1.0 if progress > 1 else progress
filled = int(round(progress * width))
filled = min(width, max(0, filled))
print("#" * filled + "-" * (width - filled))
PY
}

terminal_state() {
    case "$1" in
        COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

while true; do
    STATE=$(job_state)
    ELAPSED=$(job_elapsed)
    IFS=$'\t' read -r SIM_TIME MAX_TIME WALL_SECONDS AGENTS PHASE <<< "$(parse_progress)"

    PROGRESS=$(python3 - "$SIM_TIME" "$MAX_TIME" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
print(0.0 if max_time <= 0 else sim_time / max_time)
PY
)

    ETA_SECONDS=$(python3 - "$SIM_TIME" "$MAX_TIME" "$WALL_SECONDS" <<'PY'
import sys
sim_time = float(sys.argv[1])
max_time = float(sys.argv[2])
wall_seconds = float(sys.argv[3])
if sim_time <= 0 or wall_seconds <= 0 or max_time <= 0:
    print(-1)
else:
    total_est = wall_seconds * (max_time / sim_time)
    print(max(0.0, total_est - wall_seconds))
PY
)

    BAR=$(render_bar "$PROGRESS" 48)
    SNAP_COUNT=$(find "$OUT_DIR" -maxdepth 1 -name 'output*.xml' 2>/dev/null | wc -l | tr -d ' ')
    SIM_DAY=$(python3 - "$SIM_TIME" <<'PY'
import sys
print(f"{float(sys.argv[1]) / 1440.0:.2f}")
PY
)

    if python3 - "$ETA_SECONDS" <<'PY'
import sys
eta = float(sys.argv[1])
raise SystemExit(0 if eta >= 0 else 1)
PY
    then
        ETA_FMT=$(format_seconds "$ETA_SECONDS")
        FINISH_TS=$(python3 - "$ETA_SECONDS" <<'PY'
import sys
import time
print(int(time.time() + float(sys.argv[1])))
PY
)
        FINISH_AT=$(date -u -d "@${FINISH_TS}" '+%Y-%m-%d %H:%M:%S UTC')
    else
        ETA_FMT="estimating..."
        FINISH_AT="unknown"
    fi

    if [[ "$ONCE" -eq 0 ]]; then
        if [[ -n "${TERM:-}" ]]; then
            clear 2>/dev/null || printf '\033[2J\033[H'
        else
            printf '\033[2J\033[H'
        fi
    fi
    printf "RC2 Monitor\n\n"
    printf "Job ID:        %s\n" "$JOB_ID"
    printf "State:         %s\n" "$STATE"
    printf "SLURM elapsed: %s\n" "$ELAPSED"
    printf "Run dir:       %s\n" "$RUN_DIR"
    printf "Logs:          %s\n" "$OUT_LOG"
    printf "\n"
    printf "Phase:         %s\n" "$PHASE"
    printf "Sim time:      %.0f / %.0f min (day %s / 42.00)\n" "$SIM_TIME" "$MAX_TIME" "$SIM_DAY"
    printf "Agents:        %s\n" "$AGENTS"
    printf "Snapshots:     %s\n" "$SNAP_COUNT"
    printf "Progress:      [%s] %5.1f%%\n" "$BAR" "$(python3 - "$PROGRESS" <<'PY'
import sys
print(float(sys.argv[1]) * 100.0)
PY
)"
    printf "ETA:           %s\n" "$ETA_FMT"
    printf "Finish (UTC):  %s\n" "$FINISH_AT"

    if [[ -f "$ERR_LOG" ]]; then
        printf "\nLast stderr line:\n"
        tail -n 1 "$ERR_LOG" 2>/dev/null || true
    fi

    if terminal_state "$STATE"; then
        printf "\nJob reached terminal state: %s\n" "$STATE"
        exit 0
    fi

    if [[ "$ONCE" -eq 1 ]]; then
        exit 0
    fi

    sleep "$POLL_SECONDS"
done
