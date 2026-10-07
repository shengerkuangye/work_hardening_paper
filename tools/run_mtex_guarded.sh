#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  tools/run_mtex_guarded.sh [--timeout SECONDS] [--grace SECONDS] [--matlab PATH] "MATLAB expression"

Runs one MTEX/MATLAB batch in an isolated Linux session. The wrapper:
  - prevents overlapping project MTEX runs with flock;
  - records the exact run session and process inventory;
  - writes a completion marker after the requested MATLAB expression returns;
  - terminates only processes belonging to that run if MATLAB does not exit;
  - enforces a hard timeout and preserves stdout/stderr under .codex_tmp/mtex_runs/.

Example:
  tools/run_mtex_guarded.sh --timeout 7200 \
    "run_comprehensive_ebsd_analysis(pwd,fullfile(pwd,'results','mtex_ebsd_comprehensive'))"
EOF
}

PROJECT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
MATLAB_BIN=${MATLAB_BIN:-/usr/local/MATLAB/R2024a/bin/matlab}
RUN_TIMEOUT_SECONDS=${MTEX_TIMEOUT_SECONDS:-7200}
COMPLETION_GRACE_SECONDS=${MTEX_COMPLETION_GRACE_SECONDS:-20}

while (($# > 0)); do
  case "$1" in
    --timeout)
      (($# >= 2)) || { usage >&2; exit 2; }
      RUN_TIMEOUT_SECONDS=$2
      shift 2
      ;;
    --grace)
      (($# >= 2)) || { usage >&2; exit 2; }
      COMPLETION_GRACE_SECONDS=$2
      shift 2
      ;;
    --matlab)
      (($# >= 2)) || { usage >&2; exit 2; }
      MATLAB_BIN=$2
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
    *)
      break
      ;;
  esac
done

(($# > 0)) || { usage >&2; exit 2; }
[[ $RUN_TIMEOUT_SECONDS =~ ^[1-9][0-9]*$ ]] || {
  printf 'Timeout must be a positive integer.\n' >&2
  exit 2
}
[[ $COMPLETION_GRACE_SECONDS =~ ^[1-9][0-9]*$ ]] || {
  printf 'Grace period must be a positive integer.\n' >&2
  exit 2
}
[[ -x $MATLAB_BIN ]] || {
  printf 'MATLAB executable not found or not executable: %s\n' "$MATLAB_BIN" >&2
  exit 127
}
command -v flock >/dev/null || {
  printf 'Required command not found: flock\n' >&2
  exit 127
}
command -v setsid >/dev/null || {
  printf 'Required command not found: setsid\n' >&2
  exit 127
}

MATLAB_EXPRESSION=$*
RUN_ROOT="$PROJECT_ROOT/.codex_tmp/mtex_runs"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_$$"
RUN_DIR="$RUN_ROOT/$RUN_ID"
LOCK_FILE="$PROJECT_ROOT/.codex_tmp/mtex.lock"
mkdir -p "$RUN_DIR"

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  printf 'Another guarded MTEX run holds %s.\n' "$LOCK_FILE" >&2
  exit 75
fi

export MTEX_PROJECT_ROOT="$PROJECT_ROOT"
export MTEX_RUN_DIR="$RUN_DIR"
DONE_MARKER="$RUN_DIR/complete.marker"
FAILURE_MARKER="$RUN_DIR/failure.marker"
STDOUT_LOG="$RUN_DIR/stdout.log"
STDERR_LOG="$RUN_DIR/stderr.log"
PROCESS_BEFORE="$RUN_DIR/processes_before.txt"
PROCESS_AFTER="$RUN_DIR/processes_after.txt"
RUN_METADATA="$RUN_DIR/run_metadata.txt"

ps -eo pid=,ppid=,pgid=,sid=,comm=,args= > "$PROCESS_BEFORE"
{
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'project_root=%s\n' "$PROJECT_ROOT"
  printf 'matlab_bin=%s\n' "$MATLAB_BIN"
  printf 'timeout_seconds=%s\n' "$RUN_TIMEOUT_SECONDS"
  printf 'completion_grace_seconds=%s\n' "$COMPLETION_GRACE_SECONDS"
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'matlab_expression=%s\n' "$MATLAB_EXPRESSION"
} > "$RUN_METADATA"

MATLAB_BATCH="try; projectRoot=string(getenv('MTEX_PROJECT_ROOT')); runDir=string(getenv('MTEX_RUN_DIR')); cd(projectRoot); addpath(fullfile(projectRoot,'tools','mtex')); ${MATLAB_EXPRESSION}; runDir=string(getenv('MTEX_RUN_DIR')); fid=fopen(fullfile(runDir,'complete.marker'),'w'); assert(fid>=0); fprintf(fid,'success\\n'); fclose(fid); exit(0); catch ME; fprintf(2,'%s\\n',getReport(ME,'extended','hyperlinks','off')); runDir=string(getenv('MTEX_RUN_DIR')); fid=fopen(fullfile(runDir,'failure.marker'),'w'); if fid>=0, fprintf(fid,'%s\\n',ME.identifier); fclose(fid); end; exit(1); end"

RUN_PID=''
RUN_SID=''

session_pids() {
  [[ -n $RUN_SID ]] || return 0
  ps -eo pid=,sid= | awk -v target_sid="$RUN_SID" '$2 == target_sid {print $1}'
}

terminate_run_session() {
  local signal=${1:-TERM}
  local -a owned_pids=()
  mapfile -t owned_pids < <(session_pids)
  if ((${#owned_pids[@]} > 0)); then
    kill -s "$signal" -- "${owned_pids[@]}" 2>/dev/null || true
  fi
}

on_interrupt() {
  printf 'Guard interrupted; terminating run-owned session %s.\n' "${RUN_SID:-unknown}" >&2
  terminate_run_session TERM
  sleep 2
  terminate_run_session KILL
  exit 130
}
trap on_interrupt INT TERM HUP

setsid "$MATLAB_BIN" -batch "$MATLAB_BATCH" >"$STDOUT_LOG" 2>"$STDERR_LOG" &
RUN_PID=$!
RUN_SID=$(ps -o sid= -p "$RUN_PID" | tr -d ' ')
RUN_PGID=$(ps -o pgid= -p "$RUN_PID" | tr -d ' ')
{
  printf 'launcher_pid=%s\n' "$RUN_PID"
  printf 'run_pgid=%s\n' "$RUN_PGID"
  printf 'run_sid=%s\n' "$RUN_SID"
} >> "$RUN_METADATA"

START_EPOCH=$(date +%s)
DONE_EPOCH=''
TIMED_OUT=0
FORCED_AFTER_COMPLETION=0

while kill -0 "$RUN_PID" 2>/dev/null; do
  NOW_EPOCH=$(date +%s)
  if [[ -f $DONE_MARKER ]]; then
    if [[ -z $DONE_EPOCH ]]; then
      DONE_EPOCH=$NOW_EPOCH
    elif ((NOW_EPOCH - DONE_EPOCH >= COMPLETION_GRACE_SECONDS)); then
      printf 'Completion marker exists but MATLAB is still running; terminating run-owned session %s.\n' "$RUN_SID" >&2
      FORCED_AFTER_COMPLETION=1
      terminate_run_session TERM
      sleep 2
      terminate_run_session KILL
      break
    fi
  fi
  if ((NOW_EPOCH - START_EPOCH >= RUN_TIMEOUT_SECONDS)); then
    printf 'MTEX run exceeded %s seconds; terminating run-owned session %s.\n' "$RUN_TIMEOUT_SECONDS" "$RUN_SID" >&2
    TIMED_OUT=1
    terminate_run_session TERM
    sleep 2
    terminate_run_session KILL
    break
  fi
  sleep 1
done

set +e
wait "$RUN_PID"
MATLAB_EXIT=$?
set -e

terminate_run_session TERM
sleep 1
terminate_run_session KILL
ps -eo pid=,ppid=,pgid=,sid=,comm=,args= > "$PROCESS_AFTER"

REMAINING_PIDS=$(session_pids | paste -sd, -)
{
  printf 'matlab_exit=%s\n' "$MATLAB_EXIT"
  printf 'timed_out=%s\n' "$TIMED_OUT"
  printf 'forced_after_completion=%s\n' "$FORCED_AFTER_COMPLETION"
  printf 'remaining_run_session_pids=%s\n' "$REMAINING_PIDS"
  printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} >> "$RUN_METADATA"

if [[ -n $REMAINING_PIDS ]]; then
  printf 'Run-owned processes remain after cleanup: %s\n' "$REMAINING_PIDS" >&2
  exit 70
fi
if ((TIMED_OUT)); then
  exit 124
fi
if [[ -f $FAILURE_MARKER ]]; then
  printf 'MATLAB reported failure. See %s and %s.\n' "$STDOUT_LOG" "$STDERR_LOG" >&2
  exit "${MATLAB_EXIT:-1}"
fi
if [[ ! -f $DONE_MARKER ]]; then
  printf 'MATLAB exited without a completion marker. See %s and %s.\n' "$STDOUT_LOG" "$STDERR_LOG" >&2
  exit "${MATLAB_EXIT:-1}"
fi

printf 'MTEX run completed: %s\n' "$RUN_DIR"
if ((FORCED_AFTER_COMPLETION)); then
  printf 'The calculation completed, and the guard closed a run-owned MATLAB session that did not exit within the grace period.\n'
fi
