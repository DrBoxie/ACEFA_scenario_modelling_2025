#!/usr/bin/env bash

set -u

###############################################################################
# Log directory
###############################################################################

LOG_DIR="/pvol/Log files"

rm -rf "$LOG_DIR"
mkdir -p "$LOG_DIR"

RUN_START=$(date +%s)

###############################################################################
# Helper functions
###############################################################################

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

show_runtime() {
    local start_time=$1
    local end_time
    end_time=$(date +%s)

    local elapsed=$((end_time - start_time))
    local total_elapsed=$((end_time - RUN_START))

    log_msg "Block runtime: ${elapsed}s ($(printf '%02d:%02d:%02d' \
        $((elapsed/3600)) \
        $(((elapsed%3600)/60)) \
        $((elapsed%60))))"

    log_msg "Total runtime: ${total_elapsed}s ($(printf '%02d:%02d:%02d' \
        $((total_elapsed/3600)) \
        $(((total_elapsed%3600)/60)) \
        $((total_elapsed%60))))"
}

handle_failure() {
    local step="$1"

    echo
    echo "================================================================="
    echo "ERROR DETECTED"
    echo "Step: $step"
    echo "Time: $(date)"
    echo "================================================================="
    echo

    echo "Memory usage:"
    free -h || true

    echo
    echo "Disk usage:"
    df -h || true

    echo
    echo "Top memory consumers:"
    ps aux --sort=-%mem | head -20 || true

    echo
    echo "If this was an out-of-memory issue, clean up files manually."
    echo "Press ENTER when ready to retry this step."
    read

    return 0
}

run_command() {
    local description="$1"
    local logfile="$2"

    shift 2

    while true
    do
        local block_start
        block_start=$(date +%s)

        log_msg "START: $description"

        "$@" > "$logfile" 2>&1

        local status=$?

        if [ $status -eq 0 ]; then
            log_msg "SUCCESS: $description"
            show_runtime "$block_start"
            break
        else
            echo >> "$logfile"
            echo "Exit status: $status" >> "$logfile"
            handle_failure "$description"
        fi
    done
}

###############################################################################
# Activate environment
###############################################################################

source ~/acefa_env/bin/activate

###############################################################################
# No waning
###############################################################################

echo "=== No waning: editing abc_code.py ==="

cd ~/python

sed -i '70c\waning = False' abc_code.py
sed -i '83c\    nr_particles = 15_000' abc_code.py

echo "=== No waning: running abc_code.py ==="

run_command \
    "abc_code.py (no waning)" \
    "$LOG_DIR/abc_no_waning.log" \
    python abc_code.py

echo "=== No waning: editing forward_proj.py ==="

sed -i '26c\local_waning = False' forward_proj.py

echo "=== No waning: running forward_proj.py ==="

run_command \
    "forward_proj.py (no waning)" \
    "$LOG_DIR/forward_no_waning.log" \
    python forward_proj.py

###############################################################################
# With waning
###############################################################################

echo "=== With waning: editing abc_code.py ==="

sed -i '70c\waning = True' abc_code.py
sed -i '83c\    nr_particles = 70_000' abc_code.py

echo "=== With waning: running abc_code.py ==="

run_command \
    "abc_code.py (with waning)" \
    "$LOG_DIR/abc_with_waning.log" \
    python abc_code.py

echo "=== With waning: editing forward_proj.py ==="

sed -i '26c\local_waning = True' forward_proj.py

echo "=== With waning: running forward_proj.py ==="

run_command \
    "forward_proj.py (with waning)" \
    "$LOG_DIR/forward_with_waning.log" \
    python forward_proj.py

###############################################################################
# Leave environment
###############################################################################

deactivate

###############################################################################
# R: no waning
###############################################################################

echo "=== R: no waning ==="

cd ~/R_code

sed -i '2c\source_file_type <- "no waning"' car_ihr_calc_code.R

run_command \
    "car_ihr_calc_code.R (no waning)" \
    "$LOG_DIR/car_ihr_no_waning.log" \
    Rscript car_ihr_calc_code.R

###############################################################################
# R: with waning
###############################################################################

echo "=== R: with waning ==="

cd ~/R_code

sed -i '2c\source_file_type <- "with waning"' car_ihr_calc_code.R

run_command \
    "car_ihr_calc_code.R (with waning)" \
    "$LOG_DIR/car_ihr_with_waning.log" \
    Rscript car_ihr_calc_code.R

###############################################################################
# Finished
###############################################################################

TOTAL_RUNTIME=$(( $(date +%s) - RUN_START ))

echo
echo "======================================================="
echo "ALL DONE"
echo "Total runtime: ${TOTAL_RUNTIME} seconds"
echo "Logs are in: $LOG_DIR"
echo "======================================================="