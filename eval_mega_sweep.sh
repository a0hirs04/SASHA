#!/usr/bin/env bash
# Evaluate all 100 mega sweep variants using official evaluators
# Reports: RC1 C1.4 status + RC2 all criteria for each variant
set -uo pipefail

SWEEP="/work/a0hirs04/PROJECT-NORTHSTAR/build/mega_sweep"
PROJECT="/home/a0hirs04/PROJECT-NORTHSTAR"

echo "════════════════════════════════════════════════════════════════"
echo "  MEGA SWEEP EVALUATION — 100 VARIANTS × (RC1 + RC2)"
echo "  $(date)"
echo "════════════════════════════════════════════════════════════════"

# Header
printf "\n%-6s  %-6s %-6s %-6s %-6s %-6s %-6s  " \
    "Var" "uptake" "kill" "hif" "cap" "off" "tgfb"
printf "%-5s %-5s %-5s %-5s %-5s %-5s %-5s %-5s  " \
    "C1.1" "C1.2" "C1.3" "C1.4" "C1.5" "C1.6" "C1.7" "C1.8"
printf "%-5s %-5s %-5s %-5s %-5s %-5s\n" \
    "C2.1" "C2.2" "C2.3" "C2.4" "C2.5" "C2.6"
printf '%0.s─' {1..120}
echo ""

# Track winners
winners=""
rc1_pass_count=0
rc2_pass_count=0
both_pass_count=0

for v in $(seq -w 1 100); do
    vname="v${v}"
    vdir="$SWEEP/$vname"
    [[ -d "$vdir" ]] || continue

    # Get params
    uptake="?"  kill="?"  hif="?"  cap="?"  off="?"  tgfb="?"
    if [[ -f "$SWEEP/variants.tsv" ]]; then
        line=$(grep "^${vname}" "$SWEEP/variants.tsv" 2>/dev/null | head -1)
        if [[ -n "$line" ]]; then
            uptake=$(echo "$line" | cut -f2)
            kill=$(echo "$line" | cut -f3)
            hif=$(echo "$line" | cut -f4)
            cap=$(echo "$line" | cut -f5)
            off=$(echo "$line" | cut -f6)
            tgfb=$(echo "$line" | cut -f7)
        fi
    fi

    printf "%-6s  %-6s %-6s %-6s %-6s %-6s %-6s  " \
        "$vname" "$uptake" "$kill" "$hif" "$cap" "$off" "$tgfb"

    # ── RC1 evaluation ──
    rc1_work="$vdir/rc1"
    rc1_out="$rc1_work/replicate_01_seed42/output"
    rc1_result=""
    rc1_all_pass=true

    if [[ -d "$rc1_out" ]] && ls "$rc1_out"/output*.xml &>/dev/null; then
        rc1_result=$(python3 "$PROJECT/run_reality_check_1.py" \
            --work-dir "$rc1_work" --seeds 42 --quorum 1 --evaluate-only 2>&1)

        for crit in "C1.1" "C1.2" "C1.3" "C1.4" "C1.5" "C1.6" "C1.7" "C1.8"; do
            if echo "$rc1_result" | grep -q "\[PASS\].*${crit}\|${crit}.*\[PASS\]"; then
                printf "%-5s " "PASS"
            elif echo "$rc1_result" | grep -q "\[FAIL\].*${crit}\|${crit}.*\[FAIL\]"; then
                printf "%-5s " "FAIL"
                rc1_all_pass=false
            else
                # Try alternate grep
                pass_count=$(echo "$rc1_result" | grep -c "\[PASS\]" || true)
                fail_count=$(echo "$rc1_result" | grep -c "\[FAIL\]" || true)
                printf "%-5s " "?"
            fi
        done
    else
        printf "%-5s %-5s %-5s %-5s %-5s %-5s %-5s %-5s  " \
            "--" "--" "--" "--" "--" "--" "--" "--"
        rc1_all_pass=false
    fi

    # ── RC2 evaluation ──
    rc2_out="$vdir/rc2/output"
    rc2_all_pass=true

    if [[ -d "$rc2_out" ]] && ls "$rc2_out"/output*.xml &>/dev/null; then
        rc2_result=$(python3 "$PROJECT/evaluate_rc2.py" \
            --out-dir "$rc2_out" 2>&1)

        for crit in "RC2-1" "RC2-2" "RC2-3" "RC2-4" "RC2-5" "RC2-6"; do
            if echo "$rc2_result" | grep -q "\[PASS\].*${crit}\|${crit}.*\[PASS\]"; then
                printf "%-5s " "PASS"
            elif echo "$rc2_result" | grep -q "\[FAIL\].*${crit}\|${crit}.*\[FAIL\]"; then
                printf "%-5s " "FAIL"
                rc2_all_pass=false
            elif echo "$rc2_result" | grep -q "\[INFO\].*${crit}\|${crit}.*\[INFO\]"; then
                printf "%-5s " "INFO"
            else
                printf "%-5s " "?"
            fi
        done
    else
        printf "%-5s %-5s %-5s %-5s %-5s %-5s" \
            "--" "--" "--" "--" "--" "--"
        rc2_all_pass=false
    fi

    echo ""

    # Track
    $rc1_all_pass && rc1_pass_count=$((rc1_pass_count + 1))
    $rc2_all_pass && rc2_pass_count=$((rc2_pass_count + 1))
    if $rc1_all_pass && $rc2_all_pass; then
        both_pass_count=$((both_pass_count + 1))
        winners="$winners $vname"
    fi
done

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  SUMMARY"
echo "════════════════════════════════════════════════════════════════"
echo "  RC1 all-pass: $rc1_pass_count / 100"
echo "  RC2 all-pass: $rc2_pass_count / 100"
echo "  BOTH pass:    $both_pass_count / 100"
if [[ -n "$winners" ]]; then
    echo ""
    echo "  WINNERS (pass RC1 + RC2):$winners"
    echo ""
    echo "  Next: Run full 5-seed gate on winners, then RC3-RC5"
fi
echo "════════════════════════════════════════════════════════════════"
