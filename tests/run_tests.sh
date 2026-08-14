#!/usr/bin/env bash
#
# LazypipeX test driver — runs the tiers in order and summarises them.
# See docs/testing_roadmap.md §10 (harness layout) and §11 (when to run what).
#
# Each tier is an independent script emitting TAP on stdout and exiting with its
# failure count.  This driver adds ordering, the blocking rules from §2, and one
# summary table at the end; it does not reinterpret any tier's verdict.
#
# Usage:
#     tests/run_tests.sh                 # every implemented tier, in order
#     tests/run_tests.sh --quick         # tiers 0 and 1 only (no databases, <3 min)
#     tests/run_tests.sh --tier 0,1,2    # an explicit subset
#     tests/run_tests.sh --keep-going    # do not skip tiers behind a failed blocker
#     tests/run_tests.sh --list          # show the tiers and their cost
#     tests/run_tests.sh -q              # summary only; per-tier TAP goes to the log dir
#
# Blocking (docs/testing_roadmap.md §2): a Tier 0 failure makes everything after
# it meaningless, a Tier 1 failure makes Tiers 2-4 meaningless, and a Tier 2
# failure makes Tiers 3-4 meaningless.  Tiers 3 and 4 block nothing.  Rather than
# let one broken environment produce four screens of consequential failures, the
# driver skips what cannot be trusted and says so.  --keep-going overrides.
#
# Exit status: total failed tests across the tiers that ran (0 = all good), or
# 99 for a usage error.  Skips and TODOs are not failures.

set -uo pipefail

TESTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
REPO=$( cd "$TESTS_DIR/.." && pwd )

TIERS_ALL="0 1 2 3 4"
TIERS="$TIERS_ALL"
KEEP_GOING=0
QUIET=0
LOGDIR="${RUN_TESTS_LOGDIR:-${TMPDIR:-/tmp}/lazytest-run.$$}"

tier_script() {
	case "$1" in
		0) printf '%s' "$TESTS_DIR/t0_environment.sh" ;;
		1) printf '%s' "$TESTS_DIR/t1_units.sh" ;;
		2) printf '%s' "$TESTS_DIR/t2_databases.sh" ;;
		3) printf '%s' "$TESTS_DIR/t3_steps.sh" ;;
		4) printf '%s' "$TESTS_DIR/t4_e2e_bio.sh" ;;
		*) printf '' ;;
	esac
}

tier_name() {
	case "$1" in
		0) printf '%s' "environment and prerequisites" ;;
		1) printf '%s' "unit and smoke tests" ;;
		2) printf '%s' "database installation" ;;
		3) printf '%s' "step-wise pipeline" ;;
		4) printf '%s' "end-to-end and biological" ;;
	esac
}

tier_cost() {
	case "$1" in
		0) printf '%s' "no databases, <30 s" ;;
		1) printf '%s' "no databases, <2 min" ;;
		2) printf '%s' "databases, 1-3 min" ;;
		3) printf '%s' "databases, 3-5 min" ;;
		4) printf '%s' "databases, 15-25 min" ;;
	esac
}

# Which tiers a failure here invalidates (§2's Blocking column).
tier_blocks() {
	case "$1" in
		0) printf '%s' "1 2 3 4" ;;
		1) printf '%s' "2 3 4" ;;
		2) printf '%s' "3 4" ;;
		*) printf '' ;;
	esac
}

usage() {
	sed -n '2,26p' "$0" | sed 's/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
	case "$1" in
		--tier|-t)
			shift
			[ $# -gt 0 ] || { echo "ERROR: --tier needs a list, e.g. --tier 0,1,2" >&2; exit 99; }
			TIERS=$( printf '%s' "$1" | tr ',' ' ' )
			;;
		--quick)      TIERS="0 1" ;;
		--keep-going) KEEP_GOING=1 ;;
		-q|--quiet)   QUIET=1 ;;
		--list)
			printf 'tier  script                cost                    name\n'
			for t in $TIERS_ALL; do
				s=$( tier_script "$t" )
				printf '%-5s %-21s %-23s %s%s\n' "$t" "$( basename "$s" )" \
					"$( tier_cost "$t" )" "$( tier_name "$t" )" \
					"$( [ -x "$s" ] || [ -f "$s" ] && printf '' || printf '   [MISSING]' )"
			done
			exit 0
			;;
		-h|--help) usage; exit 0 ;;
		*) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 99 ;;
	esac
	shift
done

for t in $TIERS; do
	case " $TIERS_ALL " in
		*" $t "*) ;;
		*) echo "ERROR: no such tier: $t (have: $TIERS_ALL)" >&2; exit 99 ;;
	esac
done

mkdir -p "$LOGDIR" || exit 99

printf '# LazypipeX test run\n'
printf '# host    : %s\n' "$( hostname )"
printf '# date    : %s\n' "$( date -Is )"
printf '# repo    : %s\n' "$REPO"
printf '# tiers   : %s\n' "$TIERS"
printf '# logs    : %s\n' "$LOGDIR"
printf '#\n'

TOTAL_FAIL=0
SUMMARY=""
SKIP_TIERS=""

for t in $TIERS; do
	script=$( tier_script "$t" )
	label="tier $t ($( tier_name "$t" ))"

	if [ ! -f "$script" ]; then
		SUMMARY="$SUMMARY$( printf '%-5s %-32s %s' "$t" "$( tier_name "$t" )" "NOT IMPLEMENTED" )\n"
		printf '# --- %s: not implemented, skipping\n' "$label"
		continue
	fi

	case " $SKIP_TIERS " in
		*" $t "*)
			SUMMARY="$SUMMARY$( printf '%-5s %-32s %s' "$t" "$( tier_name "$t" )" "SKIPPED (blocker failed)" )\n"
			printf '# --- %s: skipped, an earlier blocking tier failed\n' "$label"
			continue
			;;
	esac

	printf '# === %s: %s\n' "$label" "$( tier_cost "$t" )"
	log="$LOGDIR/tier$t.tap"
	t0=$( date +%s )
	bash "$script" > "$log" 2>&1
	rc=$?
	secs=$(( $( date +%s ) - t0 ))

	[ "$QUIET" -eq 1 ] || sed 's/^/  /' "$log"

	# Every tier ends with "# passed N, failed N, skipped N, todo N (of N)".
	line=$( grep -m1 '^# passed ' "$log" )
	if [ -z "$line" ]; then
		# The tier died before its summary — a broken harness, not a test verdict.
		SUMMARY="$SUMMARY$( printf '%-5s %-32s %s' "$t" "$( tier_name "$t" )" "ERROR (no summary; rc=$rc, ${secs}s)" )\n"
		printf '# --- %s: produced no summary line (rc=%s); see %s\n' "$label" "$rc" "$log"
		TOTAL_FAIL=$(( TOTAL_FAIL + 1 ))
		for b in $( tier_blocks "$t" ); do SKIP_TIERS="$SKIP_TIERS $b"; done
		continue
	fi

	stats=$( printf '%s' "$line" | sed 's/^# //' )
	nfail=$( printf '%s' "$line" | sed -n 's/.*failed \([0-9]*\).*/\1/p' )
	nfail=${nfail:-0}
	TOTAL_FAIL=$(( TOTAL_FAIL + nfail ))
	SUMMARY="$SUMMARY$( printf '%-5s %-32s %s (%ss)' "$t" "$( tier_name "$t" )" "$stats" "$secs" )\n"

	if [ "$nfail" -gt 0 ] && [ "$KEEP_GOING" -eq 0 ]; then
		for b in $( tier_blocks "$t" ); do
			case " $TIERS " in *" $b "*) SKIP_TIERS="$SKIP_TIERS $b" ;; esac
		done
		[ -n "$( tier_blocks "$t" )" ] && \
			printf '# --- %s failed; skipping tiers %s (--keep-going to run anyway)\n' \
				"$label" "$( tier_blocks "$t" )"
	fi
	printf '#\n'
done

printf '#\n'
printf '# ================================ summary ================================\n'
printf '%b' "$SUMMARY" | sed 's/^/# /'
printf '# %s\n' "-----------------------------------------------------------------------"
if [ "$TOTAL_FAIL" -eq 0 ]; then
	printf '# ALL GOOD — 0 failed tests\n'
else
	printf '# %d failed test(s); per-tier TAP in %s\n' "$TOTAL_FAIL" "$LOGDIR"
fi

[ "$TOTAL_FAIL" -gt 255 ] && exit 255
exit "$TOTAL_FAIL"
