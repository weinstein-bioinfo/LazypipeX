# TAP-emitting assertion helpers shared by the LazypipeX test tiers.
#
# Source this file; do not execute it.  See docs/testing_roadmap.md §10.
#
#   tap_init "title"        emit the TAP header
#   pass  ID "desc"         report a passing test
#   fail  ID "desc" [diag]  report a failing test (extra args become diagnostics)
#   skipt ID "desc" reason  report a skipped test (counted separately)
#   todof ID "desc" reason  report an expected failure (known defect; TAP TODO)
#   diag  "text"            emit a comment line (multi-line safe)
#   try   cmd args...       run a command under a timeout, capturing $_OUT/$_RC
#   try_sh "shell string"   same, for a command that needs shell parsing
#   tap_done                emit the plan + summary; returns the failure count
#
# Tests never abort the run: pass/fail/skipt always return 0 so a failing check
# cannot short-circuit the rest of the tier.

_TAP_N=0
_TAP_PASS=0
_TAP_FAIL=0
_TAP_SKIP=0
_TAP_TODO=0
_TAP_FAILED_IDS=""

# Seconds any single command may run before it is killed.  Generous, because
# R and some bio-tools are slow to start on a shared filesystem.
: "${TEST_TIMEOUT:=120}"

tap_init() {
	printf 'TAP version 13\n'
	if [ $# -gt 0 ]; then
		diag "$*"
	fi
	return 0
}

# Multi-line safe.  Piped rather than fed by heredoc so that captured command
# output containing $( ) or backticks is never re-expanded.
diag() {
	printf '%s\n' "$*" | while IFS= read -r line; do
		printf '# %s\n' "$line"
	done
	return 0
}

pass() {
	_TAP_N=$(( _TAP_N + 1 ))
	_TAP_PASS=$(( _TAP_PASS + 1 ))
	printf 'ok %d - %s %s\n' "$_TAP_N" "$1" "$2"
	return 0
}

fail() {
	local id="$1" desc="$2"
	shift 2
	_TAP_N=$(( _TAP_N + 1 ))
	_TAP_FAIL=$(( _TAP_FAIL + 1 ))
	_TAP_FAILED_IDS="$_TAP_FAILED_IDS $id"
	printf 'not ok %d - %s %s\n' "$_TAP_N" "$id" "$desc"
	# One diag block per argument: callers pass a captured message and a hint as
	# separate arguments and expect them on separate lines, not run together.
	for _d in "$@"; do
		[ -n "$_d" ] && diag "$_d"
	done
	return 0
}

skipt() {
	_TAP_N=$(( _TAP_N + 1 ))
	_TAP_SKIP=$(( _TAP_SKIP + 1 ))
	printf 'ok %d - %s %s # SKIP %s\n' "$_TAP_N" "$1" "$2" "$3"
	return 0
}

# A check that SHOULD pass but does not, because of a known defect in the code
# under test rather than a broken installation.  Emitted as a TAP TODO so that
# `prove` and CI treat it as an expected failure and the tier's exit status
# still reflects only real, actionable failures.  $3 is the known-drift note.
todof() {
	local id="$1" desc="$2" reason="$3"
	shift 3
	_TAP_N=$(( _TAP_N + 1 ))
	_TAP_TODO=$(( _TAP_TODO + 1 ))
	printf 'not ok %d - %s %s # TODO %s\n' "$_TAP_N" "$id" "$desc" "$reason"
	for _d in "$@"; do
		[ -n "$_d" ] && diag "$_d"
	done
	return 0
}

# Run a command with a timeout.  Combined output lands in $_OUT, exit status in
# $_RC.  Always returns 0 so callers can branch on $_RC instead of on set -e.
try() {
	_OUT=$( timeout "$TEST_TIMEOUT" "$@" 2>&1 )
	_RC=$?
	return 0
}

try_sh() {
	_OUT=$( timeout "$TEST_TIMEOUT" bash -c "$1" 2>&1 )
	_RC=$?
	return 0
}

# Best-effort version string for a tool, for the record in the log.  Tools
# disagree wildly on how to be asked, so the awkward ones are named here and
# everything else falls back to --version.
tool_version() {
	local t="$1" v=""
	case "$t" in
		bwa)          v=$( timeout 20 bwa 2>&1 | grep -i -m1 '^Version' ) ;;
		mga)          v="(no version flag)" ;;
		ktImportText) v=$( timeout 20 ktImportText 2>&1 | grep -i -m1 'krona' ) ;;
		blastn|blastp|blastx|blastdbcmd)
		              v=$( timeout 20 "$t" -version 2>&1 | head -1 ) ;;
		hmmscan|hmmstat)
		              v=$( timeout 20 "$t" -h 2>&1 | grep -i -m1 'HMMER' ) ;;
		seqkit|csvtk|taxonkit)
		              v=$( timeout 20 "$t" version 2>&1 | head -1 ) ;;
		java)         v=$( timeout 20 java -version 2>&1 | head -1 ) ;;
		prodigal)     v=$( timeout 20 prodigal -v 2>&1 | grep -i -m1 'prodigal' ) ;;
		perl)         v=$( timeout 20 perl -e 'print "perl $]"' 2>&1 ) ;;
		*)            v=$( timeout 20 "$t" --version 2>&1 | head -1 ) ;;
	esac
	[ -n "$v" ] || v="(version unknown)"
	printf '%s' "$v"
	return 0
}

tap_done() {
	printf '1..%d\n' "$_TAP_N"
	diag "passed $_TAP_PASS, failed $_TAP_FAIL, skipped $_TAP_SKIP, todo $_TAP_TODO (of $_TAP_N)"
	if [ "$_TAP_FAIL" -gt 0 ]; then
		diag "failed:$_TAP_FAILED_IDS"
	fi
	if [ "$_TAP_FAIL" -gt 255 ]; then
		return 255
	fi
	return "$_TAP_FAIL"
}
