#!/usr/bin/env bash
#
# LazypipeX Tier 1 — unit and smoke tests.
# Implements UNIT-01 … UNIT-32 of docs/testing_roadmap.md §4.
#
# Needs no reference databases and no network: everything here is either a
# manual/argument-parsing check, a pure-Perl unit test, or a compiled helper run
# against a fixture under tests/fixtures/.
#
# Usage:
#     module use /projappl/project_2003755/Lazypipe-db/modulefiles/projects
#     module load lazypipe/3.1
#     tests/t1_units.sh                # TAP on stdout
#     tests/t1_units.sh | grep -v ^#   # results only
#
# Exit status = number of failed tests (0 = all good).  Neither skips nor TODOs
# are failures: skips mark absent optional tooling, TODOs mark known defects in
# LazypipeX itself that are documented in docs/testing_roadmap.md §12.
#
# The script writes nothing outside $TMPDIR — in particular UNIT-30 compiles to
# the scratch directory rather than over the repo's bin/.

set -uo pipefail

TESTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# shellcheck source=lib/assert.sh
. "$TESTS_DIR/lib/assert.sh"

REPO=$( cd "$TESTS_DIR/.." && pwd )

# Mirror lazypipe.pl's own install-dir resolution so the tests exercise the
# installation the user will actually run (see t0_environment.sh).
if [ -n "${LAZYPIPE_INSTALL_DIR:-}" ] && [ -f "$LAZYPIPE_INSTALL_DIR/lazypipe.pl" ]; then
	INSTALL="$LAZYPIPE_INSTALL_DIR"
else
	INSTALL="$REPO"
fi

FIX="$TESTS_DIR/fixtures"
LZ="$INSTALL/lazypipe.pl"

WORK=$( mktemp -d "${TMPDIR:-/tmp}/lazytest-t1.XXXXXX" ) || exit 99
trap 'rm -rf "$WORK"' EXIT

# lazypipe.pl reads config.yaml from the CURRENT directory if present, so run
# every invocation from the install dir for a predictable config.
cd "$INSTALL" || exit 99

# The shipped tmpdir default is $LOCAL_SCRATCH/wrkdir, which is unset outside a
# Slurm job and then resolves to /wrkdir.  Every invocation below therefore
# passes --tmpdir explicitly; that is a test-harness decision, not a workaround
# for a bug the tier should hide (ENV-09 and HPC-03 cover the default itself).
TMPD="$WORK/tmp"
mkdir -p "$TMPD"

# A tool counts as available only if it is also executable.  `command -v`
# alone is not enough: bash returns the path of a non-executable file found
# on PATH, so a downloaded-but-not-chmod+x binary was reported as installed
# and then failed at run time with "Permission denied".
have() {
	local p
	p=$( type -P "$1" 2>/dev/null ) || return 1
	[ -n "$p" ] && [ -x "$p" ]
}

# Snapshot for UNIT-99: whatever the working tree looks like before any test
# runs is the baseline this tier must leave untouched.
TREE_BEFORE=$( git -C "$REPO" status --porcelain 2>/dev/null )

# Number of lines in $_OUT matching a pattern.
count_matches() { printf '%s\n' "$_OUT" | grep -c -- "$1"; }

# ----------------------------------------------------------------- report ---

tap_init "LazypipeX Tier 1 — unit and smoke tests"
diag "host        : $( hostname )"
diag "date        : $( date -Is )"
diag "repo        : $REPO"
diag "install dir : $INSTALL"
diag "scratch     : $WORK"
if [ "$INSTALL" != "$REPO" ]; then
	diag "NOTE: \$LAZYPIPE_INSTALL_DIR differs from this checkout; testing the installed copy."
fi
diag ""

# ================================================== UNIT-01 main manual ======

# The 11 --pipe step names the manual documents.  UNIT-06 checks that the
# parser agrees with this list.
PIPE_STEPS="pre flt ass rea ann1 ann2 rep rgrep sta pack clean"

try perl "$LZ" -h
if [ "$_RC" -ne 0 ]; then
	fail UNIT-01 "main manual prints" "perl lazypipe.pl -h exited $_RC" "$_OUT"
else
	u01_missing=""
	printf '%s' "$_OUT" | grep -q 'USAGE' || u01_missing="$u01_missing USAGE"
	for s in $PIPE_STEPS; do
		# Steps appear in the manual as "pre|preprocess :" / "pack           :".
		printf '%s' "$_OUT" | grep -Eq "^ +($s\||$s +:)" || u01_missing="$u01_missing $s"
	done
	if [ -z "$u01_missing" ]; then
		pass UNIT-01 "main manual prints USAGE and all 11 --pipe step names"
	else
		fail UNIT-01 "main manual prints USAGE and all 11 --pipe step names" \
			"absent from the manual:$u01_missing"
	fi
fi

# The roadmap expects the manual to carry the version.  It does not: the string
# lives only in $PIPELINE_VERSION and is never printed, so a user cannot tell
# which LazypipeX they are running from the manual alone.
try perl "$LZ" -h
if printf '%s' "$_OUT" | grep -q '3\.1'; then
	pass UNIT-01a "main manual names the pipeline version"
else
	todof UNIT-01a "main manual names the pipeline version" \
		"lazypipe.pl never prints \$PIPELINE_VERSION" \
		"add the version to \$usage; there is currently no way to read it off the CLI"
fi

# ============================================== UNIT-02 install_db manual ====

try perl "$INSTALL/perl/install_db.pl" -h
if printf '%s' "$_OUT" | grep -qi 'usage'; then
	pass UNIT-02 "install_db.pl prints its manual"
else
	fail UNIT-02 "install_db.pl prints its manual" "rc=$_RC, no usage block" "$_OUT"
fi

# ============================================ UNIT-03 helper script manuals ==

# Every shipped script must answer -h with a usage block rather than a stack
# trace.  Compiled helpers print usage with no arguments and are covered by
# ENV-10; here they are re-checked for the usage block itself.
u03_bad=""
u03_n=0
for f in "$INSTALL"/perl/*.pl "$INSTALL"/scripts/*.pl; do
	[ -e "$f" ] || continue
	u03_n=$(( u03_n + 1 ))
	rel="${f#$INSTALL/}"
	try perl "$f" -h
	if printf '%s' "$_OUT" | grep -qi 'usage'; then
		continue
	fi
	if printf '%s' "$_OUT" | grep -q "Can't locate"; then
		miss=$( printf '%s' "$_OUT" | sed -n "s|.*Can't locate \([A-Za-z0-9_/]*\)\.pm.*|\1|p" \
			| head -1 | sed 's|/|::|g' )
		u03_bad="$u03_bad\n  $rel: does not load — missing Perl module $miss"
	elif printf '%s' "$_OUT" | grep -q "BEGIN failed\|syntax error"; then
		u03_bad="$u03_bad\n  $rel: does not load ($( printf '%s' "$_OUT" | head -1 | cut -c1-90 ))"
	else
		u03_bad="$u03_bad\n  $rel: rc=$_RC, no usage block"
	fi
done
for b in retrieve_reads get_contigs filtfa filtfq; do
	[ -x "$INSTALL/bin/$b" ] || continue
	u03_n=$(( u03_n + 1 ))
	try "$INSTALL/bin/$b"
	printf '%s' "$_OUT" | grep -qi 'usage' \
		|| u03_bad="$u03_bad\n  bin/$b: rc=$_RC, no usage block"
done
if [ -z "$u03_bad" ]; then
	pass UNIT-03 "all $u03_n shipped scripts and helpers print a usage block"
else
	fail UNIT-03 "all $u03_n shipped scripts and helpers print a usage block" \
		"$( printf '%b' "$u03_bad" )"
fi

# ================================================= UNIT-04 compile-only ======

u04_bad=""
u04_n=0
for f in "$LZ" "$INSTALL"/perl/*.pl "$INSTALL"/scripts/*.pl "$INSTALL"/perl/Lazypipe/*.pm; do
	[ -e "$f" ] || continue
	u04_n=$(( u04_n + 1 ))
	rel="${f#$INSTALL/}"
	try perl -I "$INSTALL/perl" -c "$f"
	if [ "$_RC" -ne 0 ]; then
		# List every module the file needs and cannot find, not just the first:
		# fixing one only to rerun and hit the next wastes a cycle each time.
		missing=$( perl -e '
			my %seen;
			for my $m ( $ARGV[0] =~ /Can.t locate ([A-Za-z0-9_\/]+)\.pm/g ) {
				$m =~ s{/}{::}g;
				$seen{$m} = 1;
			}
			print join( ", ", sort keys %seen );
		' "$_OUT" )
		if [ -n "$missing" ]; then
			# Re-probe with each failing use-line stubbed out, to report every
			# missing module rather than only the one that aborted compilation.
			allmiss=$( perl -ne 'print "$1\n" if /^\s*use\s+([A-Za-z0-9_:]+)/' "$f" \
				| while IFS= read -r m; do
					perl -M"$m" -e1 >/dev/null 2>&1 || printf '%s ' "$m"
				done )
			u04_bad="$u04_bad\n  $rel: missing Perl modules: ${allmiss:-$missing}"
		else
			u04_bad="$u04_bad\n  $rel: $( printf '%s' "$_OUT" | grep -v '^$' | head -1 | cut -c1-100 )"
		fi
	fi
done
if [ -z "$u04_bad" ]; then
	pass UNIT-04 "all $u04_n Perl files compile (syntax OK)"
else
	fail UNIT-04 "all $u04_n Perl files compile (syntax OK)" \
		"$( printf '%b' "$u04_bad" )" \
		"a file that will not compile cannot run, whichever --pipe step reaches it"
fi

# ========================================== UNIT-05 unknown option rejected ==

# Getopt::Long::Configure("pass_through") is set for the config-file pre-scan at
# lazypipe.pl:113 and never reset, so it is still in force for the main
# GetOptions call.  An unknown option is therefore silently ignored instead of
# being rejected.  The run below is otherwise completely valid, so anything but
# a clean exit-0 would be the unknown option being noticed.
mkdir -p "$WORK/res05"
try_sh "perl '$LZ' --se -1 '$INSTALL/data/samples/M15small_R1.fastq' --nosuchoption \
	-p clean --res '$WORK/res05' --tmpdir '$TMPD' -s u05"
if [ "$_RC" -ne 0 ] && printf '%s' "$_OUT" | grep -qi 'unknown option\|USAGE'; then
	pass UNIT-05 "unknown command-line option is rejected"
else
	todof UNIT-05 "unknown command-line option is rejected" \
		"pass_through leaks from the --config pre-scan into the main GetOptions" \
		"--nosuchoption was accepted (rc=$_RC, no usage printed)" \
		"a mistyped option silently runs with defaults; see docs/testing_roadmap.md §12"
fi

# ============================================ UNIT-06 invalid pipe step ======

try_sh "perl '$LZ' -1 '$INSTALL/data/samples/M15small_R1.fastq' -p bogus \
	--res '$WORK/res06' --tmpdir '$TMPD'"
if [ "$_RC" -ne 0 ] && printf '%s' "$_OUT" | grep -q 'invalid argument --pipe'; then
	pass UNIT-06 "invalid --pipe step is rejected"
else
	fail UNIT-06 "invalid --pipe step is rejected" "rc=$_RC" "$_OUT"
fi

# Every documented step name must be accepted by the parser.  Run each one
# through the same early-exit probe used by UNIT-09: a nonexistent read1 with
# -p pre appended makes the run die *after* --pipe has been parsed, so a step
# name the parser does not know still surfaces as 'invalid argument --pipe'.
u06_bad=""
for s in $PIPE_STEPS main all; do
	try_sh "perl '$LZ' --se -1 /nonexistent_R1.fastq -p '$s,pre' \
		--res '$WORK/res06' --tmpdir '$TMPD'"
	printf '%s' "$_OUT" | grep -q 'invalid argument --pipe' \
		&& u06_bad="$u06_bad $s"
done
if [ -z "$u06_bad" ]; then
	pass UNIT-06a "every documented --pipe step name is accepted by the parser"
else
	fail UNIT-06a "every documented --pipe step name is accepted by the parser" \
		"rejected:$u06_bad" \
		"the manual documents a step the parser does not implement"
fi

# ============================================== UNIT-07 missing reads ========

try_sh "perl '$LZ' -1 /nonexistent_R1.fastq -2 /nonexistent_R2.fastq -p pre \
	--res '$WORK/res07' --tmpdir '$TMPD'"
if [ "$_RC" -ne 0 ] && printf '%s' "$_OUT" | grep -q 'check read1 file'; then
	pass UNIT-07 "missing read1 is rejected with a clear message"
else
	fail UNIT-07 "missing read1 is rejected with a clear message" "rc=$_RC" "$_OUT"
fi

# ========================================== UNIT-08 reverse-read guessing ====

# _R1 -> _R2.  The guessed sibling deliberately does not exist, so the run dies
# naming the file it guessed — which is what proves the substitution happened.
cp "$INSTALL/data/samples/M15small_R1.fastq" "$WORK/guess_R1.fastq"
try_sh "perl '$LZ' -1 '$WORK/guess_R1.fastq' -p pre --res '$WORK/res08' --tmpdir '$TMPD'"
if printf '%s' "$_OUT" | grep -q "check read2 file: $WORK/guess_R2.fastq"; then
	u08_R1=ok
else
	u08_R1="did not resolve _R1 -> _R2 (rc=$_RC)"
fi

# lower-case _r1 -> _r2
cp "$INSTALL/data/samples/M15small_R1.fastq" "$WORK/guess_r1.fastq"
try_sh "perl '$LZ' -1 '$WORK/guess_r1.fastq' -p pre --res '$WORK/res08' --tmpdir '$TMPD'"
if printf '%s' "$_OUT" | grep -q "check read2 file: $WORK/guess_r2.fastq"; then
	u08_r1=ok
else
	u08_r1="did not resolve _r1 -> _r2 (rc=$_RC)"
fi

# A name with no recognisable mate marker must die with a clear message rather
# than silently reusing read1 as read2.
cp "$INSTALL/data/samples/M15small_R1.fastq" "$WORK/sample.fastq"
try_sh "perl '$LZ' -1 '$WORK/sample.fastq' -p pre --res '$WORK/res08' --tmpdir '$TMPD'"
if [ "$_RC" -ne 0 ] && printf '%s' "$_OUT" | grep -qi 'filename for reverse reads'; then
	u08_bad_name=ok
else
	u08_bad_name="unguessable name did not die cleanly (rc=$_RC)"
fi

if [ "$u08_R1" = ok ] && [ "$u08_r1" = ok ] && [ "$u08_bad_name" = ok ]; then
	pass UNIT-08 "reverse-read filename guessing (_R1/_r1) and its failure mode"
else
	fail UNIT-08 "reverse-read filename guessing (_R1/_r1) and its failure mode" \
		"$( [ "$u08_R1"       = ok ] || printf '  %s\n' "$u08_R1" )" \
		"$( [ "$u08_r1"       = ok ] || printf '  %s\n' "$u08_r1" )" \
		"$( [ "$u08_bad_name" = ok ] || printf '  %s\n' "$u08_bad_name" )"
fi

# ======================================= UNIT-09 annotation strategies =======

# Dry probe: options_format() resolves --anns *before* it validates read1, so a
# nonexistent read1 with -p pre exercises strategy resolution and then exits
# without running any search.  A key that fails to resolve prints
# 'undefined annotation strategy'; a key that resolves to a malformed value
# dies with 'invalid annotation.strategy'.
STRATEGIES=$( perl -MYAML::Tiny -e '
	my $y = YAML::Tiny->read($ARGV[0]) or exit 1;
	my $s = $y->[0]->{"ann.strategies"} or exit 0;
	print "$_\n" for sort keys %$s;
' "$INSTALL/config.yaml" 2>/dev/null )

if [ -z "$STRATEGIES" ]; then
	fail UNIT-09 "every ann.strategies key resolves" "no ann.strategies section in config.yaml"
else
	u09_bad=""
	u09_n=0
	while IFS= read -r key; do
		[ -n "$key" ] || continue
		u09_n=$(( u09_n + 1 ))
		try_sh "perl '$LZ' --anns '$key' -1 /nonexistent_R1.fastq -p pre \
			--res '$WORK/res09' --tmpdir '$TMPD'"
		if printf '%s' "$_OUT" | grep -q 'undefined annotation strategy'; then
			u09_bad="$u09_bad\n  $key: not resolved (undefined annotation strategy)"
		elif printf '%s' "$_OUT" | grep -q 'invalid annotation.strategy'; then
			u09_bad="$u09_bad\n  $key: resolved to a value with no --ann1"
		fi
	done <<< "$STRATEGIES"
	if [ -z "$u09_bad" ]; then
		pass UNIT-09 "all $u09_n ann.strategies keys resolve to a non-empty --ann1"
	else
		fail UNIT-09 "all $u09_n ann.strategies keys resolve to a non-empty --ann1" \
			"$( printf '%b' "$u09_bad" )" \
			"a key that does not resolve leaves ann1 at the config default and the run produces no annotation"
	fi
fi

# ======================================= UNIT-10 unknown strategy is loud ====

# An unknown --anns key only prints to STDERR and lets the run continue with
# ann1 at its config default (0), so pipe_annotation_round1 is skipped and the
# run "succeeds" with no annotation at all.  Probed the same dry way, so the
# exit status below comes from the read1 check; what is being tested is whether
# the strategy error is fatal on its own.
try_sh "perl '$LZ' --anns no.such.strategy --se -1 '$INSTALL/data/samples/M15small_R1.fastq' \
	-p clean --res '$WORK/res10' --tmpdir '$TMPD' -s u10"
if [ "$_RC" -ne 0 ]; then
	pass UNIT-10 "an unknown --anns key is fatal"
else
	todof UNIT-10 "an unknown --anns key is fatal" \
		"options_format() warns on STDERR and continues" \
		"the run exited 0 with no annotation performed; see docs/testing_roadmap.md §12 item 2"
fi

# ============================================= UNIT-11 History logging =======

mkdir -p "$WORK/res11"
try_sh "perl '$LZ' --se -1 '$INSTALL/data/samples/M15small_R1.fastq' -p clean \
	--res '$WORK/res11' --tmpdir '$TMPD' -s u11"
HIST="$WORK/res11/u11/History.log"
if [ "$_RC" -ne 0 ]; then
	fail UNIT-11 "each invocation appends to History.log" "the probe run itself failed (rc=$_RC)" "$_OUT"
elif [ ! -s "$HIST" ]; then
	fail UNIT-11 "each invocation appends to History.log" "no History.log at $HIST"
else
	n1=$( wc -l < "$HIST" )
	try_sh "perl '$LZ' --se -1 '$INSTALL/data/samples/M15small_R1.fastq' -p clean \
		--res '$WORK/res11' --tmpdir '$TMPD' -s u11"
	n2=$( wc -l < "$HIST" )
	u11_bad=""
	[ "$n2" -gt "$n1" ] || u11_bad="$u11_bad\n  second run did not append (was $n1 lines, now $n2)"
	grep -q '[0-9]\{4\}/[0-9]\{2\}/[0-9]\{2\}' "$HIST" || u11_bad="$u11_bad\n  no timestamp on the logged line"
	grep -q 'lazypipe.pl' "$HIST" || u11_bad="$u11_bad\n  no command line on the logged line"
	if [ -z "$u11_bad" ]; then
		pass UNIT-11 "each invocation appends a timestamped command line to History.log"
	else
		fail UNIT-11 "each invocation appends a timestamped command line to History.log" \
			"$( printf '%b' "$u11_bad" )"
	fi
fi

# ======================================= UNIT-20..25 library unit tests ======

# The Test::More suites carry their own assertions; this tier reports one TAP
# line per file plus a roll-up, and forwards the failing subtest names.
if ! have prove; then
	skipt UNIT-20 "Lazypipe::Utils unit tests" "prove not on PATH"
	skipt UNIT-25 "CIGAR helper unit tests"    "prove not on PATH"
else
	for t in utils cigar; do
		tfile="$TESTS_DIR/perl/$t.t"
		case "$t" in
			utils) tid=UNIT-20; tdesc="Lazypipe::Utils unit tests (UNIT-20..24)" ;;
			cigar) tid=UNIT-25; tdesc="CIGAR helper unit tests (UNIT-25)" ;;
		esac
		if [ ! -f "$tfile" ]; then
			fail "$tid" "$tdesc" "missing test file: $tfile"
			continue
		fi
		try perl -I "$INSTALL/perl" "$tfile"
		nok=$( printf '%s\n' "$_OUT" | grep -c '^ok ' )
		nnok=$( printf '%s\n' "$_OUT" | grep '^not ok ' | grep -vc '# TODO' )
		ntodo=$( printf '%s\n' "$_OUT" | grep -c '^not ok .*# TODO' )
		if [ "$_RC" -eq 0 ] && [ "$nnok" -eq 0 ]; then
			pass "$tid" "$tdesc — $nok assertions, $ntodo todo"
		else
			fail "$tid" "$tdesc" \
				"$nnok failing assertions (rc=$_RC):" \
				"$( printf '%s\n' "$_OUT" | grep '^not ok ' | grep -v '# TODO' | head -8 )"
		fi
	done
fi

# ========================================= UNIT-30 compiled helper build =====

# Compiled to $WORK, never over the repo's bin/: a test run must not change the
# binaries the user is testing.
if ! have g++; then
	skipt UNIT-30 "retrieve_reads compiles" "no g++ on PATH"
else
	try_sh "g++ -Wall -O3 -std=c++11 -I'$INSTALL/cpp' '$INSTALL/cpp/retrieve_reads.cpp' \
		-o '$WORK/retrieve_reads' 2>&1"
	if [ "$_RC" -eq 0 ] && [ -x "$WORK/retrieve_reads" ]; then
		nwarn=$( printf '%s\n' "$_OUT" | grep -c 'warning:' )
		pass UNIT-30 "retrieve_reads compiles from cpp/ ($nwarn compiler warnings)"
	else
		fail UNIT-30 "retrieve_reads compiles from cpp/" "rc=$_RC" "$_OUT"
	fi
fi

# The SeqAn-based helpers need the SeqAn headers, which are not part of this
# repo; they are built out-of-band and only smoke-tested here (UNIT-31).
if ! have g++; then
	skipt UNIT-30a "SeqAn-based helpers compile" "no g++ on PATH"
elif [ -z "${seqan:-}" ]; then
	skipt UNIT-30a "SeqAn-based helpers compile" "\$seqan not set (see Makefile)"
else
	u30_bad=""
	for h in get_contigs filtfa filtfq; do
		try_sh "g++ -Wall -O3 -DNDEBUG -std=c++14 -I'$INSTALL/cpp' -I'$seqan/include' \
			'$INSTALL/cpp/$h.cpp' -o '$WORK/$h' 2>&1"
		[ "$_RC" -eq 0 ] || u30_bad="$u30_bad\n  $h: $( printf '%s' "$_OUT" | head -1 | cut -c1-90 )"
	done
	if [ -z "$u30_bad" ]; then
		pass UNIT-30a "SeqAn-based helpers compile"
	else
		fail UNIT-30a "SeqAn-based helpers compile" "$( printf '%b' "$u30_bad" )"
	fi
fi

# ========================================== UNIT-31 filtfa / filtfq smoke ====

# 10 records in, 3 ids listed: filter mode must keep 7, select mode exactly the
# 3 named — and the two modes must partition the input with nothing lost.
if [ ! -x "$INSTALL/bin/filtfa" ]; then
	skipt UNIT-31 "filtfa filters and selects by id" "bin/filtfa not built"
else
	try_sh "'$INSTALL/bin/filtfa' -i '$FIX/fa/seqs.fa' -o '$WORK/fa.filter.fa' \
		-f '$FIX/fa/ids.txt' -m filter"
	rc_f=$_RC
	try_sh "'$INSTALL/bin/filtfa' -i '$FIX/fa/seqs.fa' -o '$WORK/fa.select.fa' \
		-f '$FIX/fa/ids.txt' -m select"
	rc_s=$_RC
	nf=$( grep -c '^>' "$WORK/fa.filter.fa" 2>/dev/null || echo 0 )
	ns=$( grep -c '^>' "$WORK/fa.select.fa" 2>/dev/null || echo 0 )
	u31_bad=""
	[ "$rc_f" -eq 0 ] || u31_bad="$u31_bad\n  filter mode exited $rc_f"
	[ "$rc_s" -eq 0 ] || u31_bad="$u31_bad\n  select mode exited $rc_s"
	[ "$nf" -eq 7 ]   || u31_bad="$u31_bad\n  filter kept $nf records, expected 7"
	[ "$ns" -eq 3 ]   || u31_bad="$u31_bad\n  select kept $ns records, expected 3"
	if [ -s "$WORK/fa.select.fa" ]; then
		got=$( grep '^>' "$WORK/fa.select.fa" | sed 's/^>//; s/ .*//' | sort | tr '\n' ',' )
		[ "$got" = "c2,c5,c9," ] || u31_bad="$u31_bad\n  select returned [$got], expected [c2,c5,c9,]"
	fi
	if [ -z "$u31_bad" ]; then
		pass UNIT-31 "filtfa: filter keeps 7/10, select keeps the 3 named records"
	else
		fail UNIT-31 "filtfa: filter keeps 7/10, select keeps the 3 named records" \
			"$( printf '%b' "$u31_bad" )"
	fi
fi

if [ ! -x "$INSTALL/bin/filtfq" ]; then
	skipt UNIT-31a "filtfq filters and selects read pairs" "bin/filtfq not built"
else
	try_sh "'$INSTALL/bin/filtfq' -1 '$FIX/fq/reads_R1.fastq' -2 '$FIX/fq/reads_R2.fastq' \
		-o '$WORK/fq.flt_R1.fastq' -O '$WORK/fq.flt_R2.fastq' -f '$FIX/fq/ids.txt' -m filter"
	rc_f=$_RC
	try_sh "'$INSTALL/bin/filtfq' -1 '$FIX/fq/reads_R1.fastq' -2 '$FIX/fq/reads_R2.fastq' \
		-o '$WORK/fq.sel_R1.fastq' -O '$WORK/fq.sel_R2.fastq' -f '$FIX/fq/ids.txt' -m select"
	rc_s=$_RC
	nf1=$(( $( wc -l < "$WORK/fq.flt_R1.fastq" 2>/dev/null || echo 0 ) / 4 ))
	ns1=$(( $( wc -l < "$WORK/fq.sel_R1.fastq" 2>/dev/null || echo 0 ) / 4 ))
	ns2=$(( $( wc -l < "$WORK/fq.sel_R2.fastq" 2>/dev/null || echo 0 ) / 4 ))
	u31a_bad=""
	[ "$rc_f" -eq 0 ] || u31a_bad="$u31a_bad\n  filter mode exited $rc_f"
	[ "$rc_s" -eq 0 ] || u31a_bad="$u31a_bad\n  select mode exited $rc_s"
	[ "$nf1" -eq 6 ]  || u31a_bad="$u31a_bad\n  filter kept $nf1 forward reads, expected 6"
	[ "$ns1" -eq 4 ]  || u31a_bad="$u31a_bad\n  select kept $ns1 forward reads, expected 4"
	[ "$ns2" -eq 4 ]  || u31a_bad="$u31a_bad\n  select kept $ns2 reverse reads, expected 4 (pairs must stay in step)"
	if [ -z "$u31a_bad" ]; then
		pass UNIT-31a "filtfq: filter keeps 6/10 pairs, select keeps 4/10, mates stay paired"
	else
		fail UNIT-31a "filtfq: filter keeps 6/10 pairs, select keeps 4/10, mates stay paired" \
			"$( printf '%b' "$u31a_bad" )"
	fi
fi

# ============================================ UNIT-32 retrieve_reads modes ===

# Fixture: c1 -> r1,r2,r3   c2 -> r4   c3 -> r5,r6, with c1/c2 both assigned to
# taxid 1239574 (Mamastrovirus 10) and c3 to 3048202.  The read files carry a
# seventh read that no contig claims, so a mode that ignores the map entirely
# would return 7 and be caught.
RR="$INSTALL/bin/retrieve_reads"
[ -x "$WORK/retrieve_reads" ] && RR="$WORK/retrieve_reads"
if [ ! -x "$RR" ]; then
	skipt UNIT-32 "retrieve_reads -c / -s / -t" "bin/retrieve_reads not built"
else
	RES="$WORK/res32"
	mkdir -p "$RES/reads"
	cp "$FIX/retrieve/annot_table.tsv" "$FIX/retrieve/readid_contigid.tsv" "$RES/"

	# mode, argument, expected read count
	u32_bad=""
	while IFS='|' read -r mode arg expect label; do
		[ -n "$mode" ] || continue
		out_prefix="probe_${label}"
		try_sh "'$RR' $mode '$arg' -r '$RES' -1 '$FIX/retrieve/read1.fq' \
			-2 '$FIX/retrieve/read2.fq' -p '$out_prefix'"
		rc=$_RC
		f1="$RES/reads/${out_prefix}_r1.fq"
		got=0
		[ -s "$f1" ] && got=$(( $( wc -l < "$f1" ) / 4 ))
		[ "$rc" -eq 0 ] || u32_bad="$u32_bad\n  $mode $arg: exited $rc"
		[ "$got" -eq "$expect" ] \
			|| u32_bad="$u32_bad\n  $mode $arg: retrieved $got reads, expected $expect"
	done <<-'CASES'
		-c|c1|3|contig
		-t|1239574|4|taxid
		-t|3048202|2|taxid2
		-s|Mamastrovirus 10|4|species
	CASES

	# Retrieved ids must be a subset of the input ids, never invented.
	if [ -s "$RES/reads/probe_taxid_r1.fq" ]; then
		extra=$( grep '^@r' "$RES/reads/probe_taxid_r1.fq" | sed 's/^@//; s#/.*##' | sort -u \
			| comm -23 - <( cut -f1 "$RES/readid_contigid.tsv" | sort -u ) | tr '\n' ',' )
		[ -z "$extra" ] || u32_bad="$u32_bad\n  -t 1239574 returned ids absent from the map: $extra"
	fi

	if [ -z "$u32_bad" ]; then
		pass UNIT-32 "retrieve_reads returns the expected reads for -c, -t and -s"
	else
		fail UNIT-32 "retrieve_reads returns the expected reads for -c, -t and -s" \
			"$( printf '%b' "$u32_bad" )"
	fi
fi

# ------------------------------------------------------------------ done ---

# Nothing in this tier may write into the repository working tree.  Compared
# against the snapshot taken before the first test ran, so pre-existing local
# edits are ignored and only changes this run caused are reported.
if [ -z "$TREE_BEFORE" ]; then
	skipt UNIT-99 "the tier wrote nothing into the repository working tree" "not a git checkout"
else
	TREE_AFTER=$( git -C "$REPO" status --porcelain 2>/dev/null )
	if [ "$TREE_BEFORE" = "$TREE_AFTER" ]; then
		pass UNIT-99 "the tier wrote nothing into the repository working tree"
	else
		fail UNIT-99 "the tier wrote nothing into the repository working tree" \
			"$( diff <( printf '%s\n' "$TREE_BEFORE" ) <( printf '%s\n' "$TREE_AFTER" ) | head -10 )" \
			"tests must write only under \$TMPDIR (docs/testing_roadmap.md §1)"
	fi
fi

diag ""
tap_done
