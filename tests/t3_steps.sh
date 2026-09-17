#!/usr/bin/env bash
#
# LazypipeX Tier 3 — step-wise pipeline tests on sample data.
# Implements STEP-01 … STEP-16 of docs/testing_roadmap.md §6.
#
# Every pipeline step is invoked separately against the bundled toy library
# data/samples/M15small_R{1,2}.fastq (9 842 read pairs from a mink faecal
# sample), so a failure localises to one step instead of one long run.
#
# The steps are a chain: each consumes the previous step's output.  A step whose
# prerequisite did not pass is reported as a SKIP rather than run, because a
# cascade of failures from one broken step tells you nothing you did not already
# know from the first one.
#
# Needs installed databases (Tier 2 green) and the full tool chain.  Wall time is
# 10-30 min depending on the databases chosen below.
#
# Usage:
#     module use /projappl/project_2003755/Lazypipe-db/modulefiles/projects
#     module load lazypipe/3.1
#     tests/t3_steps.sh                # TAP on stdout
#     tests/t3_steps.sh | grep -v ^#   # results only
#
# Which databases the annotation steps use is site-specific, so they come from
# the environment.  The defaults are the small virus-only sets: they exercise the
# same code paths as the large ones and keep the tier inside its time budget —
# §6 names minimap.refseq.abv for STEP-06, but that is a 28 GB FASTA and
# minimap2 indexes its target on the fly, which alone would exceed the budget.
#
#     T3_ANN1=minimap.refseq.abv tests/t3_steps.sh
#
# Set T3_KEEP=1 to leave the results tree behind for inspection.
#
# The ICTV step (--pipe ictv, STEP-09b/09c) is off by default: it is supported
# but still in development and is not part of the User Guide, so it must not
# decide the exit status of a routine post-installation run.  Turn it on with
#
#     T3_ICTV=1 tests/t3_steps.sh
#
# It needs blastn.ictv, blastp.ictv and the ICTV.VMR table installed; when any of
# them is absent the two tests skip rather than fail.
#
# Exit status = number of failed tests (0 = all good).  Neither skips nor TODOs
# are failures: skips mark absent optional tooling or an unrun prerequisite,
# TODOs mark known defects documented in docs/testing_roadmap.md §12.

set -uo pipefail

TESTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# shellcheck source=lib/assert.sh
. "$TESTS_DIR/lib/assert.sh"

REPO=$( cd "$TESTS_DIR/.." && pwd )

if [ -n "${LAZYPIPE_INSTALL_DIR:-}" ] && [ -f "$LAZYPIPE_INSTALL_DIR/lazypipe.pl" ]; then
	INSTALL="$LAZYPIPE_INSTALL_DIR"
else
	INSTALL="$REPO"
fi

LZP="$INSTALL/lazypipe.pl"
R1="$INSTALL/data/samples/M15small_R1.fastq"

: "${T3_SAMPLE:=M15test}"
: "${T3_HOSTGEN:=Neovison_vison}"
: "${T3_HOSTGEN_OTHER:=Ixodes_scapularis}"	# a filter that should NOT match mink
: "${T3_ANN1:=minimap.refseq.vi}"		# virus-only and small: seconds, not minutes
# Second engine, for --append and chaining.  The virus-only set is 0.55 GB and
# answers in seconds; uniref100.abv is 33 GB and takes about nine minutes, which
# would be the bulk of this tier's runtime for no extra coverage — .vi already
# grows annot1.tsv and puts a second engine in the 'search' column, which is what
# STEP-07 and STEP-08 assert.  Point this at .abv for a deeper pass.
: "${T3_ANN1_ALT:=diamondp.uniref100.vi}"
# Round 2 against the virus-only index too.  STEP-09 asserts the mechanics — that
# round 2 re-searches a subset of round 1 — and the viral index proves that in a
# second rather than the 95 s blastn.refseq.abv took, which was 42 % of the whole
# tier.  The classical minimap.vi -> blastn.abv two-round strategy is a biological
# question, not a mechanical one, and is covered by Tier 4.
: "${T3_ANN2:=vi:blastn.refseq.vi}"
: "${T3_NUMTH:=8}"
: "${T3_STEP_TIMEOUT:=1800}"
: "${T3_KEEP:=0}"

# Results must never land in the repository.
WORK=$( mktemp -d "${TMPDIR:-/tmp}/lazytest-t3.XXXXXX" ) || exit 99
if [ "$T3_KEEP" = "1" ]; then
	trap 'printf "# results kept at %s\n" "$WORK"' EXIT
else
	trap 'rm -rf "$WORK"' EXIT
fi

RES="$WORK/res"
TMPD="$WORK/tmp"
LOGS="$WORK/logs"
OUT="$RES/$T3_SAMPLE"		# lazypipe.pl appends the sample name to --res
mkdir -p "$RES" "$TMPD" "$LOGS"

cd "$INSTALL" || exit 99

# A tool counts as available only if it is also executable.  `command -v`
# alone is not enough: bash returns the path of a non-executable file found
# on PATH, so a downloaded-but-not-chmod+x binary was reported as installed
# and then failed at run time with "Permission denied".
have() {
	local p
	p=$( type -P "$1" 2>/dev/null ) || return 1
	[ -n "$p" ] && [ -x "$p" ]
}

TREE_BEFORE=$( git -C "$REPO" status --porcelain 2>/dev/null )

# Records of what passed, so a dependent step can skip instead of cascading.
PASSED=""
mark_ok() { PASSED="$PASSED $1"; }

# need ID "desc" PREREQ...  -> 1 (and emits the SKIP) if any prerequisite is missing
need() {
	local id="$1" desc="$2" p
	shift 2
	for p in "$@"; do
		case " $PASSED " in
			*" $p "*) ;;
			*) skipt "$id" "$desc" "prerequisite $p did not pass"; return 1 ;;
		esac
	done
	return 0
}

# Undefined-variable warnings seen across the whole tier, reported once at the
# end.  $TM is expected: config.yaml's par_trimm names it, and it is unset unless
# trimmomatic is installed — harmless while --pre is fastp, which is the default.
UNDEF_VARS=""
STEP_TIMES=""
TOTAL_SECS=0
_SECS=""
N_LZ=0		# successful lazypipe.pl invocations; each must append one provenance.txt section

# Run one pipeline step.  Sets $_RC/$_OUT (via try_sh) and $_XCUT to any
# cross-cutting problem: §6 requires exit 0, no ERROR: on stderr, no undefined
# environment variable, and a new History.log line for every invocation.
lz_step() {
	local args="$1" hist_before=0 hist_after=0 v t0
	[ -s "$OUT/History.log" ] && hist_before=$( wc -l < "$OUT/History.log" )

	t0=$( date +%s )
	TEST_TIMEOUT="$T3_STEP_TIMEOUT" try_sh "perl '$LZP' -1 '$R1' --res '$RES' \
		-s '$T3_SAMPLE' -t $T3_NUMTH --tmpdir '$TMPD' --logs '$LOGS' -v $args"
	# Recorded by report_step, never asserted on: §1 says flag order-of-magnitude
	# drift, not absolute runtimes, which vary with databases and filesystem.
	_SECS=$(( $( date +%s ) - t0 ))

	_XCUT=""
	# a run that dies (e.g. in option checking) never reaches write_provenance()
	[ "$_RC" -eq 0 ] && N_LZ=$(( N_LZ + 1 ))
	[ "$_RC" -eq 0 ] || _XCUT="$_XCUT\n  exit status $_RC"
	if printf '%s' "$_OUT" | grep -q '^ERROR:'; then
		_XCUT="$_XCUT\n  ERROR: on stderr: $( printf '%s' "$_OUT" | grep -m1 '^ERROR:' | cut -c1-90 )"
	fi
	for v in $( printf '%s' "$_OUT" | sed -n 's/.*undefined environment variable "\([A-Za-z_][A-Za-z0-9_]*\)".*/\1/p' | sort -u ); do
		case " $UNDEF_VARS " in *" $v "*) ;; *) UNDEF_VARS="$UNDEF_VARS $v" ;; esac
		[ "$v" = "TM" ] || _XCUT="$_XCUT\n  undefined environment variable \$$v"
	done
	[ -s "$OUT/History.log" ] && hist_after=$( wc -l < "$OUT/History.log" )
	if [ "$hist_after" -le "$hist_before" ]; then
		_XCUT="$_XCUT\n  no new line in History.log"
	fi
	return 0
}

# Reads in a fastq(.gz).
nreads() {
	[ -s "$1" ] || { printf '0'; return 0; }
	case "$1" in
		*.gz) printf '%s' "$(( $( zcat "$1" | wc -l ) / 4 ))" ;;
		*)    printf '%s' "$(( $( wc -l < "$1" ) / 4 ))" ;;
	esac
}

nrecords() { [ -s "$1" ] && grep -c '^>' "$1" || printf '0'; }

# Shortest sequence in a FASTA, for the length-threshold assertions.
min_seqlen() {
	[ -s "$1" ] || { printf '0'; return 0; }
	awk '/^>/ { if(n) print l; l=0; n=1; next } { l+=length($0) } END { if(n) print l }' "$1" \
		| sort -n | head -1
}

# Report a step: $3 empty means the step is good.
report_step() {
	local id="$1" desc="$2" bad="$3"
	# Steps that only inspect earlier output (STEP-05, STEP-11) never call
	# lz_step, so there is no time to attribute to them.
	if [ -n "${_SECS:-}" ]; then
		STEP_TIMES="$STEP_TIMES  $id ${_SECS}s\n"
		TOTAL_SECS=$(( TOTAL_SECS + _SECS ))
		_SECS=""
	fi
	if [ -z "$bad" ]; then
		pass "$id" "$desc"
		mark_ok "$id"
	else
		fail "$id" "$desc" "$( printf '%b' "$bad" )"
	fi
}

# ----------------------------------------------------------------- report ---

tap_init "LazypipeX Tier 3 — step-wise pipeline tests"
diag "host        : $( hostname )"
diag "date        : $( date -Is )"
diag "install dir : $INSTALL"
diag "results     : $RES"
diag "sample      : $T3_SAMPLE  ($( nreads "$R1" ) read pairs)"
diag "hostgen     : $T3_HOSTGEN"
diag "ann1        : $T3_ANN1   (alt: $T3_ANN1_ALT)"
diag "ann2        : $T3_ANN2"
diag ""

if [ ! -s "$R1" ]; then
	fail STEP-00 "the bundled sample library is present" "missing $R1"
	diag ""
	tap_done
	exit $?
fi

# ================================================== STEP-01 preprocess =======

lz_step "-p pre"
s01="$_XCUT"
IN_PAIRS=$( nreads "$R1" )
TRIM1="$OUT/reads/read1.trim.fq.gz"
TRIM2="$OUT/reads/read2.trim.fq.gz"
for f in "$TRIM1" "$TRIM2" "$OUT/reports/fastp.report.html" "$OUT/reports/fastp.json"; do
	[ -s "$f" ] || s01="$s01\n  missing or empty ${f#$OUT/}"
done
TRIM_PAIRS=$( nreads "$TRIM1" )
if [ "$TRIM_PAIRS" -gt 0 ] && [ "$IN_PAIRS" -gt 0 ]; then
	pct=$(( TRIM_PAIRS * 100 / IN_PAIRS ))
	[ "$pct" -ge 80 ] || s01="$s01\n  kept $TRIM_PAIRS/$IN_PAIRS pairs (${pct}%), below the 80% floor"
else
	s01="$s01\n  no trimmed pairs"
fi
# The JSON must parse: a truncated fastp report is a silent corruption.
if [ -s "$OUT/reports/fastp.json" ] && have perl; then
	perl -MJSON::PP -e 'JSON::PP->new->decode(do{local $/; open my $f,"<",$ARGV[0] or die; <$f>})' \
		"$OUT/reports/fastp.json" >/dev/null 2>&1 \
		|| perl -e 'my $s=do{local $/; open my $f,"<",$ARGV[0] or die; <$f>}; die unless $s=~/^\s*\{.*\}\s*$/s' \
			"$OUT/reports/fastp.json" >/dev/null 2>&1 \
		|| s01="$s01\n  reports/fastp.json does not parse as JSON"
fi
report_step STEP-01 "preprocess: trimmed pairs $TRIM_PAIRS/$IN_PAIRS, fastp report written" "$s01"

# ================================================ STEP-02 host filtering =====

if need STEP-02 "host filtering removes host reads" STEP-01; then
	lz_step "-p flt --hostgen $T3_HOSTGEN"
	s02="$_XCUT"
	HFLT1="$OUT/reads/read1.trim.hflt.fq.gz"
	HFLT2="$OUT/reads/read2.trim.hflt.fq.gz"
	for f in "$HFLT1" "$HFLT2"; do
		[ -s "$f" ] || s02="$s02\n  missing or empty ${f#$OUT/}"
	done
	HFLT_PAIRS=$( nreads "$HFLT1" )
	if [ "$HFLT_PAIRS" -le 0 ]; then
		s02="$s02\n  host filtering left 0 reads"
	elif [ "$HFLT_PAIRS" -ge "$TRIM_PAIRS" ]; then
		s02="$s02\n  filtered ($HFLT_PAIRS) is not less than input ($TRIM_PAIRS) — did the filter run?"
	fi
	report_step STEP-02 "host filtering: $HFLT_PAIRS/$TRIM_PAIRS pairs kept against $T3_HOSTGEN" "$s02"
fi

# ==================================================== STEP-03 assembly =======

if need STEP-03 "assembly with megahit" STEP-02; then
	lz_step "-p ass --ass megahit"
	s03="$_XCUT"
	CONTIGS="$OUT/contigs.fa"
	NCONTIG=$( nrecords "$CONTIGS" )
	if [ "$NCONTIG" -lt 1 ]; then
		s03="$s03\n  contigs.fa missing or has no records"
	else
		minlen=$( min_seqlen "$CONTIGS" )
		minreq=$( perl -MYAML::Tiny -e '
			my $y = YAML::Tiny->read("config.yaml");
			print $y->[0]{"general.parameters"}{min_contig_length} // 300;' 2>/dev/null )
		minreq=${minreq:-300}
		[ "$minlen" -ge "$minreq" ] \
			|| s03="$s03\n  shortest contig is ${minlen}nt, below min_contig_length ($minreq)"
	fi
	report_step STEP-03 "assembly: $NCONTIG contigs, shortest $( min_seqlen "$OUT/contigs.fa" )nt" "$s03"
fi

# =============================================== STEP-03b spades assembly ====

if ! have spades.py; then
	skipt STEP-03b "assembly with spades" "spades.py not on PATH"
elif [ "${T3_SPADES:-0}" != "1" ]; then
	# Off by default: a second full assembly doubles the tier's runtime and
	# overwrites contigs.fa, which every later step depends on.
	skipt STEP-03b "assembly with spades" "set T3_SPADES=1 to run (overwrites contigs.fa)"
else
	lz_step "-p ass --ass spades"
	s03b="$_XCUT"
	[ -s "$OUT/contigs.fa" ] || s03b="$s03b\n  missing contigs.fa"
	report_step STEP-03b "assembly with spades" "$s03b"
fi

# ===================================================== STEP-04 realign =======

if need STEP-04 "realign reads to contigs" STEP-03; then
	lz_step "-p rea"
	s04="$_XCUT"
	MAP="$OUT/readid_contigid.tsv"
	if [ ! -s "$MAP" ]; then
		s04="$s04\n  missing or empty readid_contigid.tsv"
	else
		nmap=$( wc -l < "$MAP" )
		[ "$nmap" -gt 0 ] || s04="$s04\n  readid_contigid.tsv has no rows"
		# Every contig named in the map must exist in contigs.fa, or the map is
		# stale relative to the assembly it claims to describe.
		# readid_contigid.tsv is headerless: readid<TAB>contigid, one pair per line.
		orphan=$( awk -F'\t' '{ print $2 }' "$MAP" | sort -u > "$WORK/map_ctg.txt"
			grep '^>' "$OUT/contigs.fa" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$WORK/fa_ctg.txt"
			comm -23 "$WORK/map_ctg.txt" "$WORK/fa_ctg.txt" | head -3 | tr '\n' ' ' )
		[ -z "$orphan" ] || s04="$s04\n  contig ids in the map are absent from contigs.fa: $orphan"
	fi
	report_step STEP-04 "realign: $( [ -s "$MAP" ] && wc -l < "$MAP" || echo 0 ) read-to-contig rows" "$s04"
fi

# ============================================== STEP-05 ORF prediction =======

# ORFs are produced by the realign step, not by a --pipe step of their own, so
# this asserts on what STEP-04 already wrote rather than running the pipeline
# again.  Only the default predictor (--gen mga) is exercised: prodigal is being
# discontinued in favour of orfipy, so pinning a test to --gen prod would pin it
# to a code path on its way out.
if need STEP-05 "ORF prediction" STEP-04; then
	s05=""
	AA="$OUT/contigs.orfs.aa.fa"
	NT="$OUT/contigs.orfs.nt.fa"
	for f in "$AA" "$NT"; do
		[ -s "$f" ] || s05="$s05\n  missing or empty ${f#$OUT/}"
	done
	if [ -s "$AA" ] && [ -s "$NT" ]; then
		naa=$( nrecords "$AA" )
		nnt=$( nrecords "$NT" )
		[ "$naa" -eq "$nnt" ] || s05="$s05\n  aa has $naa records, nt has $nnt — they must match"
		minorf=$( min_seqlen "$NT" )
		minreq=$( perl -MYAML::Tiny -e '
			my $y = YAML::Tiny->read("config.yaml");
			print $y->[0]{"general.parameters"}{min_orf_length} // 72;' 2>/dev/null )
		minreq=${minreq:-72}
		[ "$minorf" -ge "$minreq" ] \
			|| s05="$s05\n  shortest ORF is ${minorf}nt, below min_orf_length ($minreq)"
	fi
	report_step STEP-05 "ORF prediction: $( nrecords "$AA" ) ORFs, shortest $( min_seqlen "$NT" )nt" "$s05"
fi

# ======================================= STEP-06 annotation round 1 ==========

# The column list in §6 is from an older generation of the file: it names
# contig/clen/species, while annot1.tsv carries qseqid/qseqlen and no species
# column.  Asserted here against what the pipeline actually writes.
ANNOT1_COLS="search db dbtype qseqid orf qseqlen sseqid bitscore alen pident qlen qcov slen scov staxid sname bphage division"

if need STEP-06 "annotation round 1 with $T3_ANN1" STEP-05; then
	lz_step "-p ann1 --ann1 $T3_ANN1"
	s06="$_XCUT"
	A1="$OUT/annot1.tsv"
	if [ ! -s "$A1" ]; then
		s06="$s06\n  missing or empty annot1.tsv"
	else
		hdr=$( head -1 "$A1" )
		for c in $ANNOT1_COLS; do
			printf '%s' "$hdr" | tr '\t' '\n' | grep -qx -- "$c" \
				|| s06="$s06\n  annot1.tsv has no '$c' column"
		done
		rows=$(( $( wc -l < "$A1" ) - 1 ))
		[ "$rows" -ge 1 ] || s06="$s06\n  annot1.tsv has a header but no rows"
		# The sample is a mink faecal library with known viral content, so a run
		# that annotates nothing as viral has found nothing worth finding.
		dcol=$( printf '%s' "$hdr" | tr '\t' '\n' | grep -nx division | cut -d: -f1 )
		if [ -n "$dcol" ]; then
			nvi=$( awk -F'\t' -v c="$dcol" 'NR>1 && $c=="Viruses" { n++ } END { print n+0 }' "$A1" )
			[ "$nvi" -ge 1 ] || s06="$s06\n  no row with division=Viruses; the viral contigs went unannotated"
		fi
		# and the viral contigs must have been written out for round 2
		[ -s "$OUT/contigs.ann1.vi.fa" ] \
			|| s06="$s06\n  contigs.ann1.vi.fa is missing or empty"
		# Every staxid must resolve, or the downstream binning is built on sand.
		if have taxonkit && [ -n "${taxonomy_ncbi:-}" ] && [ "$rows" -ge 1 ]; then
			col=$( printf '%s' "$hdr" | tr '\t' '\n' | grep -nx staxid | cut -d: -f1 )
			awk -F'\t' -v c="$col" 'NR>1 && $c ~ /^[0-9]+$/ { print $c }' "$A1" | sort -u > "$WORK/staxids.txt"
			if [ -s "$WORK/staxids.txt" ]; then
				unres=$( taxonkit lineage --data-dir "$taxonomy_ncbi" < "$WORK/staxids.txt" 2>/dev/null \
					| awk -F'\t' '$2 == "" { print $1 }' | head -3 | tr '\n' ' ' )
				[ -z "$unres" ] || s06="$s06\n  staxids that do not resolve in the taxonomy: $unres"
			fi
		fi
	fi
	report_step STEP-06 "annotation round 1: $( [ -s "$A1" ] && echo $(( $( wc -l < "$A1" ) - 1 )) || echo 0 ) rows from $T3_ANN1" "$s06"
fi

# ======================================= STEP-07 round 1, --append ===========

if need STEP-07 "round 1 --append grows annot1.tsv" STEP-06; then
	rows_before=$(( $( wc -l < "$OUT/annot1.tsv" ) - 1 ))
	lz_step "-p ann1 --ann1 $T3_ANN1_ALT --append"
	s07="$_XCUT"
	rows_after=$(( $( wc -l < "$OUT/annot1.tsv" ) - 1 ))
	if [ "$rows_after" -le "$rows_before" ]; then
		s07="$s07\n  row count did not increase: $rows_before -> $rows_after"
	fi
	# Both engines must be represented, or --append silently replaced instead.
	nsearch=$( awk -F'\t' 'NR>1 { print $1 }' "$OUT/annot1.tsv" | sort -u | grep -c . )
	[ "$nsearch" -ge 2 ] || s07="$s07\n  only $nsearch distinct 'search' value(s) after --append; expected both engines"
	report_step STEP-07 "round 1 --append: $rows_before -> $rows_after rows, $nsearch engines" "$s07"
fi

# ======================================= STEP-08 round 1, chained ============

if need STEP-08 "round 1 chaining is disjoint" STEP-06; then
	lz_step "-p ann1 --ann1 $T3_ANN1,$T3_ANN1_ALT"
	s08="$_XCUT"
	A1="$OUT/annot1.tsv"
	if [ ! -s "$A1" ]; then
		s08="$s08\n  missing annot1.tsv"
	else
		# A chain hands the second engine only what the first did not annotate,
		# so the two contig sets must not overlap.
		awk -F'\t' 'NR>1 { print $1"\t"$4 }' "$A1" | sort -u > "$WORK/chain.tsv"
		eng1=$( awk -F'\t' 'NR>1 { print $1 }' "$A1" | sort -u | head -1 )
		eng2=$( awk -F'\t' 'NR>1 { print $1 }' "$A1" | sort -u | sed -n 2p )
		if [ -z "$eng2" ]; then
			diag "  note: only engine '$eng1' produced hits, so the chain had nothing to hand on"
		else
			awk -F'\t' -v e="$eng1" '$1==e { print $2 }' "$WORK/chain.tsv" | sort -u > "$WORK/c1.txt"
			awk -F'\t' -v e="$eng2" '$1==e { print $2 }' "$WORK/chain.tsv" | sort -u > "$WORK/c2.txt"
			both=$( comm -12 "$WORK/c1.txt" "$WORK/c2.txt" | head -3 | tr '\n' ' ' )
			[ -z "$both" ] || s08="$s08\n  contigs annotated by both $eng1 and $eng2: $both"
		fi
	fi
	report_step STEP-08 "round 1 chaining: engines annotate disjoint contig sets" "$s08"
fi

# ======================================= STEP-09 annotation round 2 ==========

if need STEP-09 "annotation round 2" STEP-05; then
	lz_step "-p ann1,ann2 --ann1 $T3_ANN1 --ann2 $T3_ANN2"
	s09="$_XCUT"
	A2="$OUT/annot2.tsv"
	if [ ! -s "$A2" ]; then
		s09="$s09\n  missing or empty annot2.tsv"
	else
		# Round 2 re-searches a subset selected from round 1, so its contigs must
		# be a subset of the round-1 contigs.
		awk -F'\t' 'NR>1 { print $4 }' "$OUT/annot1.tsv" | sort -u > "$WORK/a1_ctg.txt"
		awk -F'\t' 'NR>1 { print $4 }' "$A2"            | sort -u > "$WORK/a2_ctg.txt"
		extra=$( comm -13 "$WORK/a1_ctg.txt" "$WORK/a2_ctg.txt" | head -3 | tr '\n' ' ' )
		[ -z "$extra" ] || s09="$s09\n  round-2 contigs absent from round 1: $extra"
	fi
	report_step STEP-09 "annotation round 2: $( [ -s "$A2" ] && echo $(( $( wc -l < "$A2" ) - 1 )) || echo 0 ) rows" "$s09"
fi

# ==================================== STEP-09b/09c ICTV annotation + EM ======

# --pipe ictv re-annotates the viral contigs against the ICTV exemplar and
# additional-isolate indexes, then runs EM_loop() (perl/Lazypipe/SeqAn.pm) to
# spread each contig's alignment scores over isolates and report a probability
# per isolate and per genus.  Not covered by docs/testing_roadmap.md §6: the step
# is supported but in development and undocumented in the User Guide, so it is
# opt-in and never contributes a failure unless it was asked for.
#
# Split in two so that a toy library with no ICTV hit skips the probability
# assertions instead of passing them vacuously: 09b is the plumbing, 09c is what
# the EM produced.
ICTV_VMR=$( perl -MYAML::Tiny -e '
	my $y = YAML::Tiny->read("config.yaml");
	my $p = $y->[0]{"ICTV.VMR"}{db} // "";
	$p =~ s{\$(\w+)}{ defined($ENV{$1}) ? $ENV{$1} : "\$$1" }ge;
	print $p;' 2>/dev/null )

# Installed means what it means everywhere else in the pipeline: glob("$db*")
# is non-empty, which is exactly what --databases reports.
ICTV_DBS=""
if [ "${T3_ICTV:-0}" = "1" ]; then
	try perl "$LZP" --databases
	for d in blastn.ictv blastp.ictv; do
		printf '%s\n' "$_OUT" | grep -qx -- "$d:" || ICTV_DBS="$ICTV_DBS $d"
	done
fi

ICTV_OUT=""
ICTV_D="$OUT/ictv"
if [ "${T3_ICTV:-0}" != "1" ]; then
	skipt STEP-09b "ICTV annotation writes ictv/ictv.annot.tsv" \
		"set T3_ICTV=1 to run (step is in development and not in the User Guide)"
elif [ -n "$ICTV_DBS" ]; then
	skipt STEP-09b "ICTV annotation writes ictv/ictv.annot.tsv" \
		"not installed:$ICTV_DBS (perl/install_db.pl --db blastn.ictv --db blastp.ictv)"
elif [ ! -s "$ICTV_VMR" ]; then
	skipt STEP-09b "ICTV annotation writes ictv/ictv.annot.tsv" \
		"no ICTV VMR table at ${ICTV_VMR:-<ICTV.VMR:db unset in config.yaml>}"
elif need STEP-09b "ICTV annotation writes ictv/ictv.annot.tsv" STEP-06; then
	lz_step "-p ictv"
	s09b="$_XCUT"
	ICTV_OUT="$_OUT"		# EM_loop() traces its iterations to stderr under -v
	AI="$ICTV_D/ictv.annot.tsv"
	if [ ! -s "$AI" ]; then
		s09b="$s09b\n  missing or empty ictv/ictv.annot.tsv"
	else
		# The three columns EM_q_t_logprob() requires; it dies without them.
		hdr=$( head -1 "$AI" )
		for c in qseqid staxid bitscore; do
			printf '%s' "$hdr" | tr '\t' '\n' | grep -qx -- "$c" \
				|| s09b="$s09b\n  ictv.annot.tsv has no '$c' column, which the EM requires"
		done
		# Subject taxids here are ICTV Isolate.NIDs, and lazypipe.pl derives them
		# with ^VMR([0-9]+); anything else joins against nothing in the VMR table.
		scol=$( printf '%s' "$hdr" | tr '\t' '\n' | grep -nx staxid | cut -d: -f1 )
		if [ -n "$scol" ]; then
			bad=$( awk -F'\t' -v c="$scol" 'NR>1 && $c !~ /^[0-9]+$/ { print $c; exit }' "$AI" )
			[ -z "$bad" ] || s09b="$s09b\n  non-numeric Isolate.NID in staxid: $bad"
		fi
		# The EM fits over whatever is in staxid, and every table behind it joins
		# that against ICTV.VMR on Isolate.NID.  A database keyed on anything else
		# — NCBI taxids, say — still annotates, still runs the EM, and then joins
		# to NA, so the reports come out empty with no error anywhere.  Numeric
		# staxids are therefore not enough: they have to be ICTV isolates.
		if [ -n "$scol" ] && [ -s "$ICTV_VMR" ]; then
			awk -F'\t' -v c="$scol" 'NR>1 { print $c }' "$AI" | sort -u > "$WORK/ictv_staxid.txt"
			awk -F'\t' 'NR==1 { for(i=1;i<=NF;i++) if($i=="Isolate ID") c=i; next }
			            c     { sub(/^VMR/,"",$c); print $c }' "$ICTV_VMR" | sort -u > "$WORK/ictv_vmrnid.txt"
			nhit=$( comm -12 "$WORK/ictv_staxid.txt" "$WORK/ictv_vmrnid.txt" | grep -c . )
			nsub=$( grep -c . "$WORK/ictv_staxid.txt" )
			[ "$nhit" -gt 0 ] || s09b="$s09b\n  none of the $nsub staxid(s) is an Isolate.NID in $( basename "$ICTV_VMR" )"
			[ "$nhit" -gt 0 ] || s09b="$s09b\n  the ICTV database is not keyed on ICTV isolates, so every ICTV report joins to NA and is filtered away"
		fi
	fi
	ICTV_ROWS=$( [ -s "$AI" ] && echo $(( $( wc -l < "$AI" ) - 1 )) || echo 0 )
	report_step STEP-09b "ICTV annotation: $ICTV_ROWS rows in ictv/ictv.annot.tsv" "$s09b"
fi

if [ -n "${ICTV_ROWS:-}" ] && [ "${ICTV_ROWS:-0}" -eq 0 ]; then
	skipt STEP-09c "EM assigns probabilities to ICTV isolates" \
		"the ICTV search found no hit in this library, so the EM had nothing to fit"
elif need STEP-09c "EM assigns probabilities to ICTV isolates" STEP-09b; then
	s09c=""
	PROB="$ICTV_D/ictv.taxid_prob.tsv"
	ITER="$ICTV_D/ictv.taxid_prob_byiter.tsv"

	# --- what the EM wrote ---------------------------------------------------
	if [ ! -s "$PROB" ]; then
		s09c="$s09c\n  missing or empty ictv/ictv.taxid_prob.tsv"
	else
		phdr=$( head -1 "$PROB" )
		for c in Isolate.NID Species Genus prob prob.genus; do
			printf '%s' "$phdr" | tr '\t' '\n' | grep -qx -- "$c" \
				|| s09c="$s09c\n  ictv.taxid_prob.tsv has no '$c' column"
		done
		pcol=$( printf '%s' "$phdr" | tr '\t' '\n' | grep -nx prob | cut -d: -f1 )
		ntax=$(( $( wc -l < "$PROB" ) - 1 ))
		[ "$ntax" -ge 1 ] || s09c="$s09c\n  ictv.taxid_prob.tsv has a header but no isolate"
		if [ -n "$pcol" ] && [ "$ntax" -ge 1 ]; then
			# F(t) is a distribution over the isolates the EM was given: every
			# value in [0,1] and the column summing to 1.  A sum that drifts off 1
			# means mass was lost or double counted in the M-step.
			oob=$( awk -F'\t' -v c="$pcol" 'NR>1 && ($c+0 < 0 || $c+0 > 1) { print $c; exit }' "$PROB" )
			[ -z "$oob" ] || s09c="$s09c\n  prob outside [0,1]: $oob"
			psum=$( awk -F'\t' -v c="$pcol" 'NR>1 { s+=$c } END { printf "%.6f", s+0 }' "$PROB" )
			awk -v s="$psum" 'BEGIN { exit !(s > 0.99 && s < 1.01) }' \
				|| s09c="$s09c\n  prob column sums to $psum, not 1: the EM posterior does not normalise"
			diag "  EM fitted $ntax isolate(s), prob sums to $psum"
		fi
	fi

	# --- the per-iteration trace --------------------------------------------
	if [ ! -s "$ITER" ]; then
		s09c="$s09c\n  missing or empty ictv/ictv.taxid_prob_byiter.tsv"
	elif [ -s "$PROB" ]; then
		nit=$(( $( head -1 "$ITER" | tr '\t' '\n' | grep -c . ) - 2 ))	# minus Isolate.NID, Species
		[ "$nit" -ge 2 ] || s09c="$s09c\n  the iteration trace has $nit iteration column(s); EM runs at least the initial one plus three steps"
		[ "$( wc -l < "$ITER" )" -eq "$( wc -l < "$PROB" )" ] \
			|| s09c="$s09c\n  taxid_prob_byiter.tsv and taxid_prob.tsv disagree on the number of isolates"
	fi

	# --- the EM's own trace on stderr ---------------------------------------
	# EM_loop() prints "iteration = N, logL = X, logLdiff = Y" per step under -v.
	nem=$( printf '%s\n' "$ICTV_OUT" | grep -cF 'EM_loop(): iteration' )
	if [ "$nem" -lt 2 ]; then
		s09c="$s09c\n  EM_loop() logged $nem iteration(s); it did not run"
	else
		# EM increases the likelihood at every step by construction.  A negative
		# logLdiff is a defect in the E- or M-step, not a tolerance question.
		drop=$( printf '%s\n' "$ICTV_OUT" \
			| sed -n 's/.*logLdiff = \(-\{0,1\}[0-9.e+-]*\).*/\1/p' \
			| awk '$1+0 < -1e-6 { print $1; exit }' )
		[ -z "$drop" ] || s09c="$s09c\n  log-likelihood decreased during EM (logLdiff = $drop)"
		# Converged, rather than stopped by maxiter with the fit still moving.
		printf '%s\n' "$ICTV_OUT" | grep -qF 'exiting EM' \
			|| s09c="$s09c\n  EM ran $nem iterations without converging (no 'exiting EM'); it hit maxiter"
	fi

	# --- the reports the probabilities feed ----------------------------------
	AT="$ICTV_D/ictv.annot_table.tsv"
	if [ ! -s "$AT" ]; then
		s09c="$s09c\n  missing or empty ictv/ictv.annot_table.tsv"
	else
		athdr=$( head -1 "$AT" )
		for c in Isolate.prob Genus.prob Species Genus Family; do
			printf '%s' "$athdr" | tr '\t' '\n' | grep -qx -- "$c" \
				|| s09c="$s09c\n  ictv.annot_table.tsv has no '$c' column"
		done
		# lazypipe.pl filters this table at Genus.prob >= 0.001; a row below the
		# cutoff means the filter did not run.
		gcol=$( printf '%s' "$athdr" | tr '\t' '\n' | grep -nx 'Genus.prob' | cut -d: -f1 )
		if [ -n "$gcol" ]; then
			low=$( awk -F'\t' -v c="$gcol" 'NR>1 && $c ~ /^[0-9.eE+-]+$/ && $c+0 < 0.001 { print $c; exit }' "$AT" )
			[ -z "$low" ] || s09c="$s09c\n  row kept with Genus.prob = $low, below the 0.001 cutoff"
		fi
	fi
	# Both workbooks come from "$perl_scripts/write_excel.pl" (lazypipe.pl:1546 and
	# :1590), which is not in the repository — the shell redirect still creates the
	# file, so an empty ictv.*.xlsx means the generator was never there rather than
	# that it wrote nothing.  Checked by the zip magic (PK\x03\x04); the missing
	# script also aborts the step, which the cross-cutting exit-status check sees.
	for x in ictv.annot_table.xlsx ictv.abund_table.xlsx; do
		f="$ICTV_D/$x"
		if [ ! -s "$f" ]; then
			s09c="$s09c\n  missing or empty ictv/$x — is perl/write_excel.pl installed?"
		elif [ "$( head -c 2 "$f" )" != "PK" ]; then
			s09c="$s09c\n  ictv/$x is not a workbook (no PK zip magic)"
		fi
	done

	report_step STEP-09c "EM probabilities: taxid_prob, iteration trace and ICTV tables" "$s09c"
fi

# ======================================================= STEP-10 reports =====

if need STEP-10 "reports are written" STEP-09; then
	lz_step "-p rep"
	s10="$_XCUT"
	# taxprofile.txt is deliberately absent from this list — see STEP-10a.
	for f in abund_table.tsv abund_table.xlsx annot_table.tsv annot_table.xlsx \
	         reports/krona.report.html reports/krona.data.txt; do
		[ -s "$OUT/$f" ] || s10="$s10\n  missing or empty $f"
	done
	[ -d "$OUT/contigs" ] || s10="$s10\n  missing contigs/ directory"

	# §6 asks for a readn_pc column summing to ~100 %.  There is no such column:
	# abund_table.tsv carries absolute readn, and readn_pc exists only inside
	# R/NGSlib.R, which derives it for the plots.  Assert what the table has —
	# reads are attributed, and never more than were assembled.
	if [ -s "$OUT/abund_table.tsv" ]; then
		col=$( head -1 "$OUT/abund_table.tsv" | tr '\t' '\n' | grep -nx readn | cut -d: -f1 )
		if [ -z "$col" ]; then
			s10="$s10\n  abund_table.tsv has no readn column"
		else
			sum=$( awk -F'\t' -v c="$col" 'NR>1 { s+=$c } END { printf "%d", s+0 }' "$OUT/abund_table.tsv" )
			if [ "$sum" -le 0 ]; then
				s10="$s10\n  abund_table.tsv attributes 0 reads"
			elif [ "$sum" -gt "$(( TRIM_PAIRS * 2 ))" ]; then
				s10="$s10\n  abund_table.tsv attributes $sum reads, more than the $(( TRIM_PAIRS * 2 )) that entered"
			fi
		fi
	fi

	# The annotation table is the file §6's column list actually describes.
	if [ -s "$OUT/annot_table.tsv" ]; then
		hdr=$( head -1 "$OUT/annot_table.tsv" )
		for c in contig clen staxid sname division species genus family; do
			printf '%s' "$hdr" | tr '\t' '\n' | grep -qx -- "$c" \
				|| s10="$s10\n  annot_table.tsv has no '$c' column"
		done
	fi
	report_step STEP-10 "reports: abundance and annotation tables, krona, contigs/" "$s10"
fi

# =============================================== STEP-10a taxprofile ========

# The User Guide (Table of outputs) documents taxprofile.txt as a CAMI-format
# profile, but the only line that would write it — the abundtable2taxprofile.pl
# call at lazypipe.pl:1089 — is commented out, so the file is never produced.
if need STEP-10a "taxprofile.txt is written" STEP-10; then
	if [ -s "$OUT/taxprofile.txt" ]; then
		pass STEP-10a "taxprofile.txt is written"
		mark_ok STEP-10a
	else
		todof STEP-10a "taxprofile.txt is written" \
			"known defect, docs/testing_roadmap.md §12 item 16" \
			"the generator at lazypipe.pl:1089 is commented out, so the file is never produced" \
			"the User Guide documents it as a pipeline output"
	fi
fi

# ======================================= STEP-11 contig sorting =============

if need STEP-11 "sorted contig sets partition contigs.fa" STEP-10; then
	s11=""
	: > "$WORK/sorted_all.txt"
	present=0
	for g in ab ph vi un; do
		f="$OUT/contigs.ann1.$g.fa"
		[ -e "$f" ] || continue
		present=$(( present + 1 ))
		grep '^>' "$f" 2>/dev/null | sed 's/^>//; s/[[:space:]].*//' >> "$WORK/sorted_all.txt"
	done
	if [ "$present" -eq 0 ]; then
		s11="$s11\n  none of contigs.ann1.{ab,ph,vi,un}.fa exist"
	else
		sort "$WORK/sorted_all.txt" > "$WORK/sorted_sorted.txt"
		sort -u "$WORK/sorted_all.txt" > "$WORK/sorted_uniq.txt"
		dup=$( comm -23 "$WORK/sorted_sorted.txt" "$WORK/sorted_uniq.txt" | head -3 | tr '\n' ' ' )
		[ -z "$dup" ] || s11="$s11\n  contigs appearing in more than one set: $dup"
		grep '^>' "$OUT/contigs.fa" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$WORK/fa_all.txt"
		lost=$( comm -23 "$WORK/fa_all.txt" "$WORK/sorted_uniq.txt" | head -3 | tr '\n' ' ' )
		[ -z "$lost" ] || s11="$s11\n  contigs in contigs.fa missing from every set: $lost"
		alien=$( comm -13 "$WORK/fa_all.txt" "$WORK/sorted_uniq.txt" | head -3 | tr '\n' ' ' )
		[ -z "$alien" ] || s11="$s11\n  contigs in a set but not in contigs.fa: $alien"
	fi
	report_step STEP-11 "contig sorting: $present sets partition contigs.fa with no overlap or loss" "$s11"
fi

# ====================================== STEP-12 reference-genome reports =====

# The gate is the NCBI datasets CLI, not igv-reports: generate_igv_html() writes
# the IGV pages itself with filebin2uri(), while the step fetches reference
# genomes with `datasets` / `dataformat`.
if ! have datasets || ! have dataformat; then
	skipt STEP-12 "reference-genome reports" "NCBI datasets CLI (datasets/dataformat) not on PATH"
elif need STEP-12 "reference-genome reports" STEP-10; then
	lz_step "-p rgrep"
	s12="$_XCUT"
	RG="$OUT/reports/refgen.report.html"
	if [ ! -s "$RG" ]; then
		s12="$s12\n  missing or empty reports/refgen.report.html"
	else
		# The summary page is only a table of links and carries no data: URI —
		# use_data_uri is set for the per-species pages (lazypipe.pl:1271), not
		# for the summary.  Asserting on the summary alone also misses the real
		# failure mode: when the NCBI download fails the step warns, still writes
		# the summary, and leaves every link behind it dead.  So check the linked
		# pages exist, and that the tracks in them are embedded.
		nigv=$( find "$OUT/reports" -name '*.igv.html' 2>/dev/null | wc -l )
		if [ "$nigv" -lt 1 ]; then
			s12="$s12\n  summary written but no per-species *.igv.html behind its links"
		else
			igv=$( find "$OUT/reports" -name '*.igv.html' 2>/dev/null | head -1 )
			grep -q 'data:' "$igv" \
				|| s12="$s12\n  $( basename "$igv" ) embeds no data: URI"
		fi
	fi
	report_step STEP-12 "reference-genome reports ($( find "$OUT/reports" -name '*.igv.html' 2>/dev/null | wc -l ) IGV pages)" "$s12"
fi

# ================================================ STEP-13 stats + QC =========

if need STEP-13 "stats and QC plots" STEP-10; then
	lz_step "-p sta"
	s13="$_XCUT"
	npng=0
	# The pipeline writes plots under figures/ and reports/figures/ depending on
	# the plot; take either.
	for p in "$OUT"/figures/*.png "$OUT"/reports/figures/*.png; do
		[ -e "$p" ] || continue
		npng=$(( npng + 1 ))
		[ -s "$p" ] || { s13="$s13\n  empty PNG: $( basename "$p" )"; continue; }
		# PNG magic: \x89PNG
		magic=$( head -c 4 "$p" | od -An -tx1 | tr -d ' \n' )
		[ "$magic" = "89504e47" ] || s13="$s13\n  $( basename "$p" ) is not a PNG (magic $magic)"
	done
	[ "$npng" -gt 0 ] || s13="$s13\n  no PNG written to figures/"
	report_step STEP-13 "stats and QC plots: $npng PNG(s)" "$s13"
fi

# ========================================================= STEP-14 pack ======

if need STEP-14 "pack writes a tarball" STEP-10; then
	lz_step "-p pack"
	s14="$_XCUT"
	TARB=$( ls -1 "$RES/$T3_SAMPLE".tar.gz "$OUT".tar.gz "$RES"/*.tar.gz 2>/dev/null | head -1 )
	if [ -z "$TARB" ] || [ ! -s "$TARB" ]; then
		s14="$s14\n  no tarball produced under $RES"
	else
		listing=$( tar -tzf "$TARB" 2>/dev/null )
		printf '%s' "$listing" | grep -q 'reports/' || s14="$s14\n  tarball lists no reports/ entry"
		printf '%s' "$listing" | grep -q 'abund_table' || s14="$s14\n  tarball lists no abundance table"
		printf '%s' "$listing" | grep -q 'provenance.txt' || s14="$s14\n  tarball lists no provenance.txt"
	fi
	report_step STEP-14 "pack: tarball lists reports/ and the abundance tables" "$s14"
fi

# ======================================================== STEP-15 clean ======

if need STEP-15 "clean removes intermediates and keeps reports" STEP-14; then
	lz_step "-p clean"
	s15="$_XCUT"
	# The reports must survive the clean.
	for f in abund_table.tsv annot_table.tsv; do
		[ -s "$OUT/$f" ] || s15="$s15\n  clean removed $f, which must survive"
	done
	# Intermediates that clean is expected to take away.
	for f in contigs.bwa.sam contigs.fa.bwt; do
		[ -e "$OUT/$f" ] && s15="$s15\n  intermediate still present after clean: $f"
	done
	report_step STEP-15 "clean: intermediates gone, reports intact" "$s15"
fi

# =================================================== STEP-15b provenance =====

# write_provenance() appends one section per lazypipe.pl invocation, so after the
# steps above provenance.txt must hold one section per successful lz_step call,
# and the sections of the steps that used a database or tool must name it.
if need STEP-15b "provenance.txt has one section per run" STEP-01; then
	s15b=""
	PROV="$OUT/provenance.txt"
	if [ ! -s "$PROV" ]; then
		s15b="$s15b\n  missing or empty provenance.txt"
	else
		nsec=$( grep -c '^# Provenance of a LazypipeX run$' "$PROV" )
		[ "$nsec" -eq "$N_LZ" ] || s15b="$s15b\n  $nsec sections for $N_LZ successful lazypipe.pl invocations"
		for h in '^pipeline:$' '^run:$' '^databases:$' '^tools (versions at run time):$' '^R packages:$'; do
			n=$( grep -c -- "$h" "$PROV" )
			[ "$n" -eq "$nsec" ] || s15b="$s15b\n  '$h' appears $n times in $nsec sections"
		done
		grep -q "^sample: *$T3_SAMPLE$" "$PROV" || s15b="$s15b\n  no 'sample: $T3_SAMPLE' line"
		grep -q '^  lazypipe.pl:  sha256 [0-9a-f]\{64\}$' "$PROV" || s15b="$s15b\n  no sha256 of lazypipe.pl"
		# the databases this run used are listed: round 1 and the host genome
		case " $PASSED " in
			*" STEP-06 "*)
				grep -q "^  $T3_ANN1 " "$PROV" || s15b="$s15b\n  round-1 database $T3_ANN1 not listed" ;;
		esac
		grep -q "^  $T3_HOSTGEN " "$PROV" || s15b="$s15b\n  host genome $T3_HOSTGEN not listed"
		if have megahit && ! grep -q '^  megahit  .*[0-9]' "$PROV"; then
			s15b="$s15b\n  megahit is on PATH but its version is not recorded"
		fi
	fi
	report_step STEP-15b "provenance: $( [ -s "$PROV" ] && grep -c '^# Provenance of a LazypipeX run$' "$PROV" || echo 0 ) sections for $N_LZ runs" "$s15b"
fi

# =============================================== STEP-16 read retrieval ======

RR="$INSTALL/bin/retrieve_reads"
if [ ! -x "$RR" ]; then
	skipt STEP-16 "read retrieval with retrieve_reads" "bin/retrieve_reads not built"
elif need STEP-16 "read retrieval with retrieve_reads" STEP-04; then
	s16=""
	# retrieve_reads reads only uncompressed fastq (§12 item 15), so decompress
	# first exactly as the User Guide instructs.
	gunzip -kf "$OUT"/reads/read1.trim.fq.gz "$OUT"/reads/read2.trim.fq.gz 2>/dev/null
	TOPC=$( awk -F'\t' '{ print $2 }' "$OUT/readid_contigid.tsv" | sort | uniq -c | sort -rn | head -1 | awk '{print $2}' )
	if [ -z "$TOPC" ]; then
		s16="$s16\n  could not pick a contig from readid_contigid.tsv"
	else
		TEST_TIMEOUT="$T3_STEP_TIMEOUT" try_sh "'$RR' -r '$OUT' -c '$TOPC' -p t3probe"
		[ "$_RC" -eq 0 ] || s16="$s16\n  retrieve_reads -c $TOPC exited $_RC"
		got="$OUT/reads/t3probe_r1.fq"
		if [ ! -s "$got" ]; then
			s16="$s16\n  no reads written for contig $TOPC"
		else
			nret=$( nreads "$got" )
			nexp=$( awk -F'\t' -v c="$TOPC" '$2==c { n++ } END { print n+0 }' "$OUT/readid_contigid.tsv" )
			[ "$nret" -gt 0 ] || s16="$s16\n  retrieved 0 reads for $TOPC"
			# Retrieved ids must be a subset of the library.
			sed -n '1~4p' "$got" | sed 's/^@//; s#[/[:space:]].*##' | sort -u > "$WORK/got_ids.txt"
			sed -n '1~4p' "$OUT/reads/read1.trim.fq" | sed 's/^@//; s#[/[:space:]].*##' | sort -u > "$WORK/lib_ids.txt"
			alien=$( comm -13 "$WORK/lib_ids.txt" "$WORK/got_ids.txt" | head -3 | tr '\n' ' ' )
			[ -z "$alien" ] || s16="$s16\n  retrieved read ids absent from the library: $alien"
			diag "  retrieved $nret reads for contig $TOPC (map lists $nexp)"
		fi
	fi
	report_step STEP-16 "read retrieval: reads for the top contig are a subset of the library" "$s16"
fi

# ------------------------------------------------------------------ done ---

if [ -n "$STEP_TIMES" ]; then
	diag "step wall time (recorded, never asserted on):"
	diag "$( printf '%b' "$STEP_TIMES" )"
	diag "  total in pipeline steps: ${TOTAL_SECS}s"
fi

if [ -n "$UNDEF_VARS" ]; then
	diag "undefined environment variables seen during the run:$UNDEF_VARS"
	diag "  \$TM is expected while --pre is fastp: config.yaml names it only in par_trimm"
fi

if [ -z "$TREE_BEFORE" ]; then
	skipt STEP-99 "the tier wrote nothing into the repository working tree" "not a git checkout"
else
	TREE_AFTER=$( git -C "$REPO" status --porcelain 2>/dev/null )
	if [ "$TREE_BEFORE" = "$TREE_AFTER" ]; then
		pass STEP-99 "the tier wrote nothing into the repository working tree"
	else
		fail STEP-99 "the tier wrote nothing into the repository working tree" \
			"$( diff <( printf '%s\n' "$TREE_BEFORE" ) <( printf '%s\n' "$TREE_AFTER" ) | head -10 )" \
			"tests must write only under \$TMPDIR (docs/testing_roadmap.md §1)"
	fi
fi

diag ""
tap_done
