#!/usr/bin/env bash
#
# LazypipeX Tier 2 — database installation tests.
# Implements DB-01 … DB-06 of docs/testing_roadmap.md §5.1 and DB-10 … DB-27 of §5.2.
#
# Needs no network and downloads nothing.  §5.1 reads config.yaml, asks
# lazypipe.pl what it considers installed and lints the two database sections;
# §5.2 opens each installed database with its own engine and queries a small
# virus-only one per engine.
#
# install_db.pl is run twice and both calls are safe: once with a key asserted
# absent from config.yaml (DB-06), and once with a key asserted already
# installed (DB-27), where install_db() returns before reaching its wget.
#
# DB-06 covers only the part of "download integrity" that is checkable without
# a network — install_db.pl's failure contract.  The per-engine half of it is
# DB-10 … DB-27 below.
#
# Usage:
#     module use /projappl/project_2003755/Lazypipe-db/modulefiles/projects
#     module load lazypipe/3.1
#     tests/t2_databases.sh                # TAP on stdout
#     tests/t2_databases.sh | grep -v ^#   # results only
#
# Which databases a site must have depends on the strategies it intends to run,
# so DB-01 and DB-02 take their required set from the environment.  The
# defaults are the RefSeq-only minimum of §5.3, which is what Tiers 3 and 4
# need; set either to the empty string to check only that the listing works:
#
#     T2_REQUIRE_DBS="minimap.refseq.vi blastn.refseq.vi" tests/t2_databases.sh
#     T2_REQUIRE_FILTERS="" tests/t2_databases.sh
#
# Exit status = number of failed tests (0 = all good).  Neither skips nor TODOs
# are failures: skips mark absent optional tooling, TODOs mark known defects in
# LazypipeX itself that are documented in docs/testing_roadmap.md §12.

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

LZ="$INSTALL/lazypipe.pl"
INSTALL_DB="$INSTALL/perl/install_db.pl"
LINT="$TESTS_DIR/lib/config_lint.pl"

# Query fixtures for the smoke tests.  Both were cut from databases the pipeline
# ships against — a phiX174 fragment (in every nucleotide set) and two viral
# major-capsid proteins from UniRef100 — so a zero-hit result means the database
# is wrong, not the query.
FIXVI="$TESTS_DIR/fixtures/tiny.vi.fa"
FIXAA="$TESTS_DIR/fixtures/tiny.orfs.aa.fa"

# A tool counts as available only if it is also executable.  `command -v`
# alone is not enough: bash returns the path of a non-executable file found
# on PATH, so a downloaded-but-not-chmod+x binary was reported as installed
# and then failed at run time with "Permission denied".
have() {
	local p
	p=$( type -P "$1" 2>/dev/null ) || return 1
	[ -n "$p" ] && [ -x "$p" ]
}

# lazypipe.pl and install_db.pl both prefer ./config.yaml over the installed
# one, so run from the install dir and lint exactly the file they will read.
cd "$INSTALL" || exit 99
CFG="$INSTALL/config.yaml"

# The RefSeq-only minimum of docs/testing_roadmap.md §5.3.  Taxonomy is not
# listed here: it is not an ann.databases key and --databases cannot report it
# (DB-10 covers it).
: "${T2_REQUIRE_DBS:=minimap.refseq.vi minimap.refseq.abv blastn.refseq.vi blastn.refseq.abv}"
: "${T2_REQUIRE_FILTERS:=Homo_sapiens Neovison_vison}"

# Query smoke tests (DB-16, DB-19, DB-21, DB-23, DB-25) run only against
# databases at or below this on-disk size.  They are not optional in principle —
# they are the only checks that prove a database answers a query — but their cost
# scales with the database: minimap2 given a 910 GB core_nt FASTA builds an index
# before it aligns anything, and blastn against core_nt did not finish inside
# three minutes here, while the same query against blastn.refseq.vi took one
# second.  The default admits one database per engine (the RefSeq/viral ones plus
# a host genome) and keeps the tier inside its 2-10 min budget; raise it to cover
# the large sets, or set it to 0 to skip every query test.
: "${T2_SMOKE_MAX_MB:=8192}"
: "${T2_SMOKE_TIMEOUT:=300}"

# nodes.dmp older than this many days is reported but does not fail (DB-13).
: "${T2_TAXONOMY_MAX_AGE_DAYS:=183}"

WORK=$( mktemp -d "${TMPDIR:-/tmp}/lazytest-t2.XXXXXX" ) || exit 99
trap 'rm -rf "$WORK"' EXIT

# ----------------------------------------------------------------- report ---

tap_init "LazypipeX Tier 2 — database installation tests"
diag "host        : $( hostname )"
diag "date        : $( date -Is )"
diag "repo        : $REPO"
diag "install dir : $INSTALL"
diag "config      : $CFG"
for _v in databases hostgenomes data taxonomy_ncbi; do
	diag "\$$_v = ${!_v:-<unset>}"
done
if [ "$INSTALL" != "$REPO" ]; then
	diag "NOTE: \$LAZYPIPE_INSTALL_DIR differs from this checkout; testing the installed copy."
fi
diag ""

if [ ! -f "$CFG" ]; then
	# Without a config there is nothing any of these tests can say.  Report it
	# once, as a failure, rather than six times.
	fail DB-00 "config.yaml is readable" "no config file at $CFG"
	diag ""
	tap_done
	exit $?
fi

# Keys printed by --databases / --filters: one per line, at column 0, with a
# trailing colon.  The indented detail lines are skipped by the anchor.
listed_keys() { printf '%s\n' "$_OUT" | sed -n 's/^\([A-Za-z][A-Za-z0-9._+-]*\):$/\1/p'; }

# Names in $2 that listing $1 did not report.
missing_from() {
	local have="$1" want="$2" miss="" w
	for w in $want; do
		printf '%s\n' "$have" | grep -qx -- "$w" || miss="$miss $w"
	done
	printf '%s' "${miss# }"
}

# ============================================ DB-01 reference DB listing =====

try perl "$LZ" --databases
if [ "$_RC" -ne 0 ]; then
	fail DB-01 "installed reference databases are listed" \
		"perl lazypipe.pl --databases exited $_RC" \
		"$( printf '%s' "$_OUT" | head -5 )"
elif ! printf '%s' "$_OUT" | grep -q 'Listing Installed Reference Databases'; then
	# The listing is guarded by `defined($opt{'ann.databases'})`; if that section
	# is missing lazypipe.pl does not list, does not complain, and carries on
	# into a normal run.  Exit 0 alone is therefore not evidence of anything.
	fail DB-01 "installed reference databases are listed" \
		"--databases printed no listing header" \
		"ann.databases is probably missing from $CFG"
else
	db_keys=$( listed_keys )
	db_n=$( printf '%s' "$db_keys" | grep -c . )
	if [ "$db_n" -eq 0 ]; then
		fail DB-01 "installed reference databases are listed" \
			"--databases reported no installed database at all" \
			"every ann.databases path globbed empty; check \$databases and the install"
	else
		db_miss=$( missing_from "$db_keys" "$T2_REQUIRE_DBS" )
		if [ -z "$db_miss" ]; then
			pass DB-01 "installed reference databases are listed ($db_n installed)"
		else
			fail DB-01 "installed reference databases are listed ($db_n installed)" \
				"required but not installed: $db_miss" \
				"required set is \$T2_REQUIRE_DBS (default: the §5.3 RefSeq minimum)"
		fi
	fi
fi

# =============================================== DB-02 filter listing ========

try perl "$LZ" --filters
if [ "$_RC" -ne 0 ]; then
	fail DB-02 "installed background filters are listed" \
		"perl lazypipe.pl --filters exited $_RC" \
		"$( printf '%s' "$_OUT" | head -5 )"
elif ! printf '%s' "$_OUT" | grep -q 'Listing Installed Background Filters'; then
	fail DB-02 "installed background filters are listed" \
		"--filters printed no listing header" \
		"host.databases is probably missing from $CFG"
else
	flt_keys=$( listed_keys )
	flt_n=$( printf '%s' "$flt_keys" | grep -c . )
	if [ "$flt_n" -eq 0 ]; then
		fail DB-02 "installed background filters are listed" \
			"--filters reported no installed filter at all" \
			"check \$hostgenomes and that bwa index ran"
	else
		flt_miss=$( missing_from "$flt_keys" "$T2_REQUIRE_FILTERS" )
		if [ -z "$flt_miss" ]; then
			pass DB-02 "installed background filters are listed ($flt_n installed)"
			# --filters accepts a filter on .amb/.ann/.bwt alone, so a listed
			# filter is not a complete BWA index.  DB-24 is what proves it.
			diag "  note: --filters checks only .amb/.ann/.bwt; .pac/.sa are DB-24's job"
		else
			fail DB-02 "installed background filters are listed ($flt_n installed)" \
				"required but not installed: $flt_miss" \
				"required set is \$T2_REQUIRE_FILTERS"
		fi
	fi
fi

# ================================================= DB-03..05 config lint =====

# Run one lint mode.  Sets $_FAILS and $_DRIFTS to the findings of each
# severity; returns 1 if the linter itself could not run, in which case $_OUT
# holds its error.
run_lint() {
	_FAILS=""; _DRIFTS=""
	try perl "$LINT" "$1" "$CFG"
	[ "$_RC" -eq 0 ] || return 1
	_FAILS=$(  printf '%s\n' "$_OUT" | grep '^FAIL'  | cut -f2- )
	_DRIFTS=$( printf '%s\n' "$_OUT" | grep '^DRIFT' | cut -f2- )
	return 0
}

# A lint result becomes a failure, a TODO or a pass: FAIL findings mean this
# installation is broken, DRIFT findings mean LazypipeX itself is (see §12).
report_lint() {
	local id="$1" desc="$2" note="$3"
	if [ -n "$_FAILS" ]; then
		# Drift findings are reported alongside, not swallowed: a broken install
		# does not make the known config defects go away.
		fail "$id" "$desc" "$_FAILS" "$_DRIFTS"
	elif [ -n "$_DRIFTS" ]; then
		todof "$id" "$desc" "$note" "$_DRIFTS"
	else
		pass "$id" "$desc"
	fi
}

if [ ! -f "$LINT" ]; then
	skipt DB-03 "config lint: database paths"    "missing $LINT"
	skipt DB-04 "config lint: required fields"   "missing $LINT"
	skipt DB-05 "config lint: strategies resolve" "missing $LINT"
else
	if run_lint paths; then
		report_lint DB-03 "config lint: database paths" \
			"config.yaml drift, see docs/testing_roadmap.md §12"
	else
		fail DB-03 "config lint: database paths" "$_OUT"
	fi

	if run_lint fields; then
		report_lint DB-04 "config lint: required fields" \
			"config.yaml drift, see docs/testing_roadmap.md §12"
	else
		fail DB-04 "config lint: required fields" "$_OUT"
	fi

	if run_lint strategies; then
		report_lint DB-05 "config lint: strategies resolve" \
			"config.yaml drift, see docs/testing_roadmap.md §12 item 9"
	else
		fail DB-05 "config lint: strategies resolve" "$_OUT"
	fi
fi

# ============================================== DB-06 install_db contract ====

# Download integrity proper needs the network and the databases themselves; what
# is checkable offline is that a failed install is reportable at all.  A scripted
# install loop can only detect a bad key from the exit status.
#
# The key below is asserted absent from config.yaml before install_db.pl is run
# with it: this tier must never be able to start a real download.

BOGUS_KEY="lazypipe.tier2.no.such.database"

key_defined() {
	perl -MYAML::Tiny -e '
		my $y = YAML::Tiny->read($ARGV[0]) or exit 2;
		my %o = %{ $y->[0] };
		for my $sec ( "ann.databases", "host.databases" ){
			exit 0 if ref($o{$sec}) eq "HASH" && exists $o{$sec}->{ $ARGV[1] };
		}
		exit 1;
	' "$CFG" "$1" 2>/dev/null
}

if [ ! -f "$INSTALL_DB" ]; then
	skipt DB-06 "install_db.pl reports an unmatched --db as a failure" "missing $INSTALL_DB"
elif key_defined "$BOGUS_KEY"; then
	# Refuse to run rather than risk invoking a real install.
	skipt DB-06 "install_db.pl reports an unmatched --db as a failure" \
		"'$BOGUS_KEY' unexpectedly exists in $CFG"
else
	try perl "$INSTALL_DB" --db "$BOGUS_KEY"
	db06_rc="$_RC"
	db06_out="$_OUT"

	if ! printf '%s' "$db06_out" | grep -q 'did not match any database'; then
		fail DB-06 "install_db.pl reports an unmatched --db as a failure" \
			"expected a 'did not match any database' message, got (rc=$db06_rc):" \
			"$( printf '%s' "$db06_out" | head -5 )"
	elif [ "$db06_rc" -ne 0 ]; then
		pass DB-06 "install_db.pl reports an unmatched --db as a failure (rc=$db06_rc)"
	else
		todof DB-06 "install_db.pl reports an unmatched --db as a failure" \
			"known defect, docs/testing_roadmap.md §12 item 3" \
			"install_db.pl printed 'did not match any database' and exited 0" \
			"a scripted install loop cannot detect the failure and installs nothing"
	fi

fi

# ######################################################################### #
#                  §5.2 — per-engine integrity checks                       #
#                          DB-10 … DB-27                                    #
# ######################################################################### #

# Everything below runs against the databases this site actually has.  The
# inventory resolves config.yaml's $ENV paths and marks a database installed by
# lazypipe.pl's own rule (glob("$db*") is non-empty), so these checks cover
# exactly the set DB-01/DB-02 reported and never invent a database to fail on.

INVENTORY="$TESTS_DIR/lib/db_inventory.pl"
INV="$WORK/inventory.tsv"

if [ ! -f "$INVENTORY" ]; then
	inv_ok=0
	diag "NOTE: missing $INVENTORY — DB-10..DB-27 skipped"
elif ! perl "$INVENTORY" "$CFG" > "$INV" 2>"$WORK/inv.err"; then
	inv_ok=0
	diag "NOTE: $INVENTORY failed — DB-10..DB-27 skipped"
	diag "$( head -3 "$WORK/inv.err" )"
else
	inv_ok=1
fi

# "key<TAB>path" for every installed database of one engine.
inv_of() { awk -F'\t' -v e="$1" '$3==e && $5==1 { print $2"\t"$4 }' "$INV"; }

# Resolved path of a single inventory row, by key.
inv_path() { awk -F'\t' -v k="$1" '$2==k && $5==1 { print $4; exit }' "$INV"; }

# On-disk megabytes of a database: every file the install produced, which is
# what glob("$db*") means everywhere else in the pipeline.
db_size_mb() {
	local s
	s=$( du -scm -- "$1"* 2>/dev/null | tail -1 | awk '{ print $1 }' )
	printf '%s' "${s:-0}"
}

# Whether a query smoke test may run against this database (see T2_SMOKE_MAX_MB).
smoke_allowed() {
	[ "$T2_SMOKE_MAX_MB" -gt 0 ] || return 1
	[ "$( db_size_mb "$1" )" -le "$T2_SMOKE_MAX_MB" ]
}

# Report an aggregate check: $3 empty means every database passed.
report_agg() {
	local id="$1" desc="$2" bad="$3" hint="${4:-}"
	if [ -z "$bad" ]; then
		pass "$id" "$desc"
	else
		fail "$id" "$desc" "$( printf '%b' "$bad" )" "$hint"
	fi
}

TAXDIR=""
[ "$inv_ok" -eq 1 ] && TAXDIR=$( inv_path taxonomy )

# ================================================== DB-10..13 taxonomy =======

if [ "$inv_ok" -eq 0 ]; then
	skipt DB-10 "NCBI taxonomy dump files are present and non-empty" "no inventory"
elif [ -z "$TAXDIR" ] || [ ! -d "$TAXDIR" ]; then
	fail DB-10 "NCBI taxonomy dump files are present and non-empty" \
		"taxonomy.db does not resolve to a directory: ${TAXDIR:-<unset>}" \
		"every taxid-aware step depends on it; check \$taxonomy_ncbi"
else
	t10_bad=""
	for f in nodes.dmp names.dmp merged.dmp delnodes.dmp division.dmp; do
		if [ ! -e "$TAXDIR/$f" ]; then
			t10_bad="$t10_bad\n  $f: missing"
		elif [ ! -s "$TAXDIR/$f" ]; then
			t10_bad="$t10_bad\n  $f: present but empty"
		fi
	done
	report_agg DB-10 "NCBI taxonomy dump files are present and non-empty" "$t10_bad" \
		"re-run: perl perl/install_db.pl --db taxonomy -v"
fi

# 1239574 is Mamastrovirus 10 / mink astrovirus — the taxid Tier 4 expects to
# find in M15small, so a taxonomy that cannot resolve it cannot pass Tier 4.
TAXID_PROBE=1239574

if ! have taxonkit; then
	skipt DB-11 "taxonomy resolves a known taxid to a full lineage" "taxonkit not on PATH"
	skipt DB-12 "taxonomy reformat yields the ranks the pipeline asks for" "taxonkit not on PATH"
elif [ -z "$TAXDIR" ] || [ ! -s "$TAXDIR/nodes.dmp" ]; then
	skipt DB-11 "taxonomy resolves a known taxid to a full lineage" "no usable taxonomy dir"
	skipt DB-12 "taxonomy reformat yields the ranks the pipeline asks for" "no usable taxonomy dir"
else
	try_sh "echo $TAXID_PROBE | taxonkit lineage --data-dir '$TAXDIR'"
	lineage=$( printf '%s\n' "$_OUT" | head -1 | cut -f2 )
	if [ "$_RC" -ne 0 ]; then
		fail DB-11 "taxonomy resolves a known taxid to a full lineage" \
			"taxonkit lineage exited $_RC" "$( printf '%s' "$_OUT" | head -3 )"
	elif [ -z "$lineage" ]; then
		fail DB-11 "taxonomy resolves a known taxid to a full lineage" \
			"taxid $TAXID_PROBE resolved to an empty lineage" \
			"the dump is readable but does not contain this taxid; it is probably truncated or stale"
	elif ! printf '%s' "$lineage" | grep -q ';'; then
		fail DB-11 "taxonomy resolves a known taxid to a full lineage" \
			"taxid $TAXID_PROBE resolved to a single rank: $lineage"
	else
		pass DB-11 "taxonomy resolves a known taxid to a full lineage"
	fi

	# The pipeline asks for {s} {g} {f} {o} {c} {k}; species, genus and family
	# are the ones the reports are keyed on and must not come back empty or NA.
	try_sh "echo $TAXID_PROBE | taxonkit reformat --data-dir '$TAXDIR' -I 1 -f '{s}\t{g}\t{f}\t{o}\t{c}\t{k}'"
	if [ "$_RC" -ne 0 ]; then
		fail DB-12 "taxonomy reformat yields the ranks the pipeline asks for" \
			"taxonkit reformat exited $_RC" "$( printf '%s' "$_OUT" | head -3 )"
	else
		t12_bad=""
		i=2
		for rank in species genus family; do
			v=$( printf '%s\n' "$_OUT" | head -1 | cut -f$i )
			case "$v" in ""|NA) t12_bad="$t12_bad\n  $rank: '$v'" ;; esac
			i=$(( i + 1 ))
		done
		report_agg DB-12 "taxonomy reformat yields the ranks the pipeline asks for" "$t12_bad" \
			"empty species/genus/family make every annotation report lose its binning"
	fi
fi

if [ -z "$TAXDIR" ] || [ ! -e "$TAXDIR/nodes.dmp" ]; then
	fail DB-13 "taxonomy is present and its age is known" \
		"no nodes.dmp under ${TAXDIR:-<unset>}"
else
	# Age is reported, never fatal: an old taxonomy still works, it just stops
	# resolving recently added taxids.  config.yaml sets update_time 5 with
	# update 0, so nothing refreshes it automatically.
	age_days=$( perl -e 'printf "%d", -M $ARGV[0]' "$TAXDIR/nodes.dmp" 2>/dev/null )
	age_days=${age_days:-0}
	if [ "$age_days" -gt "$T2_TAXONOMY_MAX_AGE_DAYS" ]; then
		pass DB-13 "taxonomy is present and its age is known ($age_days days)"
		diag "  note: nodes.dmp is $age_days days old (> $T2_TAXONOMY_MAX_AGE_DAYS)"
		diag "        refresh with: perl perl/install_db.pl --db taxonomy --force"
	else
		pass DB-13 "taxonomy is present and its age is known ($age_days days)"
	fi
fi

# ================================================= DB-14..16 BLAST ===========

blast_list=""
[ "$inv_ok" -eq 1 ] && blast_list=$( { inv_of blastn; inv_of blastp; inv_of blastx; } | sort -u )

if [ "$inv_ok" -eq 0 ]; then
	skipt DB-14 "BLAST indices report a positive sequence count" "no inventory"
	skipt DB-15 "BLAST indices resolve taxids" "no inventory"
	skipt DB-16 "BLAST answers a query with at least one hit" "no inventory"
elif ! have blastdbcmd; then
	skipt DB-14 "BLAST indices report a positive sequence count" "blastdbcmd not on PATH"
	skipt DB-15 "BLAST indices resolve taxids" "blastdbcmd not on PATH"
	skipt DB-16 "BLAST answers a query with at least one hit" "blastdbcmd not on PATH"
elif [ -z "$blast_list" ]; then
	skipt DB-14 "BLAST indices report a positive sequence count" "no BLAST database installed"
	skipt DB-15 "BLAST indices resolve taxids" "no BLAST database installed"
	skipt DB-16 "BLAST answers a query with at least one hit" "no BLAST database installed"
else
	b14_bad=""; b15_bad=""; b16_bad=""; b16_ran=0; b16_skipped=""
	n_blast=0
	while IFS="$( printf '\t' )" read -r k p; do
		[ -n "${k:-}" ] || continue
		n_blast=$(( n_blast + 1 ))

		# DB-14 — the index opens and reports sequences.  A multi-volume set is
		# addressed through its .nal/.pal alias, and -info fails if any volume
		# named by the alias is missing, so this covers volume resolution too.
		try blastdbcmd -db "$p" -info
		if [ "$_RC" -ne 0 ]; then
			b14_bad="$b14_bad\n  $k: blastdbcmd -info exited $_RC: $( printf '%s' "$_OUT" | head -1 | cut -c1-80 )"
			continue
		fi
		nseq=$( printf '%s\n' "$_OUT" | sed -n 's/.*[^0-9,]\([0-9,]\{1,\}\) sequences;.*/\1/p' | head -1 | tr -d ',' )
		if [ -z "$nseq" ]; then
			b14_bad="$b14_bad\n  $k: -info printed no sequence count"
		elif [ "$nseq" -le 0 ]; then
			b14_bad="$b14_bad\n  $k: $nseq sequences"
		fi

		# DB-15 — without taxdb.btd/.bti next to the index every staxid comes
		# back 0 and downstream binning degrades silently rather than failing.
		try_sh "blastdbcmd -db '$p' -entry all -outfmt '%T' 2>/dev/null | head -20"
		taxids=$( printf '%s\n' "$_OUT" | grep -c '^[0-9]\{1,\}$' )
		nonzero=$( printf '%s\n' "$_OUT" | grep -c '^[1-9][0-9]*$' )
		if [ "$taxids" -eq 0 ]; then
			b15_bad="$b15_bad\n  $k: -entry all returned no taxids"
		elif [ "$nonzero" -eq 0 ]; then
			b15_bad="$b15_bad\n  $k: every sampled taxid is 0 (taxdb.btd/.bti missing beside the index?)"
		fi

		# DB-16 — the database answers a real query.  The nucleotide fixture is a
		# phiX174 fragment cut from a shipped viral set and the protein fixture
		# comes from UniRef100, so no hit means the index is broken rather than
		# the query being a poor match.
		if [ ! -s "$FIXVI" ]; then
			b16_skipped="$b16_skipped $k"
		elif ! smoke_allowed "$p"; then
			b16_skipped="$b16_skipped $k($( db_size_mb "$p" )MB)"
		else
			case "$k" in
				blastp.*|blastx.*) prog=blastp; q="$FIXAA" ;;
				*)                 prog=blastn; q="$FIXVI" ;;
			esac
			if [ ! -s "$q" ]; then
				b16_skipped="$b16_skipped $k"
			else
				b16_ran=$(( b16_ran + 1 ))
				TEST_TIMEOUT="$T2_SMOKE_TIMEOUT" \
					try_sh "$prog -db '$p' -query '$q' -max_target_seqs 5 -outfmt 6 2>/dev/null"
				nhit=$( printf '%s\n' "$_OUT" | grep -c '	' )
				if [ "$_RC" -ne 0 ]; then
					b16_bad="$b16_bad\n  $k: $prog exited $_RC"
				elif [ "$nhit" -lt 1 ]; then
					b16_bad="$b16_bad\n  $k: $prog returned 0 hits for $( basename "$q" )"
				fi
			fi
		fi
	done <<EOF
$blast_list
EOF

	report_agg DB-14 "BLAST indices report a positive sequence count ($n_blast checked)" "$b14_bad" \
		"a index that will not open cannot be searched, whichever strategy names it"
	report_agg DB-15 "BLAST indices resolve taxids ($n_blast checked)" "$b15_bad" \
		"taxid 0 empties the staxid column and silently degrades all downstream binning"
	if [ "$b16_ran" -eq 0 ]; then
		skipt DB-16 "BLAST answers a query with at least one hit" \
			"no BLAST database within \$T2_SMOKE_MAX_MB=${T2_SMOKE_MAX_MB}MB"
	else
		report_agg DB-16 "BLAST answers a query with at least one hit ($b16_ran of $n_blast queried)" "$b16_bad"
		[ -n "$b16_skipped" ] && diag "  not queried (over \$T2_SMOKE_MAX_MB):$b16_skipped"
	fi
fi

# =============================================== DB-17..19 minimap2 ==========

mm_list=""
[ "$inv_ok" -eq 1 ] && mm_list=$( inv_of minimap )

if [ "$inv_ok" -eq 0 ] || [ -z "$mm_list" ]; then
	skipt DB-17 "minimap2 databases have their acc2taxid sidecar" "no minimap2 database installed"
	skipt DB-18 "acc2taxid maps are well-formed" "no minimap2 database installed"
	skipt DB-19 "minimap2 answers a query with at least one alignment" "no minimap2 database installed"
else
	m17_bad=""; m18_bad=""; m19_bad=""; m19_ran=0; m19_skipped=""
	n_mm=0
	while IFS="$( printf '\t' )" read -r k p; do
		[ -n "${k:-}" ] || continue
		n_mm=$(( n_mm + 1 ))

		# DB-17 — SeqAn.pm derives the sidecar name from the db name: X.mmi ->
		# X.acc2taxid, otherwise X.acc2taxid, and prefers X.mmi when it exists.
		case "$p" in
			*.mmi) acc="${p%.mmi}.acc2taxid" ;;
			*)     acc="$p.acc2taxid" ;;
		esac
		target="$p"
		[ -e "$p.mmi" ] && target="$p.mmi"
		if [ ! -e "$target" ]; then
			m17_bad="$m17_bad\n  $k: no database file at $target"
		fi
		if [ ! -e "$acc" ]; then
			m17_bad="$m17_bad\n  $k: no acc2taxid sidecar at $acc"
			continue
		elif [ ! -s "$acc" ]; then
			m17_bad="$m17_bad\n  $k: acc2taxid sidecar is empty"
			continue
		fi

		# DB-18 — two columns, accession then a numeric taxid, one row per
		# sequence.  The row count is compared against seqkit's own .stats
		# sidecar when the install wrote one; reading a 900 GB FASTA to count
		# records is not something a test tier can do.
		ncol=$( head -1 "$acc" | awk -F'\t' '{ print NF }' )
		if [ "${ncol:-0}" -lt 2 ]; then
			m18_bad="$m18_bad\n  $k: acc2taxid has $ncol tab-separated column(s), expected 2"
			continue
		fi
		badrows=$( head -1000 "$acc" | awk -F'\t' '$1=="" || $2 !~ /^[0-9]+$/ { n++ } END { print n+0 }' )
		if [ "$badrows" -gt 0 ]; then
			m18_bad="$m18_bad\n  $k: $badrows of the first 1000 acc2taxid rows are not accession<TAB>taxid"
		fi
		rows=$( wc -l < "$acc" )
		if [ -s "$p.stats" ]; then
			nseq=$( awk 'NR==2 { gsub(/,/,"",$4); print $4 }' "$p.stats" )
			if [ -n "$nseq" ] && [ "$nseq" -gt 0 ]; then
				# 99 % of the FASTA's accessions must be in the map.
				need=$(( nseq * 99 / 100 ))
				if [ "$rows" -lt "$need" ]; then
					m18_bad="$m18_bad\n  $k: acc2taxid has $rows rows for $nseq sequences (< 99 %)"
				fi
			fi
		fi

		# DB-19 — minimap2 given a FASTA indexes it before aligning, so this is
		# the most size-sensitive check in the tier.
		if [ ! -s "$FIXVI" ]; then
			m19_skipped="$m19_skipped $k"
		elif ! have minimap2; then
			m19_skipped="$m19_skipped $k"
		elif ! smoke_allowed "$target"; then
			m19_skipped="$m19_skipped $k($( db_size_mb "$target" )MB)"
		else
			m19_ran=$(( m19_ran + 1 ))
			TEST_TIMEOUT="$T2_SMOKE_TIMEOUT" \
				try_sh "minimap2 -x asm20 --secondary=yes --cs -s 200 '$target' '$FIXVI' 2>/dev/null"
			npaf=$( printf '%s\n' "$_OUT" | grep -c '	' )
			if [ "$_RC" -ne 0 ]; then
				m19_bad="$m19_bad\n  $k: minimap2 exited $_RC"
			elif [ "$npaf" -lt 1 ]; then
				m19_bad="$m19_bad\n  $k: no PAF line for $( basename "$FIXVI" )"
			fi
		fi
	done <<EOF
$mm_list
EOF

	report_agg DB-17 "minimap2 databases have their acc2taxid sidecar ($n_mm checked)" "$m17_bad" \
		"annotate_minimap dies with 'no acc2taxid map' (SeqAn.pm:906) when it is absent"
	report_agg DB-18 "acc2taxid maps are well-formed ($n_mm checked)" "$m18_bad"
	if [ "$m19_ran" -eq 0 ]; then
		skipt DB-19 "minimap2 answers a query with at least one alignment" \
			"no minimap2 database within \$T2_SMOKE_MAX_MB=${T2_SMOKE_MAX_MB}MB, or minimap2 absent"
	else
		report_agg DB-19 "minimap2 answers a query with at least one alignment ($m19_ran of $n_mm queried)" "$m19_bad"
		[ -n "$m19_skipped" ] && diag "  not queried (over \$T2_SMOKE_MAX_MB):$m19_skipped"
	fi
fi

# ================================================= DB-20..21 DIAMOND ========

# diamondp and diamondx point at the same .dmnd file, so the list is deduped on
# the path: testing one file twice would only make the report longer.
dmnd_list=""
[ "$inv_ok" -eq 1 ] && dmnd_list=$( { inv_of diamondp; inv_of diamondx; } | sort -u -k2 )

if [ "$inv_ok" -eq 0 ] || [ -z "$dmnd_list" ]; then
	skipt DB-20 "DIAMOND databases report a positive sequence count" "no DIAMOND database installed"
	skipt DB-21 "DIAMOND answers a query with at least one hit" "no DIAMOND database installed"
elif ! have diamond; then
	skipt DB-20 "DIAMOND databases report a positive sequence count" "diamond not on PATH"
	skipt DB-21 "DIAMOND answers a query with at least one hit" "diamond not on PATH"
else
	d20_bad=""; d21_bad=""; d21_ran=0; d21_skipped=""
	n_dmnd=0
	while IFS="$( printf '\t' )" read -r k p; do
		[ -n "${k:-}" ] || continue
		n_dmnd=$(( n_dmnd + 1 ))

		# dbinfo also fails when the .dmnd format version predates the installed
		# binary, which is the usual symptom after a diamond upgrade.
		try diamond dbinfo --db "$p"
		if [ "$_RC" -ne 0 ]; then
			d20_bad="$d20_bad\n  $k: diamond dbinfo exited $_RC: $( printf '%s' "$_OUT" | grep -i error | head -1 | cut -c1-80 )"
			continue
		fi
		nseq=$( printf '%s\n' "$_OUT" | sed -n 's/^ *Sequences  *\([0-9]\{1,\}\).*/\1/p' | head -1 )
		if [ -z "$nseq" ]; then
			d20_bad="$d20_bad\n  $k: dbinfo printed no sequence count"
		elif [ "$nseq" -le 0 ]; then
			d20_bad="$d20_bad\n  $k: $nseq sequences"
		fi

		if [ ! -s "$FIXAA" ]; then
			d21_skipped="$d21_skipped $k"
		elif ! smoke_allowed "$p"; then
			d21_skipped="$d21_skipped $k($( db_size_mb "$p" )MB)"
		else
			d21_ran=$(( d21_ran + 1 ))
			# No -o: diamond puts its temporary file next to the output path, so
			# "-o /dev/stdout" makes it try to write into /dev and exit 1.
			TEST_TIMEOUT="$T2_SMOKE_TIMEOUT" \
				try_sh "diamond blastp --db '$p' --query '$FIXAA' --very-sensitive --quiet 2>/dev/null"
			nhit=$( printf '%s\n' "$_OUT" | grep -c '	' )
			if [ "$_RC" -ne 0 ]; then
				d21_bad="$d21_bad\n  $k: diamond blastp exited $_RC"
			elif [ "$nhit" -lt 1 ]; then
				d21_bad="$d21_bad\n  $k: 0 hits for $( basename "$FIXAA" )"
			fi
		fi
	done <<EOF
$dmnd_list
EOF

	report_agg DB-20 "DIAMOND databases report a positive sequence count ($n_dmnd checked)" "$d20_bad"
	if [ "$d21_ran" -eq 0 ]; then
		skipt DB-21 "DIAMOND answers a query with at least one hit" \
			"no DIAMOND database within \$T2_SMOKE_MAX_MB=${T2_SMOKE_MAX_MB}MB"
	else
		report_agg DB-21 "DIAMOND answers a query with at least one hit ($d21_ran of $n_dmnd queried)" "$d21_bad"
		[ -n "$d21_skipped" ] && diag "  not queried (over \$T2_SMOKE_MAX_MB):$d21_skipped"
	fi
fi

# ==================================================== DB-22..23 HMM =========

hmm_list=""
[ "$inv_ok" -eq 1 ] && hmm_list=$( inv_of hmmscan )

if [ "$inv_ok" -eq 0 ] || [ -z "$hmm_list" ]; then
	skipt DB-22 "HMM databases are pressed and report a positive profile count" "no HMM database installed"
	skipt DB-23 "hmmscan answers a query" "no HMM database installed"
elif ! have hmmstat; then
	skipt DB-22 "HMM databases are pressed and report a positive profile count" "hmmstat not on PATH"
	skipt DB-23 "hmmscan answers a query" "hmmscan not on PATH"
else
	h22_bad=""; h23_bad=""; h23_ran=0; h23_skipped=""
	n_hmm=0
	while IFS="$( printf '\t' )" read -r k p; do
		[ -n "${k:-}" ] || continue
		n_hmm=$(( n_hmm + 1 ))

		# hmmscan needs the pressed binaries, not just the .hmm text file.
		for ext in h3f h3i h3m h3p; do
			[ -s "$p.$ext" ] || h22_bad="$h22_bad\n  $k: missing or empty $( basename "$p" ).$ext (run hmmpress)"
		done

		try hmmstat "$p"
		if [ "$_RC" -ne 0 ]; then
			h22_bad="$h22_bad\n  $k: hmmstat exited $_RC"
			continue
		fi
		nprof=$( printf '%s\n' "$_OUT" | grep -c '^ *[0-9]' )
		if [ "$nprof" -lt 1 ]; then
			h22_bad="$h22_bad\n  $k: hmmstat reported no profiles"
		fi

		if [ ! -s "$FIXAA" ] || ! have hmmscan; then
			h23_skipped="$h23_skipped $k"
		elif ! smoke_allowed "$p"; then
			h23_skipped="$h23_skipped $k($( db_size_mb "$p" )MB)"
		else
			h23_ran=$(( h23_ran + 1 ))
			TEST_TIMEOUT="$T2_SMOKE_TIMEOUT" \
				try_sh "hmmscan -E 0.01 '$p' '$FIXAA' > /dev/null 2>&1"
			# Zero hits is a legitimate result for a small fixture against a
			# specialised profile set, so only the exit status is asserted.
			[ "$_RC" -eq 0 ] || h23_bad="$h23_bad\n  $k: hmmscan exited $_RC"
		fi
	done <<EOF
$hmm_list
EOF

	report_agg DB-22 "HMM databases are pressed and report a positive profile count ($n_hmm checked)" "$h22_bad"
	if [ "$h23_ran" -eq 0 ]; then
		skipt DB-23 "hmmscan answers a query" "no HMM database within \$T2_SMOKE_MAX_MB, or hmmscan absent"
	else
		report_agg DB-23 "hmmscan answers a query ($h23_ran of $n_hmm queried)" "$h23_bad"
		[ -n "$h23_skipped" ] && diag "  not queried:$h23_skipped"
	fi
fi

# ============================================ DB-24..25 background filters ===

host_list=""
[ "$inv_ok" -eq 1 ] && host_list=$( inv_of hostgen )

if [ "$inv_ok" -eq 0 ] || [ -z "$host_list" ]; then
	skipt DB-24 "every installed filter has a complete BWA index" "no background filter installed"
	skipt DB-25 "BWA aligns sample reads against the matching host" "no background filter installed"
else
	f24_bad=""; n_host=0
	while IFS="$( printf '\t' )" read -r k p; do
		[ -n "${k:-}" ] || continue
		n_host=$(( n_host + 1 ))
		# --filters checks only .amb/.ann/.bwt, so a filter missing .pac or .sa
		# is listed as installed and then fails inside bwa at run time.
		for ext in amb ann bwt pac sa; do
			[ -s "$p.$ext" ] || f24_bad="$f24_bad\n  $k: missing or empty .$ext"
		done
	done <<EOF
$host_list
EOF
	report_agg DB-24 "every installed filter has a complete BWA index ($n_host checked)" "$f24_bad" \
		"a set with only .amb/.ann/.bwt passes --filters and then fails inside bwa"

	# DB-25 uses the bundled mink library against the mink filter: a host filter
	# that aligns nothing is indistinguishable from no filter at all.
	SMOKE_HOST="${T2_SMOKE_HOST:-Neovison_vison}"
	hostp=$( inv_path "$SMOKE_HOST" )
	reads="$REPO/data/samples/M15small_R1.fastq"
	if ! have bwa; then
		skipt DB-25 "BWA aligns sample reads against the matching host" "bwa not on PATH"
	elif [ -z "$hostp" ]; then
		skipt DB-25 "BWA aligns sample reads against the matching host" "$SMOKE_HOST not installed"
	elif [ ! -s "$reads" ]; then
		skipt DB-25 "BWA aligns sample reads against the matching host" "missing $reads"
	elif ! smoke_allowed "$hostp"; then
		skipt DB-25 "BWA aligns sample reads against the matching host" \
			"$SMOKE_HOST is $( db_size_mb "$hostp" )MB, over \$T2_SMOKE_MAX_MB=${T2_SMOKE_MAX_MB}MB"
	else
		head -4000 "$reads" > "$WORK/smoke_R1.fastq"	# 1000 reads
		TEST_TIMEOUT="$T2_SMOKE_TIMEOUT" \
			try_sh "bwa mem -t 2 '$hostp' '$WORK/smoke_R1.fastq' 2>/dev/null | awk '\$1 !~ /^@/ && \$3 != \"*\"' | head -5"
		naln=$( printf '%s\n' "$_OUT" | grep -c '	' )
		if [ "$_RC" -ne 0 ]; then
			fail DB-25 "BWA aligns sample reads against the matching host" \
				"bwa mem exited $_RC for $SMOKE_HOST"
		elif [ "$naln" -lt 1 ]; then
			fail DB-25 "BWA aligns sample reads against the matching host" \
				"no aligned record from 1000 M15small reads against $SMOKE_HOST" \
				"the index opens but matches nothing — suspect a truncated or wrong genome"
		else
			pass DB-25 "BWA aligns sample reads against the matching host ($SMOKE_HOST)"
		fi
	fi
fi

# ================================================== DB-26 size sanity ========

# The User Guide's tables give *download* sizes (the column is "Size *.gz")
# while what is on disk is the extracted database, and the tables describe an
# older database generation than the one config.yaml now installs (§12 item 9).
# Both make the roadmap's ~20 % band unusable: measured here, several current
# databases sit 10-20 % under their documented archive, which is generation
# drift and not a fault.  What survives both problems is the failure the check
# was written for — "a 10x shortfall means a truncated download" — so only a
# gross shortfall fails, and everything else is reported for the record.
#
# This is the one check that reads the documentation, because its reference
# values live only there.  Scraping a prose table is fragile by nature: if the
# tables are reformatted the parse yields nothing and the check skips rather
# than passing silently.  Moving the sizes into a data file would remove the
# scrape entirely and is the better long-term fix.
DOCS="$REPO/docs/UserGuide.v3.1.md"

if [ "$inv_ok" -eq 0 ]; then
	skipt DB-26 "no installed database falls 10x short of its documented size" "no inventory"
elif [ ! -f "$DOCS" ]; then
	skipt DB-26 "no installed database falls 10x short of its documented size" "missing $DOCS"
else
	# "|*key* | 5.6 GB | description" -> "key<TAB>megabytes"
	sed -n 's/^|\*\([A-Za-z0-9._]\{1,\}\)\* *| *\([0-9.]\{1,\}\) *\(GB\|MB\).*/\1 \2 \3/p' "$DOCS" \
		| awk '{ mb = ($3 == "GB") ? $2 * 1024 : $2; printf "%s\t%d\n", $1, mb }' \
		| sort -u > "$WORK/doc_sizes.tsv"

	s26_bad=""; s26_cmp=0; s26_nodoc=0; s26_under=""
	while IFS="$( printf '\t' )" read -r sec k eng p inst; do
		[ "${inst:-0}" = "1" ] || continue
		[ "$sec" = "ann" ] || continue
		doc_mb=$( awk -F'\t' -v k="$k" '$1==k { print $2; exit }' "$WORK/doc_sizes.tsv" )
		if [ -z "$doc_mb" ]; then
			s26_nodoc=$(( s26_nodoc + 1 ))
			continue
		fi
		s26_cmp=$(( s26_cmp + 1 ))
		on_mb=$( db_size_mb "$p" )
		if [ $(( on_mb * 10 )) -lt "$doc_mb" ]; then
			s26_bad="$s26_bad\n  $k: ${on_mb}MB on disk vs ${doc_mb}MB documented — a 10x shortfall, download truncated?"
		elif [ "$on_mb" -lt "$doc_mb" ]; then
			s26_under="$s26_under\n    $k: ${on_mb}MB on disk, ${doc_mb}MB documented"
		fi
	done < "$INV"

	if [ "$s26_cmp" -eq 0 ]; then
		skipt DB-26 "no installed database falls 10x short of its documented size" \
			"no installed database has a documented size (docs describe an older generation, §12 item 9)"
	else
		report_agg DB-26 "no installed database falls 10x short of its documented size ($s26_cmp compared)" "$s26_bad"
		if [ -n "$s26_under" ]; then
			diag "  under the documented archive size, within the drift the tables allow:"
			diag "$( printf '%b' "$s26_under" )"
		fi
		[ "$s26_nodoc" -gt 0 ] && diag "  $s26_nodoc installed database(s) have no size in the User Guide (§12 item 9)"
	fi
fi

# =============================================== DB-27 idempotent install ====

# Re-running an install that is already on disk must be a no-op.  install_db()
# globs the target and returns before reaching its wget when anything matches,
# so this is safe to run against a real key — but only one that is already
# installed, which is asserted from the inventory before the call.
if [ "$inv_ok" -eq 0 ] || [ ! -f "$INSTALL_DB" ]; then
	skipt DB-27 "re-installing an installed database changes nothing" "no inventory or no install_db.pl"
else
	idem_key=$( awk -F'\t' '$1=="ann" && $5==1 { print $2; exit }' "$INV" )
	idem_path=$( inv_path "$idem_key" )
	if [ -z "$idem_key" ] || [ -z "$idem_path" ]; then
		skipt DB-27 "re-installing an installed database changes nothing" "no installed database to re-install"
	else
		before=$( ls -l -- "$idem_path"* 2>/dev/null | md5sum )
		try perl "$INSTALL_DB" --db "$idem_key" -v
		after=$( ls -l -- "$idem_path"* 2>/dev/null | md5sum )
		if [ "$before" != "$after" ]; then
			fail DB-27 "re-installing an installed database changes nothing" \
				"$idem_key: files on disk changed during a --force-less re-install"
		elif ! printf '%s' "$_OUT" | grep -q 'skipping: found files on disk'; then
			fail DB-27 "re-installing an installed database changes nothing" \
				"$idem_key: expected 'skipping: found files on disk', got (rc=$_RC):" \
				"$( printf '%s' "$_OUT" | head -4 )"
		else
			pass DB-27 "re-installing an installed database changes nothing ($idem_key)"
		fi
	fi
fi

diag ""
tap_done
