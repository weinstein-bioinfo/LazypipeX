#!/usr/bin/env bash
#
# LazypipeX Tier 0 — environment and prerequisites.
# Implements ENV-01 … ENV-10 of docs/testing_roadmap.md §3.
#
# Needs no reference databases: it checks that the software is runnable and
# that config.yaml agrees with the environment it has been loaded into.
#
# Usage:
#     module use /projappl/project_2003755/Lazypipe-db/modulefiles/projects
#     module load lazypipe/3.1
#     tests/t0_environment.sh              # TAP on stdout
#     tests/t0_environment.sh | grep -v ^# # results only
#
# Exit status = number of failed tests (0 = all good).  Skips are not failures:
# they mark optional tools and Slurm-only checks that do not apply here.
#
# ENV-08 is reported as three tests so that an install missing only optional
# tools is distinguishable from one that cannot run at all:
#     ENV-08   tools required for a default `-p main` run   (hard fail)
#     ENV-08a  tools gated behind a non-default option      (skip + note)
#     ENV-08b  tools gated behind an annotation strategy    (skip + note)
#
# The script writes nothing outside $TMPDIR.

set -uo pipefail

TESTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# shellcheck source=lib/assert.sh
. "$TESTS_DIR/lib/assert.sh"

REPO=$( cd "$TESTS_DIR/.." && pwd )

# lazypipe.pl resolves its install dir from $LAZYPIPE_INSTALL_DIR, falling back
# to its own location; mirror that here so the tests exercise the installation
# the user will actually run.
if [ -n "${LAZYPIPE_INSTALL_DIR:-}" ] && [ -f "$LAZYPIPE_INSTALL_DIR/lazypipe.pl" ]; then
	INSTALL="$LAZYPIPE_INSTALL_DIR"
else
	INSTALL="$REPO"
fi

# And config.yaml is taken from the CURRENT directory if present, otherwise
# from the install dir (lazypipe.pl:35, install_db.pl:15).
if [ -f "$PWD/config.yaml" ]; then
	CONFIG="$PWD/config.yaml"
else
	CONFIG="$INSTALL/config.yaml"
fi

WORK=$( mktemp -d "${TMPDIR:-/tmp}/lazytest-t0.XXXXXX" ) || exit 99
trap 'rm -rf "$WORK"' EXIT

# ---------------------------------------------------------------- helpers ---

# Read a (possibly nested) scalar out of config.yaml.
cfg_get() {
	perl -MYAML::Tiny -e '
		my $y = YAML::Tiny->read($ARGV[0]) or exit 1;
		my $v = $y->[0];
		for my $k (@ARGV[1 .. $#ARGV]) {
			$v = (ref($v) eq "HASH") ? $v->{$k} : undef;
			last unless defined $v;
		}
		print defined($v) ? $v : "";
	' "$CONFIG" "$@" 2>/dev/null
}

# Expand $VAR references the same way options_format() does.
expand_env() {
	perl -e '
		my $s = $ARGV[0];
		while ($s =~ /\$(\w+)/) {
			my $v = $1;
			last unless defined $ENV{$v};
			$s =~ s/\$\Q$v\E/$ENV{$v}/g;
		}
		print $s;
	' "$1"
}

# First existing ancestor of a path, for writability checks on dirs that the
# pipeline will create itself.
nearest_existing() {
	local p="$1"
	while [ -n "$p" ] && [ "$p" != "/" ] && [ ! -d "$p" ]; do
		p=$( dirname "$p" )
	done
	printf '%s' "$p"
}

have() { command -v "$1" >/dev/null 2>&1; }

# ----------------------------------------------------------------- report ---

tap_init "LazypipeX Tier 0 — environment and prerequisites"
diag "host        : $( hostname )"
diag "date        : $( date -Is )"
diag "user        : ${USER:-?}"
diag "repo        : $REPO"
diag "install dir : $INSTALL"
diag "config      : $CONFIG"
if [ "$INSTALL" != "$REPO" ]; then
	diag "NOTE: \$LAZYPIPE_INSTALL_DIR differs from this checkout; testing the installed copy."
fi
if git -C "$REPO" rev-parse --short HEAD >/dev/null 2>&1; then
	diag "git         : $( git -C "$REPO" rev-parse --short HEAD ) on $( git -C "$REPO" rev-parse --abbrev-ref HEAD )"
fi
if have module || [ -n "${LOADEDMODULES:-}" ]; then
	diag "modules     : ${LOADEDMODULES:-none}"
fi
diag ""

# ================================================== ENV-01 env variables ====

# The taxonomy variable is whatever config.yaml actually references, not a
# hard-coded name — the two have drifted before (docs say $taxonomy, the config
# says $taxonomy_ncbi).
TAXVAR=$( cfg_get taxonomy db | sed -n 's/.*\$\([A-Za-z_][A-Za-z0-9_]*\).*/\1/p' )
RESVAR=$( cfg_get res         | sed -n 's/.*\$\([A-Za-z_][A-Za-z0-9_]*\).*/\1/p' )

env01_bad=""
for v in databases hostgenomes ${TAXVAR:-}; do
	[ -n "$v" ] || continue
	val="${!v:-}"
	if [ -z "$val" ]; then
		env01_bad="$env01_bad\n  \$$v is not set"
	elif [ ! -d "$val" ]; then
		env01_bad="$env01_bad\n  \$$v=$val is not a directory"
	else
		diag "\$$v = $val"
	fi
done
# $data (or whatever res: references) has no directory requirement yet, but it
# must be set: an unset variable leaves the literal string '$data/results' as
# the results root, which lands inside the current working directory.
if [ -n "${RESVAR:-}" ]; then
	if [ -z "${!RESVAR:-}" ]; then
		env01_bad="$env01_bad\n  \$$RESVAR is not set, but config res: '$( cfg_get res )' needs it"
	else
		diag "\$$RESVAR = ${!RESVAR}"
	fi
fi

if [ -z "$env01_bad" ]; then
	pass ENV-01 "environment variables set and pointing at existing directories"
else
	fail ENV-01 "environment variables set and pointing at existing directories" "$( printf '%b' "$env01_bad" )"
fi

# ==================================================== ENV-02 config parse ====

if [ ! -f "$CONFIG" ]; then
	fail ENV-02 "config.yaml parses" "no config file at $CONFIG"
else
	try perl -MYAML::Tiny -e 'YAML::Tiny->read($ARGV[0]) or die "unparseable\n"' "$CONFIG"
	if [ "$_RC" -eq 0 ]; then
		pass ENV-02 "config.yaml parses"
	else
		fail ENV-02 "config.yaml parses" "$_OUT"
	fi
fi

# ============================================= ENV-03 config env var scan ====

# Variables that are legitimately unset outside a specific context.  Everything
# else that config.yaml references must resolve, or the pipeline will warn and
# then use an unexpanded path.
CONDITIONAL_VARS="LOCAL_SCRATCH TM"

try_sh "perl -MYAML::Tiny -e '
	my \$y = YAML::Tiny->read(\$ARGV[0]) or die;
	my %vars;
	sub walk {
		my \$n = shift;
		if    (ref \$n eq \"HASH\")  { walk(\$_) for values %\$n }
		elsif (ref \$n eq \"ARRAY\") { walk(\$_) for @\$n }
		elsif (defined \$n)         { while (\$n =~ /\\\$(\\w+)/g) { \$vars{\$1} = 1 } }
	}
	walk(\$y->[0]);
	for my \$v (sort keys %vars) {
		print( (defined \$ENV{\$v} && length \$ENV{\$v}) ? \"SET\\t\$v\\n\" : \"UNSET\\t\$v\\n\" );
	}
' '$CONFIG'"

if [ "$_RC" -ne 0 ]; then
	fail ENV-03 "every \$var referenced by config.yaml is defined" "$_OUT"
else
	env03_bad=""
	env03_cond=""
	while IFS=$'\t' read -r state var; do
		[ -n "${var:-}" ] || continue
		[ "$state" = "UNSET" ] || continue
		case " $CONDITIONAL_VARS " in
			*" $var "*) env03_cond="$env03_cond $var" ;;
			*)          env03_bad="$env03_bad\n  \$$var referenced but not set" ;;
		esac
	done <<< "$_OUT"

	if [ -n "$env03_cond" ]; then
		diag "conditional vars unset (expected off-Slurm / without Trimmomatic):$env03_cond"
	fi
	if [ -z "$env03_bad" ]; then
		pass ENV-03 "every \$var referenced by config.yaml is defined"
	else
		fail ENV-03 "every \$var referenced by config.yaml is defined" \
			"$( printf '%b' "$env03_bad" )" \
			"lazypipe.pl would warn and then use the path unexpanded"
	fi
fi

# ==================================================== ENV-04 perl modules ====

PERL_MODULES="File::Basename File::Temp Getopt::Long MIME::Base64 YAML::Tiny Sort::Naturally POSIX Cwd"
env04_missing=""
for m in $PERL_MODULES; do
	if ! perl -M"$m" -e1 >/dev/null 2>&1; then
		env04_missing="$env04_missing $m"
	fi
done
if [ -z "$env04_missing" ]; then
	pass ENV-04 "required Perl modules importable"
else
	fail ENV-04 "required Perl modules importable" "missing:$env04_missing" \
		"install with: cpan$env04_missing"
fi

# ================================================= ENV-05 Lazypipe modules ====

try perl -I "$INSTALL/perl" -MLazypipe::Utils -MLazypipe::SeqAn -e1
if [ "$_RC" -ne 0 ]; then
	fail ENV-05 "Lazypipe::Utils and Lazypipe::SeqAn load" "$_OUT"
elif [ -n "$_OUT" ]; then
	fail ENV-05 "Lazypipe::Utils and Lazypipe::SeqAn load cleanly" \
		"loaded, but emitted warnings:" "$_OUT"
else
	pass ENV-05 "Lazypipe::Utils and Lazypipe::SeqAn load cleanly"
fi

# ========================================================= ENV-06 call_R ====

CALL_R=$( cfg_get general.parameters call_R )
R_WORKING=""          # the invocation later tests should use
if [ -z "$CALL_R" ]; then
	fail ENV-06 "R reachable through config call_R" "general.parameters.call_R is not set in $CONFIG"
else
	diag "call_R      : $CALL_R"
	try_sh "$CALL_R -e 'cat(R.version.string)'"
	if [ "$_RC" -eq 0 ] && printf '%s' "$_OUT" | grep -qi 'R version'; then
		R_WORKING="$CALL_R"
		pass ENV-06 "R reachable through config call_R"
		diag "R           : $( printf '%s' "$_OUT" | tail -1 )"
	else
		# Fall back to a plain Rscript so the report says whether R itself is
		# missing or only the wrapper named in config.yaml is.
		if have Rscript; then
			try_sh "Rscript --no-save -e 'cat(R.version.string)'"
			if [ "$_RC" -eq 0 ]; then
				R_WORKING="Rscript --no-save"
				fail ENV-06 "R reachable through config call_R" \
					"call_R = '$CALL_R' does not run here," \
					"but a plain 'Rscript' does: $( printf '%s' "$_OUT" | tail -1 )" \
					"set general.parameters.call_R to 'Rscript --no-save' for this site"
			else
				fail ENV-06 "R reachable through config call_R" \
					"call_R = '$CALL_R' failed, and plain Rscript also failed:" "$_OUT"
			fi
		else
			fail ENV-06 "R reachable through config call_R" \
				"call_R = '$CALL_R' failed and no Rscript on PATH:" "$_OUT"
		fi
	fi
fi

# ===================================================== ENV-07 R libraries ====

# Taken from library()/require() calls in R/*.R, not from the User Guide.
R_LIBRARIES="reshape openxlsx ggplot2 cowplot colorspace"
if [ -z "$R_WORKING" ]; then
	skipt ENV-07 "R libraries installed" "no working R invocation (see ENV-06)"
else
	env07_missing=""
	for lib in $R_LIBRARIES; do
		try_sh "$R_WORKING -e 'suppressMessages(library($lib))'"
		if [ "$_RC" -ne 0 ]; then
			env07_missing="$env07_missing $lib"
		fi
	done
	if [ -z "$env07_missing" ]; then
		pass ENV-07 "R libraries installed ($R_LIBRARIES)"
	else
		fail ENV-07 "R libraries installed" "missing:$env07_missing" \
			"install with: install.packages(c($( echo "$env07_missing" | sed 's/^ //; s/ /", "/g; s/^/"/; s/$/"/' )))"
	fi
fi

# ========================================================== ENV-08 tools ====

# Required for a default `-p main` run (fastp preprocessing, bwa filtering,
# megahit assembly, minimap annotation, mga ORFs, krona report).
TOOLS_REQUIRED="perl fastp bwa samtools seqkit csvtk taxonkit megahit minimap2 mga ktImportText pigz wget tar"
# Needed only when a non-default option is chosen.
TOOLS_OPTIONAL="spades.py:--ass spades  trimmomatic:--pre trimm  java:--pre trimm  prodigal:--gen prod  orfipy:ORF prediction alternative"
# Needed only by particular annotation strategies / pipeline steps.
TOOLS_FEATURE="blastn:BLASTN strategies  blastp:BLASTP strategies  blastx:BLASTX strategies  blastdbcmd:database checks (Tier 2)  diamond:vi.chain3 strategies  hmmscan:HMM strategies  runsanspanz.py:SANS strategies  create_report:--pipe rgrep (igv-reports)"

env08_missing=""
for t in $TOOLS_REQUIRED; do
	if have "$t"; then
		diag "$( printf '%-14s %s  [%s]' "$t" "$( command -v "$t" )" "$( tool_version "$t" )" )"
	else
		env08_missing="$env08_missing $t"
	fi
done
if [ -z "$env08_missing" ]; then
	pass ENV-08 "tools required for a default 'main' run are on PATH"
else
	fail ENV-08 "tools required for a default 'main' run are on PATH" \
		"missing:$env08_missing" \
		"a default '-p main' run will fail at the step that calls them"
fi

report_optional_group() {   # $1 = test id, $2 = description, $3 = "tool:reason  tool:reason"
	local id="$1" desc="$2" spec="$3" missing="" t reason
	local oldifs="$IFS"
	IFS='  '
	for entry in $spec; do
		[ -n "$entry" ] || continue
		t="${entry%%:*}"
		reason="${entry#*:}"
		if have "$t"; then
			diag "$( printf '%-14s %s  [%s]' "$t" "$( command -v "$t" )" "$( tool_version "$t" )" )"
		else
			missing="$missing $t($reason)"
		fi
	done
	IFS="$oldifs"
	if [ -z "$missing" ]; then
		pass "$id" "$desc"
	else
		skipt "$id" "$desc" "not installed:$missing"
	fi
	return 0
}

report_optional_group ENV-08a "option-gated tools available" "$TOOLS_OPTIONAL"
report_optional_group ENV-08b "strategy-gated tools available" "$TOOLS_FEATURE"

# ================================================ ENV-09 writable scratch ====

env09_bad=""

RES_RAW=$( cfg_get res )
RES=$( expand_env "$RES_RAW" )
if printf '%s' "$RES" | grep -q '\$'; then
	env09_bad="$env09_bad\n  results root '$RES_RAW' did not fully expand -> '$RES'"
else
	RES_BASE=$( nearest_existing "$RES" )
	if [ -z "$RES_BASE" ] || [ ! -w "$RES_BASE" ]; then
		env09_bad="$env09_bad\n  results root '$RES' is not writable (nearest existing dir: ${RES_BASE:-none})"
	else
		probe="$RES_BASE/.lazytest.$$"
		if ( : > "$probe" ) 2>/dev/null; then
			rm -f "$probe"
			diag "results root: $RES (writable via $RES_BASE)"
		else
			env09_bad="$env09_bad\n  cannot create files under '$RES_BASE'"
		fi
	fi
fi

TMP_RAW=$( cfg_get tmpdir )
if [ -z "$TMP_RAW" ]; then
	diag "tmpdir      : not set in config; File::Temp default will be used"
elif printf '%s' "$TMP_RAW" | grep -q 'LOCAL_SCRATCH' && [ -z "${SLURM_JOB_ID:-}" ]; then
	diag "tmpdir      : '$TMP_RAW' — \$LOCAL_SCRATCH is node-local and only set inside a Slurm job;"
	diag "              checked by HPC-03, not here."
else
	TMP=$( expand_env "$TMP_RAW" )
	if printf '%s' "$TMP" | grep -q '\$'; then
		env09_bad="$env09_bad\n  tmpdir '$TMP_RAW' did not fully expand -> '$TMP'"
	else
		TMP_BASE=$( nearest_existing "$TMP" )
		if [ -z "$TMP_BASE" ] || [ ! -w "$TMP_BASE" ]; then
			env09_bad="$env09_bad\n  tmpdir '$TMP' is not writable (nearest existing dir: ${TMP_BASE:-none})"
		else
			diag "tmpdir      : $TMP (writable via $TMP_BASE)"
		fi
	fi
fi

if [ -z "$env09_bad" ]; then
	pass ENV-09 "results and temporary directories are writable"
else
	fail ENV-09 "results and temporary directories are writable" "$( printf '%b' "$env09_bad" )"
fi

# =============================================== ENV-10 compiled helpers ====

# retrieve_reads is the one the documented workflow uses directly; the SeqAn
# based helpers are built separately and may legitimately be absent.
HELPERS_REQUIRED="retrieve_reads"
HELPERS_OPTIONAL="get_contigs filtfa filtfq"

check_helper() {   # $1 = name -> prints a diagnosis, returns 0 ok / 1 broken / 2 absent
	local bin="$INSTALL/bin/$1"
	if [ ! -e "$bin" ]; then
		HELPER_MSG="not built ($bin missing)"
		return 2
	fi
	if [ ! -x "$bin" ]; then
		HELPER_MSG="not executable: $bin"
		return 1
	fi
	# No arguments: each helper prints its usage block and exits 0 or 1.
	# 126/127 or a loader error means the binary cannot run here at all.
	try "$bin"
	if [ "$_RC" -ge 126 ] || printf '%s' "$_OUT" | grep -qi 'error while loading shared libraries\|cannot execute'; then
		HELPER_MSG="cannot execute (rc=$_RC): $( printf '%s' "$_OUT" | head -2 )"
		return 1
	fi
	if ! printf '%s' "$_OUT" | grep -qi 'usage'; then
		HELPER_MSG="ran (rc=$_RC) but printed no usage block: $( printf '%s' "$_OUT" | head -2 )"
		return 1
	fi
	HELPER_MSG="ok (rc=$_RC, usage printed)"
	return 0
}

env10_bad=""
for h in $HELPERS_REQUIRED; do
	check_helper "$h"
	case $? in
		0) diag "$( printf '%-16s %s' "bin/$h" "$HELPER_MSG" )" ;;
		*) env10_bad="$env10_bad\n  bin/$h: $HELPER_MSG" ;;
	esac
done
if [ -z "$env10_bad" ]; then
	pass ENV-10 "compiled helpers present and executable"
else
	fail ENV-10 "compiled helpers present and executable" "$( printf '%b' "$env10_bad" )" \
		"rebuild with: make retrieve_reads"
fi

env10_opt_missing=""
for h in $HELPERS_OPTIONAL; do
	check_helper "$h"
	case $? in
		0) diag "$( printf '%-16s %s' "bin/$h" "$HELPER_MSG" )" ;;
		*) env10_opt_missing="$env10_opt_missing $h" ;;
	esac
done
if [ -z "$env10_opt_missing" ]; then
	pass ENV-10a "SeqAn-based helpers present and executable"
else
	skipt ENV-10a "SeqAn-based helpers present and executable" \
		"not usable:$env10_opt_missing (built from cpp/ with \$seqan set; see Makefile)"
fi

# ------------------------------------------------------------------ done ---

diag ""
tap_done
