#!/usr/bin/env bash
#
# LazypipeX Tier 4 — end-to-end runs and biological validation.
# Implements E2E-01 … E2E-05 and BIO-01 … BIO-07 of docs/testing_roadmap.md §7.
#
# Tier 3 proves each step works.  This tier proves the assembled pipeline finds
# the right viruses: one full `main` run on the bundled mink faecal library,
# then assertions about what came out of it.
#
# The library (data/samples/M15small_R{1,2}.fastq, 9 842 pairs) is documented as
# carrying mink astrovirus and mink circovirus.  The BIO checks match on taxid
# *lineage* rather than species strings, because NCBI has already renamed both —
# "Mink circovirus" is now "Circovirus mink", and the roadmap records a taxid
# change (1475143 -> 3048202) that this installation's taxonomy has not yet made.
#
# Needs installed databases (Tier 2 green) and the full tool chain.  Wall time is
# roughly 15 min with the defaults; E2E-05 and BIO-07 add more and are opt-in.
#
# Usage:
#     module use /projappl/project_2003755/Lazypipe-db/modulefiles/projects
#     module load lazypipe/3.1
#     tests/t4_e2e_bio.sh                # TAP on stdout
#     T4_MATRIX=1 tests/t4_e2e_bio.sh    # add the strategy matrix (E2E-05)
#     T4_NEGCTRL=1 tests/t4_e2e_bio.sh   # add the negative control (BIO-07)
#
# Exit status = number of failed tests (0 = all good).

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
R2="$INSTALL/data/samples/M15small_R2.fastq"
INVENTORY="$TESTS_DIR/lib/db_inventory.pl"

: "${T4_ANNS:=vi.refseq}"			# the classical two-round minimap.vi -> blastn.abv
: "${T4_HOSTGEN:=Neovison_vison,Homo_sapiens}"
: "${T4_NUMTH:=8}"
: "${T4_RUN_TIMEOUT:=3600}"
: "${T4_MATRIX:=0}"				# E2E-05: one full run per strategy
: "${T4_NEGCTRL:=0}"				# BIO-07: an extra full run on shuffled reads
: "${T4_KEEP:=0}"

# BIO-04 floors, as percentages of the RAW input reads.  The denominator is fixed
# by the fixture so the floors survive changes to trimming or --hostgen, and each
# floor sits several-fold below what a healthy run yields (observed in brackets).
: "${T4_MIN_VIRAL_PC:=5}"			# total viral            [11.67 %]
: "${T4_MIN_ASTRO_PC:=2}"			# Mamastrovirus 10       [ 8.91 %]
: "${T4_MIN_CIRCO_PC:=0.5}"			# Mink circovirus        [ 2.75 %]
: "${T4_MAX_DENGUE_PC:=0.5}"			# false-positive guard   [ 0.015 %]

# Expected agents.  Taxids are recorded for reporting; the assertions resolve
# lineages to families, which is what survives a rename.
ASTRO_TAXID=1239574
# The mink circovirus occupies two taxids, and a hit against either is correct:
# 3048202 "Circovirus mink" is the species, and 1475143 "Mink circovirus" is a
# no-rank child of it.  They are summed, not treated as alternatives — a run may
# place reads under either, or split them between the two.
CIRCO_TAXIDS="3048202 1475143"
DENGUE_TAXID=11070
EXPECTED_FAMILIES="Astroviridae Circoviridae"

WORK=$( mktemp -d "${TMPDIR:-/tmp}/lazytest-t4.XXXXXX" ) || exit 99
if [ "$T4_KEEP" = "1" ]; then
	trap 'printf "# results kept at %s\n" "$WORK"' EXIT
else
	trap 'rm -rf "$WORK"' EXIT
fi

RES="$WORK/res"
TMPD="$WORK/tmp"
LOGS="$WORK/logs"
MAIN="$RES/M15main"				# the E2E-01 result dir every BIO check reads
mkdir -p "$RES" "$TMPD" "$LOGS"

cd "$INSTALL" || exit 99

have() { command -v "$1" >/dev/null 2>&1; }
TREE_BEFORE=$( git -C "$REPO" status --porcelain 2>/dev/null )

PASSED=""
mark_ok() { PASSED="$PASSED $1"; }
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

RUN_TIMES=""
lz_run() {	# label, args...
	local label="$1" t0
	shift
	t0=$( date +%s )
	TEST_TIMEOUT="$T4_RUN_TIMEOUT" try_sh "perl '$LZP' -1 '$R1' --res '$RES' \
		-t $T4_NUMTH --tmpdir '$TMPD' --logs '$LOGS' -v $*"
	_SECS=$(( $( date +%s ) - t0 ))
	RUN_TIMES="$RUN_TIMES  $label ${_SECS}s\n"
	return 0
}

# Index of a named column, or empty.  Assertions look columns up by name because
# the report tables have gained and lost columns between generations.
colidx() { head -1 "$1" | tr '\t' '\n' | grep -nx -- "$2" | cut -d: -f1; }

# Sum of $3 over rows where column $2 equals $4.
sum_where() {	# file valcol keycol key
	local f="$1" vc="$2" kc="$3" key="$4"
	awk -F'\t' -v vc="$vc" -v kc="$kc" -v k="$key" 'NR>1 && $kc==k { s+=$vc } END { printf "%d", s+0 }' "$f"
}

pc_of_input() { awk -v n="$1" -v t="$RAW_READS" 'BEGIN { printf "%.2f", (t>0)? n/t*100 : 0 }'; }
ge() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a >= b) }'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a <  b) }'; }

IN_PAIRS=$(( $( wc -l < "$R1" ) / 4 ))
RAW_READS=$(( IN_PAIRS * 2 ))

# ----------------------------------------------------------------- report ---

tap_init "LazypipeX Tier 4 — end-to-end and biological validation"
diag "host        : $( hostname )"
diag "date        : $( date -Is )"
diag "install dir : $INSTALL"
diag "results     : $RES"
diag "library     : $IN_PAIRS pairs = $RAW_READS reads (the BIO-04 denominator)"
diag "strategy    : $T4_ANNS"
diag "hostgen     : $T4_HOSTGEN"
diag ""

# ==================================================== E2E-01 main run ========

lz_run E2E-01 "-p main --hostgen '$T4_HOSTGEN' --anns '$T4_ANNS' -s M15main"
# Per-step cross-cutting checks (no ERROR:, History.log, undefined vars) are
# Tier 3's job; this tier asks only whether the whole run succeeded and produced
# the artifact set.
e01=""
[ "$_RC" -eq 0 ] || e01="$e01\n  exit status $_RC: $( printf '%s' "$_OUT" | grep -iE '^(ERROR|Fatal)' | head -1 | cut -c1-90 )"
for f in contigs.fa annot1.tsv abund_table.tsv abund_table.xlsx annot_table.tsv \
         annot_table.xlsx taxprofile.txt reports/krona.report.html; do
	[ -s "$MAIN/$f" ] || e01="$e01\n  missing or empty $f"
done
[ -d "$MAIN/contigs" ] || e01="$e01\n  missing contigs/"
npng=$( ls "$MAIN"/figures/*.png "$MAIN"/reports/figures/*.png 2>/dev/null | grep -c . )
[ "$npng" -gt 0 ] || e01="$e01\n  no QC plots written"
if [ -z "$e01" ]; then
	pass E2E-01 "main run on $T4_ANNS: full artifact set present ($npng plots)"
	mark_ok E2E-01
else
	fail E2E-01 "main run on $T4_ANNS" "$( printf '%b' "$e01" )"
fi

# ==================================================== E2E-02 all + rgrep =====

if ! have create_report; then
	skipt E2E-02 "'all' run adds the reference-genome report" "create_report (igv-reports) not on PATH"
elif need E2E-02 "'all' run adds the reference-genome report" E2E-01; then
	lz_run E2E-02 "-p all --hostgen '$T4_HOSTGEN' --anns '$T4_ANNS' -s M15all"
	e02=""
	[ "$_RC" -eq 0 ] || e02="$e02\n  exit status $_RC"
	[ -s "$RES/M15all/reports/refgen.report.html" ] || e02="$e02\n  missing reports/refgen.report.html"
	if [ -z "$e02" ]; then
		pass E2E-02 "'all' run adds the reference-genome report"
		mark_ok E2E-02
	else
		fail E2E-02 "'all' run adds the reference-genome report" "$( printf '%b' "$e02" )"
	fi
fi

# ===================================================== E2E-03 single end =====

lz_run E2E-03 "--se -p main --anns '$T4_ANNS' -s M15se"
e03=""
[ "$_RC" -eq 0 ] || e03="$e03\n  exit status $_RC: $( printf '%s' "$_OUT" | grep -iE '^(ERROR|Fatal)' | head -1 | cut -c1-90 )"
[ -s "$RES/M15se/contigs.fa" ] || e03="$e03\n  no contigs.fa from the single-end run"
[ -s "$RES/M15se/abund_table.tsv" ] || e03="$e03\n  no abund_table.tsv from the single-end run"
if [ -z "$e03" ]; then
	pass E2E-03 "single-end run completes without a read2 ($( grep -c '^>' "$RES/M15se/contigs.fa" ) contigs)"
	mark_ok E2E-03
else
	fail E2E-03 "single-end run completes without a read2" "$( printf '%b' "$e03" )"
fi

# ==================================================== E2E-04 gzipped input ===

GZDIR="$WORK/gzin"
mkdir -p "$GZDIR"
gzip -c "$R1" > "$GZDIR/M15gz_R1.fastq.gz"
[ -s "$R2" ] && gzip -c "$R2" > "$GZDIR/M15gz_R2.fastq.gz"
TEST_TIMEOUT="$T4_RUN_TIMEOUT" try_sh "perl '$LZP' -1 '$GZDIR/M15gz_R1.fastq.gz' --res '$RES' \
	-t $T4_NUMTH --tmpdir '$TMPD' --logs '$LOGS' -v -p pre,ass -s M15gz"
e04=""
[ "$_RC" -eq 0 ] || e04="$e04\n  exit status $_RC: $( printf '%s' "$_OUT" | grep -iE '^(ERROR|Fatal)' | head -1 | cut -c1-90 )"
if [ ! -s "$RES/M15gz/contigs.fa" ]; then
	e04="$e04\n  no contigs.fa from the gzipped run"
elif [ -s "$MAIN/contigs.fa" ]; then
	# Same reads in, same assembly out: gzip must be transparent.  Compared on
	# contig count rather than byte identity, since megahit is not bit-reproducible.
	ngz=$( grep -c '^>' "$RES/M15gz/contigs.fa" )
	nplain=$( grep -c '^>' "$MAIN/contigs.fa" )
	lo=$(( nplain * 8 / 10 )); hi=$(( nplain * 12 / 10 ))
	if [ "$ngz" -lt "$lo" ] || [ "$ngz" -gt "$hi" ]; then
		e04="$e04\n  gzipped run gave $ngz contigs vs $nplain uncompressed — not equivalent"
	fi
fi
if [ -z "$e04" ]; then
	pass E2E-04 "gzipped input gives an equivalent assembly"
	mark_ok E2E-04
else
	fail E2E-04 "gzipped input gives an equivalent assembly" "$( printf '%b' "$e04" )"
fi

# =================================================== E2E-05 strategy matrix ==

if [ "$T4_MATRIX" != "1" ]; then
	skipt E2E-05 "every installed strategy runs and annotates" \
		"set T4_MATRIX=1 (one full run per strategy; the nt-based ones index 200+ GB FASTAs)"
elif [ ! -f "$INVENTORY" ]; then
	skipt E2E-05 "every installed strategy runs and annotates" "missing $INVENTORY"
else
	perl "$INVENTORY" "$INSTALL/config.yaml" > "$WORK/inv.tsv" 2>/dev/null
	# A strategy is runnable only if every database it names is installed.
	perl -MYAML::Tiny -e '
		my $y = YAML::Tiny->read($ARGV[0]) or exit 2;
		my %S = %{ $y->[0]{"ann.strategies"} // {} };
		my %inst;
		open my $i, "<", $ARGV[1] or exit 2;
		while(<$i>){ chomp; my @f = split /\t/; $inst{$f[1]} = $f[4] if @f >= 5 }
		for my $s ( sort keys %S ){
			my @need;
			my $v = $S{$s};
			while( $v =~ /--ann[12]\s+(\S+)/g ){
				for my $r ( split /,/, $1 ){ $r =~ s/^\w+://; push @need, $r if $r ne "" }
			}
			my @missing = grep { !$inst{$_} } @need;
			print @missing ? "SKIP\t$s\t@missing\n" : "RUN\t$s\n";
		}
	' "$INSTALL/config.yaml" "$WORK/inv.tsv" > "$WORK/strategies.txt" 2>/dev/null

	e05_bad=""; e05_ran=0
	while IFS=$'\t' read -r verdict strat rest; do
		[ "$verdict" = "RUN" ] || { diag "  E2E-05 skip $strat (databases not installed: ${rest:-?})"; continue; }
		e05_ran=$(( e05_ran + 1 ))
		lz_run "E2E-05/$strat" "-p ann1,ann2 --anns '$strat' -s M15main"
		if [ "$_RC" -ne 0 ]; then
			e05_bad="$e05_bad\n  $strat: exit $_RC"
		elif [ ! -s "$MAIN/annot1.tsv" ]; then
			e05_bad="$e05_bad\n  $strat: annot1.tsv empty"
		fi
	done < "$WORK/strategies.txt"
	if [ "$e05_ran" -eq 0 ]; then
		skipt E2E-05 "every installed strategy runs and annotates" "no strategy has all its databases installed"
	elif [ -z "$e05_bad" ]; then
		pass E2E-05 "every installed strategy runs and annotates ($e05_ran run)"
	else
		fail E2E-05 "every installed strategy runs and annotates" "$( printf '%b' "$e05_bad" )"
	fi
fi

# ############################################################################ #
#                       §7.2 — biological validation                           #
# ############################################################################ #

AB="$MAIN/abund_table.tsv"
AN="$MAIN/annot_table.tsv"

# Resolve every taxid in the abundance table to its family, once, and reuse it.
# Matching on family rather than species string is what survives the NCBI
# renames that already invalidated "Mink circovirus".
FAMMAP="$WORK/taxid_family.tsv"
: > "$FAMMAP"
if [ -s "$AB" ] && have taxonkit && [ -n "${taxonomy_ncbi:-}" ]; then
	tcol=$( colidx "$AB" taxid )
	if [ -n "$tcol" ]; then
		awk -F'\t' -v c="$tcol" 'NR>1 && $c ~ /^[0-9]+$/ { print $c }' "$AB" | sort -u \
			| taxonkit reformat --data-dir "$taxonomy_ncbi" -I 1 -f '{f}' -r NA -R NA 2>/dev/null \
			> "$FAMMAP"
	fi
fi
family_of() { awk -F'\t' -v t="$1" '$1==t { print $2; exit }' "$FAMMAP"; }

# ================================================ BIO-01 expected taxa =======

if need BIO-01 "both expected agents are detected" E2E-01; then
	b01=""
	if [ ! -s "$AB" ]; then
		b01="$b01\n  no abund_table.tsv"
	elif [ ! -s "$FAMMAP" ]; then
		b01="$b01\n  could not resolve taxids to families (taxonkit or \$taxonomy_ncbi missing)"
	else
		for fam in $EXPECTED_FAMILIES; do
			if ! awk -F'\t' -v f="$fam" '$2==f { found=1 } END { exit !found }' "$FAMMAP"; then
				b01="$b01\n  no taxon in abund_table resolves to $fam"
			fi
		done
		# Report which of the accepted taxids this run actually used.
		tcol=$( colidx "$AB" taxid )
		for t in $ASTRO_TAXID $CIRCO_TAXIDS; do
			awk -F'\t' -v c="$tcol" -v t="$t" 'NR>1 && $c==t { f=1 } END { exit !f }' "$AB" \
				&& diag "  taxid $t present ($( family_of "$t" ))"
		done
	fi
	if [ -z "$b01" ]; then
		pass BIO-01 "both expected agents are detected (by family lineage)"
		mark_ok BIO-01
	else
		fail BIO-01 "both expected agents are detected" "$( printf '%b' "$b01" )" \
			"the library is documented as carrying mink astrovirus and mink circovirus"
	fi
fi

# ============================================= BIO-02 family assignment ======

if need BIO-02 "expected families appear in annot_table" E2E-01; then
	b02=""
	fcol=$( colidx "$AN" family )
	if [ -z "$fcol" ]; then
		b02="$b02\n  annot_table.tsv has no family column"
	else
		for fam in $EXPECTED_FAMILIES; do
			awk -F'\t' -v c="$fcol" -v f="$fam" 'NR>1 && $c==f { n=1 } END { exit !n }' "$AN" \
				|| b02="$b02\n  $fam absent from the family column"
		done
	fi
	if [ -z "$b02" ]; then
		pass BIO-02 "expected families appear in annot_table"
		mark_ok BIO-02
	else
		fail BIO-02 "expected families appear in annot_table" "$( printf '%b' "$b02" )"
	fi
fi

# ============================================== BIO-03 contig support ========

if need BIO-03 "each expected family has contig-level support" E2E-01; then
	b03=""
	fcol=$( colidx "$AN" family ); qcol=$( colidx "$AN" qcov )
	ccol=$( colidx "$AN" contig ); bcol=$( colidx "$AN" bitscore )
	minq=$( perl -MYAML::Tiny -e '
		my $y = YAML::Tiny->read("config.yaml");
		print $y->[0]{"general.parameters"}{min_qcov_annot} // 0.20;' 2>/dev/null )
	minq=${minq:-0.20}
	if [ -z "$fcol" ] || [ -z "$qcol" ] || [ -z "$ccol" ]; then
		b03="$b03\n  annot_table.tsv is missing family/qcov/contig columns"
	else
		for fam in $EXPECTED_FAMILIES; do
			n=$( awk -F'\t' -v fc="$fcol" -v qc="$qcol" -v cc="$ccol" -v f="$fam" -v q="$minq" \
				'NR>1 && $fc==f && $qc+0 >= q+0 { c[$cc]=1 } END { print length(c) }' "$AN" )
			if [ "${n:-0}" -lt 1 ]; then
				b03="$b03\n  $fam has no contig with qcov >= $minq"
			else
				bs=$( awk -F'\t' -v fc="$fcol" -v bc="$bcol" -v f="$fam" \
					'NR>1 && $fc==f { if(m=="" || $bc+0<m) m=$bc+0 } END { print m }' "$AN" )
				diag "  $fam: $n contig(s) at qcov >= $minq, lowest bitscore $bs"
			fi
		done
	fi
	# The per-engine bitscore minimums (min_blastn_bits, min_minimap_DPpeak_score,
	# …) are applied by the pipeline before a row reaches annot_table, so
	# re-asserting them here would only re-check the pipeline's own filter.  The
	# observed minimum is reported above instead.
	if [ -z "$b03" ]; then
		pass BIO-03 "each expected family has contig-level support"
		mark_ok BIO-03
	else
		fail BIO-03 "each expected family has contig-level support" "$( printf '%b' "$b03" )"
	fi
fi

# ============================================== BIO-04 read fractions ========

if need BIO-04 "expected taxa clear their read-fraction floors" E2E-01; then
	b04=""
	rcol=$( colidx "$AB" readn ); tcol=$( colidx "$AB" taxid ); dcol=$( colidx "$AB" division )
	if [ -z "$rcol" ] || [ -z "$tcol" ]; then
		b04="$b04\n  abund_table.tsv is missing readn/taxid columns"
	else
		viral=$( awk -F'\t' -v r="$rcol" -v d="$dcol" 'NR>1 && $d=="Viruses" { s+=$r } END { printf "%d", s+0 }' "$AB" )
		astro=$( sum_where "$AB" "$rcol" "$tcol" "$ASTRO_TAXID" )
		circo=0
		for t in $CIRCO_TAXIDS; do
			circo=$(( circo + $( sum_where "$AB" "$rcol" "$tcol" "$t" ) ))
		done
		dengue=$( sum_where "$AB" "$rcol" "$tcol" "$DENGUE_TAXID" )

		vpc=$( pc_of_input "$viral" ); apc=$( pc_of_input "$astro" )
		cpc=$( pc_of_input "$circo" ); dpc=$( pc_of_input "$dengue" )
		diag "  of $RAW_READS raw reads: viral ${vpc}%, astro ${apc}%, circo ${cpc}%, dengue ${dpc}%"

		ge "$vpc" "$T4_MIN_VIRAL_PC"  || b04="$b04\n  total viral ${vpc}% is below the ${T4_MIN_VIRAL_PC}% floor"
		ge "$apc" "$T4_MIN_ASTRO_PC"  || b04="$b04\n  Mamastrovirus ${apc}% is below the ${T4_MIN_ASTRO_PC}% floor"
		ge "$cpc" "$T4_MIN_CIRCO_PC"  || b04="$b04\n  Mink circovirus ${cpc}% is below the ${T4_MIN_CIRCO_PC}% floor"
		lt "$dpc" "$T4_MAX_DENGUE_PC" || b04="$b04\n  dengue ${dpc}% exceeds the ${T4_MAX_DENGUE_PC}% false-positive guard"
	fi
	if [ -z "$b04" ]; then
		pass BIO-04 "expected taxa clear their read-fraction floors"
		mark_ok BIO-04
	else
		fail BIO-04 "expected taxa clear their read-fraction floors" "$( printf '%b' "$b04" )"
	fi
fi

# ============================================ BIO-05 retrieval round-trip ====

RR="$INSTALL/bin/retrieve_reads"
if [ ! -x "$RR" ]; then
	skipt BIO-05 "read retrieval round-trip" "bin/retrieve_reads not built"
elif need BIO-05 "read retrieval round-trip" E2E-01; then
	b05=""
	# retrieve_reads reads only uncompressed fastq (§12 item 15).
	gunzip -kf "$MAIN"/reads/read1.trim.fq.gz "$MAIN"/reads/read2.trim.fq.gz 2>/dev/null
	TEST_TIMEOUT="$T4_RUN_TIMEOUT" try_sh "'$RR' -r '$MAIN' -t $ASTRO_TAXID -p bio05"
	[ "$_RC" -eq 0 ] || b05="$b05\n  retrieve_reads -t $ASTRO_TAXID exited $_RC"
	got="$MAIN/reads/bio05_r1.fq"
	if [ ! -s "$got" ]; then
		b05="$b05\n  no FASTQ written for taxid $ASTRO_TAXID"
	else
		nret=$(( $( wc -l < "$got" ) / 4 ))
		rcol=$( colidx "$AB" readn ); tcol=$( colidx "$AB" taxid )
		nexp=$( sum_where "$AB" "$rcol" "$tcol" "$ASTRO_TAXID" )
		diag "  retrieved $nret reads, abund_table reports readn=$nexp"
		[ "$nret" -gt 0 ] || b05="$b05\n  retrieved 0 reads"
		if [ "${nexp:-0}" -gt 0 ] && [ "$nret" -ne "$nexp" ]; then
			b05="$b05\n  retrieved $nret reads but abund_table reports readn=$nexp"
		fi
	fi
	if [ -z "$b05" ]; then
		pass BIO-05 "read retrieval round-trip matches the reported readn"
		mark_ok BIO-05
	else
		fail BIO-05 "read retrieval round-trip" "$( printf '%b' "$b05" )"
	fi
fi

# ============================================ BIO-06 false-positive sweep ====

# §7 says violations are reviewed, not auto-failed, so this reports and passes.
if need BIO-06 "no unexpected viral family dominates" E2E-01; then
	rcol=$( colidx "$AB" readn ); fcol=$( colidx "$AB" family )
	dcol=$( colidx "$AB" division ); pcol=$( colidx "$AB" bphage )
	flagged=""
	if [ -n "$rcol" ] && [ -n "$fcol" ]; then
		while IFS=$'\t' read -r fam rn; do
			[ -n "$fam" ] && [ "$fam" != "NA" ] || continue
			case " $EXPECTED_FAMILIES " in *" $fam "*) continue ;; esac
			p=$( pc_of_input "$rn" )
			ge "$p" 1 && flagged="$flagged\n    $fam ${p}% ($rn reads)"
		done <<EOF
$( awk -F'\t' -v r="$rcol" -v f="$fcol" -v d="$dcol" 'NR>1 && $d=="Viruses" { s[$f]+=$r } END { for(k in s) print k"\t"s[k] }' "$AB" )
EOF
	fi
	if [ -n "$flagged" ]; then
		pass BIO-06 "no unexpected viral family dominates (review the list below)"
		diag "  viral families above 1% that were not expected:"
		diag "$( printf '%b' "$flagged" )"
		diag "  §7 says these are reviewed by hand, not auto-failed"
	else
		pass BIO-06 "no unexpected viral family dominates"
	fi
	mark_ok BIO-06
fi

# =============================================== BIO-07 negative control =====

if [ "$T4_NEGCTRL" != "1" ]; then
	skipt BIO-07 "negative control finds no viruses" "set T4_NEGCTRL=1 (adds a full pipeline run)"
elif ! need BIO-07 "negative control finds no viruses" E2E-01; then
	:
else
	# Shuffle the bases within each read: composition and quality are preserved,
	# homology is destroyed.  Fixed seed, so the control is reproducible.
	NEG="$WORK/neg_R1.fastq"
	perl -e '
		srand(42);
		my $i = 0;
		while( my $h = <> ){
			my $s = <>; my $p = <>; my $q = <>;
			chomp $s;
			my @b = split //, $s;
			for( my $j = @b; --$j; ){ my $k = int rand($j+1); @b[$j,$k] = @b[$k,$j] }
			print $h, join("",@b), "\n", $p, $q;
		}
	' "$R1" > "$NEG"
	lz_run BIO-07 "-1 '$NEG' -p pre,ass,rea,ann1 --anns '$T4_ANNS' -s M15neg"
	b07=""
	[ "$_RC" -eq 0 ] || b07="$b07\n  the control run itself failed (exit $_RC)"
	NEGA="$RES/M15neg/annot1.tsv"
	if [ -s "$NEGA" ]; then
		dcol=$( colidx "$NEGA" division )
		nvi=$( awk -F'\t' -v c="$dcol" 'NR>1 && $c=="Viruses" { n++ } END { print n+0 }' "$NEGA" )
		diag "  shuffled control produced $nvi viral annotation row(s)"
		[ "${nvi:-0}" -eq 0 ] || b07="$b07\n  $nvi viral hit(s) on shuffled reads — the index reports hits for anything"
	else
		diag "  shuffled control produced no annot1.tsv at all, which is the expected outcome"
	fi
	if [ -z "$b07" ]; then
		pass BIO-07 "negative control finds no viruses"
		mark_ok BIO-07
	else
		fail BIO-07 "negative control finds no viruses" "$( printf '%b' "$b07" )"
	fi
fi

# ------------------------------------------------------------------ done ---

if [ -n "$RUN_TIMES" ]; then
	diag "run wall time (recorded, never asserted on):"
	diag "$( printf '%b' "$RUN_TIMES" )"
fi

if [ -z "$TREE_BEFORE" ]; then
	skipt BIO-99 "the tier wrote nothing into the repository working tree" "not a git checkout"
else
	TREE_AFTER=$( git -C "$REPO" status --porcelain 2>/dev/null )
	if [ "$TREE_BEFORE" = "$TREE_AFTER" ]; then
		pass BIO-99 "the tier wrote nothing into the repository working tree"
	else
		fail BIO-99 "the tier wrote nothing into the repository working tree" \
			"$( diff <( printf '%s\n' "$TREE_BEFORE" ) <( printf '%s\n' "$TREE_AFTER" ) | head -10 )"
	fi
fi

diag ""
tap_done
