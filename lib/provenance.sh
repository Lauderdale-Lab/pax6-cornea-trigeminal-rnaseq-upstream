# shellcheck shell=bash
#
# lib/provenance.sh -- the provenance record written beside every count matrix.
#
# Sourced by slurm/count.sbatch after featureCounts succeeds. Writes
# counts/<study>/provenance/, which holds everything needed to state in a
# Methods section what was done, without reading a single log:
#
#   software_versions.tsv  every tool, the module loaded, the version it reported
#   reference.tsv          genome, annotation, index and splice sites, with MD5s
#   inputs.tsv             every library's raw FASTQ files and their MD5s
#   library_qc.tsv         per library: reads, trimming, alignment, assignment
#   strand.tsv             the strandedness measured for each dataset
#   parameters.env, modules.tsv, RUN_INFO.tsv, jobs.tsv   copies from the run
#   METHODS.md             a Methods paragraph filled in from all of the above
#
# The goal is that no fact in the Methods has to be reconstructed later.

md5_of() { md5sum "$1" 2>/dev/null | cut -d' ' -f1; }

write_provenance() {
  local study="$1" out="$2" P ds f
  P="$out/provenance"; mkdir -p "$P"
  cp "$RUN_DIR/parameters.env" "$RUN_DIR/modules.tsv" "$RUN_DIR/RUN_INFO.tsv" "$P/"
  [[ -f "$RUN_DIR/jobs.tsv" ]] && cp "$RUN_DIR/jobs.tsv" "$P/"
  mapfile -t DSS < <(study_datasets "$study")

  # --- software versions: gathered from every job of this study's datasets --
  {
    printf 'tool\tmodule\treported_version\tsteps\n'
    for ds in "${DSS[@]}" "$study"; do
      cat "$RUN_DIR"/provenance/versions/*."$ds".*.tsv 2>/dev/null || true
    done | awk -F'\t' '{ k = $1 "\t" $2 "\t" $3; if (!(k in s)) { s[k] = $4; o[++n] = k } else if (index(s[k], $4) == 0) s[k] = s[k] "," $4 }
                       END { for (i = 1; i <= n; i++) print o[i] "\t" s[o[i]] }' | sort
  } > "$P/software_versions.tsv"
  local multi
  multi=$(tail -n +2 "$P/software_versions.tsv" | cut -f1 | sort | uniq -d)
  [[ -z "$multi" ]] || warn "More than one version reported for: $(tr '\n' ' ' <<<"$multi") -- see $P/software_versions.tsv"

  # --- reference ---------------------------------------------------------------
  log "provenance: checksumming reference files"
  {
    printf 'item\tpath\tmd5\tdetail\n'
    printf 'genome_fasta\t%s\t%s\t%s sequences\n' "$GENOME_FASTA" "$(md5_of "$GENOME_FASTA")" "$(grep -c '^>' "$GENOME_FASTA")"
    printf 'annotation\t%s\t%s\t%s genes\n' "$ANNOTATION" "$(md5_of "$ANNOTATION")" \
      "$(awk -F'\t' '$3 == "gene"' "$ANNOTATION" | wc -l)"
    printf 'splice_sites\t%s\t%s\t%s junctions\n' "$RUN_DIR/reference/splice_sites.txt" \
      "$(md5_of "$RUN_DIR/reference/splice_sites.txt")" "$(wc -l < "$RUN_DIR/reference/splice_sites.txt" 2>/dev/null || echo NA)"
    for f in "$HISAT2_INDEX".*.ht2; do
      [[ -e "$f" ]] && printf 'index_file\t%s\t%s\t%s bytes\n' "$f" "$(md5_of "$f")" "$(stat -c %s "$f")"
    done
    printf 'index_type\t%s\t-\t%s\n' "$HISAT2_INDEX" "$INDEX_TYPE"
  } > "$P/reference.tsv"

  # --- raw inputs --------------------------------------------------------------
  {
    printf 'sample_id\tdataset\tbatch\tr1\tr1_md5\tr2\tr2_md5\n'
    for ds in "${DSS[@]}"; do
      while IFS=$'\t' read -r s r1 r2; do
        local m m1 m2; m="$(ds_dir "$ds")/checksums/$s.raw.md5"
        m1=$(awk -v f="$r1" '$2 == f { print $1 }' "$m" 2>/dev/null || true)
        m2=-; if [[ "$r2" != - ]]; then m2=$(awk -v f="$r2" '$2 == f { print $1 }' "$m" 2>/dev/null || true); fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$ds" "$(dataset_field "$ds" batch)" "$r1" "${m1:-NA}" "$r2" "${m2:-NA}"
      done < "$(samples_tsv "$ds")"
    done
  } > "$P/inputs.tsv"

  # --- per-library QC ------------------------------------------------------------
  {
    printf 'sample_id\tdataset\tbatch\tinput_reads_or_pairs\tsurviving_pct\tadapters\toverall_alignment_pct\tunique_alignment_pct\tassigned_fragments\tassigned_pct\tassigned_ratio_pct\n'
    for ds in "${DSS[@]}"; do
      while IFS=$'\t' read -r s _ _; do
        local tl hl inp surv ad oar uar aq
        tl="$(qc_dir "$ds")/trimmomatic/$s.trimmomatic.log"; hl="$(align_dir "$ds")/$s.hisat2.txt"
        inp=$(grep -oE '^Input (Read Pairs|Reads): [0-9]+' "$tl" 2>/dev/null | grep -oE '[0-9]+$' || true)
        surv=$(grep -oE '(Both Surviving|Surviving): [0-9]+ \([0-9.]+%\)' "$tl" 2>/dev/null | head -n1 | grep -oE '[0-9.]+%' | tr -d % || true)
        ad=$(cut -f2 "$(qc_dir "$ds")/adapters/${s}_adapter_source.tsv" 2>/dev/null || true)
        oar=$(grep -i 'overall alignment rate' "$hl" 2>/dev/null | grep -oE '[0-9.]+%' | tr -d % || true)
        uar=$(grep -E 'Aligned (concordantly )?1 time' "$hl" 2>/dev/null | head -n1 | grep -oE '[0-9.]+%' | tr -d % || true)
        aq=$(awk -F'\t' -v s="$s" '$1 == s { print $2 "\t" $4 }' "$out/assignment_qc.tsv" || true)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$ds" "$(dataset_field "$ds" batch)" \
          "${inp:-NA}" "${surv:-NA}" "${ad:-NA}" "${oar:-NA}" "${uar:-NA}" \
          "$(awk -F'\t' -v s="$s" 'NR == 1 { for (i = 2; i <= NF; i++) { n = split($i, p, "/"); c = p[n]; sub(/\.bam$/, "", c); if (c == s) col = i } }
                                    $1 == "Assigned" && col { print $col }' "$out/featureCounts.txt.summary")" \
          "${aq:-NA	NA}"
      done < "$(samples_tsv "$ds")"
    done
  } > "$P/library_qc.tsv"

  # --- strandedness --------------------------------------------------------------
  {
    printf 'dataset\t%s\n' "$(head -n1 "$(strand_dir "${DSS[0]}")/RESULT.tsv")"
    for ds in "${DSS[@]}"; do printf '%s\t%s\n' "$ds" "$(sed -n 2p "$(strand_dir "$ds")/RESULT.tsv")"; done
  } > "$P/strand.tsv"

  write_methods "$study" "$out" "$P"
}

# ver <tool> : "<Tool> <version>" for prose, from the version the tool reported
# (several versions are joined with "/", which the warning above flags);
# the module name, marked, if no version was recorded.
ver() {
  local v name
  v=$(awk -F'\t' -v t="$1" 'NR > 1 && $1 == t { print $3 }' "$P/software_versions.tsv" |
        grep -oE '[0-9]+(\.[0-9]+)+' | sort -u | paste -sd/ - || true)
  case "$1" in Subread) name="featureCounts (Subread %s)" ;; *) name="$1 %s" ;; esac
  if [[ -n "$v" ]]; then
    # shellcheck disable=SC2059  # the format string is chosen above
    printf "$name\n" "$v"
  else
    printf '%s (version not recorded)\n' "$(tsv_get "$RUN_DIR/modules.tsv" "$1" module)"
  fi
}

# pct <number> : one decimal place, or NA
pct() { [[ "$1" =~ ^[0-9.]+$ ]] && printf '%.1f' "$1" || printf NA; }

write_methods() {
  local study="$1" out="$2" P="$3" n ndet idx_desc aln asg rat surv strands
  n=$(( $(wc -l < "$P/library_qc.tsv") - 1 ))
  ndet=$(awk -F'\t' 'NR > 1 && $6 == "detected"' "$P/library_qc.tsv" | wc -l)
  read -r -a aln  < <(stats "$P/library_qc.tsv" overall_alignment_pct)
  read -r -a surv < <(stats "$P/library_qc.tsv" surviving_pct)
  read -r -a asg  < <(stats "$P/library_qc.tsv" assigned_pct)
  read -r -a rat  < <(stats "$P/library_qc.tsv" assigned_ratio_pct)
  case "$INDEX_TYPE" in
    plain)    idx_desc="a HISAT2 index without SNPs or built-in splice sites" ;;
    ss_exon)  idx_desc="a transcript-aware HISAT2 index (splice sites and exons from the annotation built in; no SNPs)" ;;
    snp_tran) idx_desc="a SNP- and transcript-aware HISAT2 index" ;;
    *)        idx_desc="a HISAT2 index of type $INDEX_TYPE" ;;
  esac
  local s_word
  case "$(awk -F'\t' '$1 == "strand" { print $2 }' "$out/PROVENANCE.tsv")" in
    0) s_word="unstranded (featureCounts -s 0)" ;;
    1) s_word="forward-stranded (featureCounts -s 1)" ;;
    2) s_word="reverse-stranded (featureCounts -s 2)" ;;
  esac
  strands=""
  local ds r a
  for ds in $(study_datasets "$study"); do
    r="$(strand_dir "$ds")/RESULT.tsv"
    a=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
          NR == 2 && c["assigned_s0"] { printf "%.1f%% / %.1f%% / %.1f%%", 100*$c["assigned_s0"], 100*$c["assigned_s1"], 100*$c["assigned_s2"] }' "$r" 2>/dev/null || true)
    [[ -n "$a" ]] && strands+="${strands:+; }${ds}: ${a}"
  done
  [[ -n "$strands" ]] && strands=" (fragments assigned at -s 0 / 1 / 2 — ${strands})"

  cat > "$P/METHODS.md" <<EOF
# Methods draft — ${study}, run ${RUN}

Generated $(date -Iseconds) from the recorded provenance of this matrix. Every
number and version below is taken from files in this folder; check the
wording, not the facts.

## RNA-seq processing

Raw reads from ${n} libraries were assessed with $(ver FastQC) before and after
trimming (both mates). Adapters and low-quality bases were removed with
$(ver Trimmomatic) (ILLUMINACLIP:<adapters>:${TRIM_ILLUMINACLIP} ${TRIM_STEPS}).
Adapter sequences were taken from each library's FastQC overrepresented
sequences where present (${ndet} of ${n} libraries) and otherwise from the
Trimmomatic TruSeq3 set. $(pct "${surv[0]}")–$(pct "${surv[1]}")% of reads or read pairs
survived trimming (median $(pct "${surv[2]}")%).

Trimmed reads were aligned to the ${ASSEMBLY} primary assembly
($(basename "$GENOME_FASTA")) with $(ver HISAT2) (${HISAT2_ARGS}), using
${idx_desc} and supplying $(awk -F'\t' '$1 == "splice_sites" { print $4 + 0 }' "$P/reference.tsv")
known splice junctions extracted from $(basename "$ANNOTATION") with
--known-splicesite-infile. Alignments were sorted and indexed with
$(ver SAMtools). Overall alignment rates were $(pct "${aln[0]}")–$(pct "${aln[1]}")% (median $(pct "${aln[2]}")%).

Library strandedness was measured, not assumed, by counting a subset of
libraries at each featureCounts strand setting${strands}. All libraries in
this study were ${s_word}.

Fragments were counted per gene with $(ver Subread)
(-s $(awk -F'\t' '$1 == "strand" { print $2 }' "$out/PROVENANCE.tsv") ${FC_PE_ARGS} ${FC_ARGS}) against
$(basename "$ANNOTATION"), all ${n} libraries in a single invocation.
$(pct "${asg[0]}")–$(pct "${asg[1]}")% of fragments were assigned to genes (median $(pct "${asg[2]}")%);
Assigned/(Assigned + Unassigned_NoFeatures) was $(pct "${rat[0]}")–$(pct "${rat[1]}")% (median $(pct "${rat[2]}")%).
QC reports were compiled with $(ver MultiQC).

Processing used pipeline commit $(run_info "$RUN" pipeline_commit | cut -c1-12)
(https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream).
Raw-file checksums: provenance/inputs.tsv. Reference checksums:
provenance/reference.tsv. Per-library values: provenance/library_qc.tsv.
EOF
}
