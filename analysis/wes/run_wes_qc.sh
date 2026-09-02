#!/usr/bin/env bash
# ==============================================================================
# WES QC and targeted-locus summaries for the VEXAS/MAS index case.
#
# Reads a BAM that is NEVER committed to this repository (see README.md) and
# writes only AGGREGATE summaries to data/wes/ and outputs/wes/.
#
# Usage:  ./run_wes_qc.sh /path/to/sample.bam
#
# Requires: samtools >= 1.17
# ==============================================================================

set -euo pipefail

BAM="${1:?usage: run_wes_qc.sh <bam>}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REF_DATA="$ROOT/analysis/wes"
OUT="$ROOT/outputs/wes"
DATA="$ROOT/data/wes"
mkdir -p "$OUT" "$DATA"

# hg19 / GRCh37 coordinates, resolved from Ensembl GRCh37 REST and checked
# against the reference sequence: codon 41 reads ATG and c.122 is T.
UBA1_REGION="chrX:47050260-47074527"
UBA1_C122="chrX:47058451"

echo "== integrity =="
samtools quickcheck -v "$BAM" || { echo "BAM failed quickcheck"; exit 1; }
[ -f "${BAM}.bai" ] || samtools index -@ 4 "$BAM"

echo "== reference build (from header) =="
samtools view -H "$BAM" | awk '$1=="@SQ"' | head -1
samtools view -H "$BAM" | awk '$1=="@SQ"' | wc -l | xargs echo "contigs:"

echo "== global counts =="
samtools flagstat -@ 4 "$BAM" > "$OUT/flagstat.txt"
samtools idxstats "$BAM" > "$OUT/idxstats.txt"

{
  echo "metric,value"
  awk '/in total/    {print "total_reads,"$1}'      "$OUT/flagstat.txt"
  awk '/ mapped \(/  {print "mapped_reads,"$1}'     "$OUT/flagstat.txt"
  awk '/duplicates/  {print "duplicate_reads,"$1; exit}' "$OUT/flagstat.txt"
} > "$DATA/qc_summary.csv"

echo "== UBA1 locus =="
# Allele counts at the p.Met41Thr position across filter settings.
{
  echo "base_quality_min,mapping_quality_min,depth,T_ref,C_alt,A,G,deletion"
  for Q in 0 13 20; do for q in 0 20; do
    L=$(samtools mpileup -r "$UBA1_C122-${UBA1_C122##*:}" -d 1000000 -Q "$Q" -q "$q" \
          --no-output-ins --no-output-del "$BAM" 2>/dev/null || true)
    D=$(echo "$L" | awk '{print $4}'); B=$(echo "$L" | awk '{print toupper($5)}')
    printf "%s,%s,%s,%s,%s,%s,%s,%s\n" "$Q" "$q" "${D:-0}" \
      "$(echo "$B" | tr -cd 'T' | wc -c | tr -d ' ')" \
      "$(echo "$B" | tr -cd 'C' | wc -c | tr -d ' ')" \
      "$(echo "$B" | tr -cd 'A' | wc -c | tr -d ' ')" \
      "$(echo "$B" | tr -cd 'G' | wc -c | tr -d ' ')" \
      "$(echo "$B" | tr -cd '*' | wc -c | tr -d ' ')"
  done; done
} > "$DATA/uba1_m41t_allele_counts.csv"

echo "== panel callability =="
# Exonic coverage per gene, over the interrogated panel.
BED="$REF_DATA/hlh_panel_exons.bed"
if [ -f "$BED" ]; then
  sort -k1,1 -k2,2n "$BED" > "$OUT/panel.sorted.bed"
  # single pass over the BAM, then attribute depths to genes by interval
  samtools depth -a -b "$OUT/panel.sorted.bed" -Q 13 -q 20 "$BAM" 2>/dev/null \
    > "$OUT/panel_depth.tsv"
  awk 'NR==FNR {n=split($0,f,"\t"); lo[FNR]=f[2]; hi[FNR]=f[3]; ch[FNR]=f[1]; gn[FNR]=f[4]; N=FNR; next}
       { for(i=1;i<=N;i++) if($1==ch[i] && $2>lo[i] && $2<=hi[i]) {
             g=gn[i]; cnt[g]++; sum[g]+=$3;
             if($3>=10) a10[g]++; if($3>=20) a20[g]++; if($3>=30) a30[g]++; if($3==0) z[g]++;
             break } }
       END { print "gene,exonic_bp,mean_depth,pct_ge10x,pct_ge20x,pct_ge30x,pct_zero,callable";
             for (g in cnt) {
               p20 = 100*a20[g]/cnt[g];
               printf "%s,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%s\n", g, cnt[g], sum[g]/cnt[g],
                 100*a10[g]/cnt[g], p20, 100*a30[g]/cnt[g], 100*z[g]/cnt[g],
                 (p20>=90 ? "yes" : (p20>=70 ? "partial" : "NO")) } }' \
      "$OUT/panel.sorted.bed" "$OUT/panel_depth.tsv" \
    | (read -r h; echo "$h"; sort -t, -k5,5gr) > "$DATA/panel_callability.csv"
  rm -f "$OUT/panel_depth.tsv"
else
  echo "  SKIP: $BED not found"
fi

echo
echo "wrote:"
echo "  $DATA/qc_summary.csv"
echo "  $DATA/uba1_m41t_allele_counts.csv"
echo "  $DATA/panel_callability.csv"
echo
echo "The BAM is not copied anywhere and is not committed."
