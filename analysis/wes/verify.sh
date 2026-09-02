#!/usr/bin/env bash
# Recompute the two published whole-exome tables from the per-base depth data
# and the UBA1 locus extract committed in data/wes/, and compare.
#
# Usage:  ./analysis/wes/verify.sh
#
# Requires: samtools >= 1.17, python3
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data/wes"
BED="$ROOT/analysis/wes/hlh_panel_exons.bed"

echo "== 1. panel_callability.csv from panel_depth_per_base.tsv.gz =="
python3 - "$DATA" "$BED" <<'PY'
import bisect
import csv
import gzip
import sys

data, bed_path = sys.argv[1], sys.argv[2]

iv = {}
for row in csv.reader(open(bed_path), delimiter="\t"):
    iv.setdefault(row[0], []).append((int(row[1]), int(row[2]), row[3]))
for c in iv:
    iv[c].sort()
starts = {c: [x[0] for x in v] for c, v in iv.items()}

depths = {}
with gzip.open(f"{data}/panel_depth_per_base.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        ch, pos, d = line.split("\t")
        pos0 = int(pos) - 1
        j = bisect.bisect_right(starts[ch], pos0) - 1
        if j >= 0 and pos0 < iv[ch][j][1]:
            depths.setdefault(iv[ch][j][2], []).append(int(d))

bad = 0
for p in csv.DictReader(open(f"{data}/panel_callability.csv")):
    v = depths[p["gene"]]
    n = len(v)
    calc = {
        "exonic_bp": n,
        "mean_depth": round(sum(v) / n, 1),
        "pct_ge10x": round(sum(x >= 10 for x in v) / n * 100, 1),
        "pct_ge20x": round(sum(x >= 20 for x in v) / n * 100, 1),
        "pct_ge30x": round(sum(x >= 30 for x in v) / n * 100, 1),
        "pct_zero": round(sum(x == 0 for x in v) / n * 100, 1),
    }
    for k, got in calc.items():
        if abs(got - float(p[k])) > 1e-9:
            print(f"   MISMATCH {p['gene']:<10} {k}: recomputed={got} published={p[k]}")
            bad += 1
print(f"   {'FAIL' if bad else 'OK'}: 27 genes x 6 columns, {bad} mismatches")
sys.exit(1 if bad else 0)
PY

echo
echo "== 2. uba1_m41t_allele_counts.csv from uba1_locus.bam =="
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
for Q in 0 13 20; do for q in 0 20; do
  printf '%s,%s,' "$Q" "$q"
  samtools mpileup -r chrX:47058451-47058451 -d 1000000 -Q "$Q" -q "$q" \
    "$DATA/uba1_locus.bam" 2>/dev/null | awk -F'\t' '{print $4"\t"toupper($5)}'
done; done > "$TMP"

python3 - "$TMP" "$DATA/uba1_m41t_allele_counts.csv" <<'PY'
import csv
import sys

pub = list(csv.DictReader(open(sys.argv[2])))
bad = 0
for line in open(sys.argv[1]):
    Qs, qs, payload = line.split(",", 2)
    depth_s, _, bases = payload.strip().partition("\t")
    Q, q, bases = int(Qs), int(qs), bases.upper()
    got = (int(depth_s), sum(bases.count(c) for c in "T.,"), bases.count("C"),
           bases.count("A"), bases.count("G"), bases.count("*"))
    p = next(x for x in pub
             if int(x["base_quality_min"]) == Q and int(x["mapping_quality_min"]) == q)
    exp = (int(p["depth"]), int(p["T_ref"]), int(p["C_alt"]),
           int(p["A"]), int(p["G"]), int(p["deletion"]))
    if got != exp:
        print(f"   MISMATCH Q{Q} q{q}: recomputed={got} published={exp}")
        bad += 1
print(f"   {'FAIL' if bad else 'OK'}: 6 filter settings, {bad} mismatches")
sys.exit(1 if bad else 0)
PY

echo
echo "Both tables reproduce from the committed data."
