# Whole-exome sequencing

Produces the exome QC and locus summary tables.

## Input

A BAM aligned to hg19 / GRCh37. Not included in this repository — raw sequence
data is deposited under controlled access (see the manuscript's Data
Availability statement). `.gitignore` blocks `*.bam`, `*.bai`, `*.cram`,
`*.crai`, `*.sam`, `*.fastq*`, `*.vcf*`.

## Reproducing the outputs without the full exome

The complete alignment carries the subject's germline genotype genome-wide and
is available only under controlled access. Derived data sufficient to recompute
both exome tables are committed here:

| File | Contents |
|---|---|
| `../../data/wes/panel_depth_per_base.tsv.gz` | Depth at each of the 100,686 exonic bases of the panel (`chrom`, `pos`, `depth`; 1-based, hg19) |
| `../../data/wes/uba1_locus_depth_per_base.tsv` | Depth at each base of chrX:47058151–47058751 |
| `../../data/wes/uba1_locus.bam`, `.bai` | The 106 reads overlapping *UBA1* c.122, footprint chrX:47058123–47058981 |

Depth was computed with `samtools depth -a -Q 13 -q 20`, the filters used for
the published table. To recompute both tables and diff them against the
committed values:

```bash
./analysis/wes/verify.sh
```

`qc_summary.csv` reports whole-BAM read counts from `samtools flagstat`, which
require the full alignment; its depth fields are recomputable from the above.

The *UBA1* extract is de-identified: the header declares only the contig
dictionary, reference build, platform and the `UPN1` pseudonym; read names are
sequential; instrument, run, chip, barcode, centre and timestamp fields and
vendor flow-signal tags are removed. `tests/test_wes.py` enforces this, caps the
file size, and keeps every other sequence-data extension blocked.

## Run

```bash
./run_wes_qc.sh /path/to/sample.bam
python -m pytest tests/test_wes.py -q
```

Requires `samtools` and `awk`. The BAM is read in place, never copied.

## Outputs

Written to `../../data/wes/`:

| File | Contents |
|---|---|
| `qc_summary.csv` | Read counts, platform, reference build, library type, depth |
| `uba1_m41t_allele_counts.csv` | Allele counts at *UBA1* c.122 across filter settings |
| `panel_callability.csv` | Per-gene exonic coverage over the interrogated panel |

## Notes

`hlh_gene_panel.txt` lists the 27 interrogated genes; `hlh_panel_exons.bed`
holds their canonical-transcript exon coordinates on GRCh37 (344 intervals),
resolved from the Ensembl GRCh37 REST API.

*UBA1* c.122 is chrX:47,058,451 on GRCh37, reference base T (p.Met41Thr = T>C).
Codon 41 is chrX:47,058,450–47,058,452.

In `panel_callability.csv`, `callable` is a threshold on `pct_ge20x`: `yes` at
≥90%, `partial` at 70–90%, `NO` below 70%. Coverage is computed at baseQ ≥13,
mapQ ≥20 against Ensembl GRCh37 canonical transcripts.

In `uba1_m41t_allele_counts.csv`, `vaf_95pct_upper_bound` is the one-sided 95%
Clopper–Pearson limit at the stated depth.
