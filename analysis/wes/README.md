# Whole-exome sequencing

Produces the exome QC and locus summary tables.

## Input

A BAM aligned to hg19 / GRCh37. Not included in this repository — raw sequence
data is deposited under controlled access (see the manuscript's Data
Availability statement). `.gitignore` blocks `*.bam`, `*.bai`, `*.cram`,
`*.crai`, `*.sam`, `*.fastq*`, `*.vcf*`.

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
