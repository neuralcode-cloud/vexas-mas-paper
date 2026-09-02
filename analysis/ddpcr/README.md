# Droplet digital PCR — *UBA1* clonal burden

Produces Figure 2A.

**Code and data not yet added to this repository.**

Assay per the manuscript Methods: QX200 droplet system (Bio-Rad), QuantaSoft
v1.7; Sanger sequencing of *UBA1* exon 3.

## Expected layout

```
analysis/ddpcr/compute_vaf.py      droplet counts -> VAF with Poisson intervals
data/ddpcr/droplet_counts.csv      one row per well
```

`droplet_counts.csv` columns, enforced by `tests/test_data_contracts.py`:
`subject` (UPN*), `sample_type`, `timepoint` (relative label, not a date),
`target`, `droplets_positive`, `droplets_negative`, `included`,
`exclusion_reason`.
