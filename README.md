# VEXAS/MAS — analysis code for the multiomic study

Code and de-identified source data for *"VEXAS and Macrophage Activation
Syndrome: clinical and multiomic approach of a unique hyperinflammatory clinical
dyad."*

| Module | Produces |
|---|---|
| `analysis/cytokines_flow/` | Tables S3–S4, Figure 2C, Figure S4 |
| `analysis/proteomics/` | Figure 2E |
| `analysis/wes/` | Exome QC and locus summary tables |
| `analysis/ddpcr/` | Figure 2A — code and data not yet added |
| `analysis/complement/` | Figure 2D, Figure S5 — code and data not yet added |

## Install

```bash
pip install -r requirements.txt
```

R modules additionally require `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggrepel`.
`analysis/wes/` requires `samtools`.

## Run

```bash
python analysis/cytokines_flow/run_analysis.py    # Tables S3-S4, Figures 2C & S4
Rscript analysis/proteomics/run_pca.R             # Figure 2E
./analysis/wes/run_wes_qc.sh <bam>                # exome summaries

python -m pytest tests/ -q                        # checks published values
```

Results are written to `outputs/`, which is not committed and regenerates on
each run. Each module's README lists its inputs and outputs.

## Layout

```
analysis/<modality>/   code and a README per assay
data/                  de-identified source data
tests/                 test_cytokines_flow.py  - Tables S3-S4, Figures 2C/S4
                       test_proteomics.py      - Figure 2E
                       test_wes.py             - exome summaries
                       test_data_contracts.py  - conventions applying to all data
outputs/               generated results (not committed)
```

## Statistical method

One VEXAS/MAS case with three early time points against four VEXAS controls
without MAS. Each analyte: Crawford–Howell modified t-test, df = 3, two-sided p,
with the standardised difference δ in control-SD units, the Crawford–Garthwaite
Bayesian percentile and 95% credible interval, and a no-overlap ratio.
Benjamini–Hochberg correction across the nine analytes. Cytokines and the two
HLH-like flow metrics on the log10 scale; the CD8+ fraction linear. Time point 1
is the primary comparison; time points 2 and 3 are reported per time point.

Implemented in `analysis/cytokines_flow/singlecase.py`. Equivalent results are
obtainable in R with the `singcar` package (`TD()`, `BTD()`).

## Exclusions

Recorded in `data/cytokines_flow.csv` with machine-readable reasons:

- IFN-γ time points 2 and 3 — anti-IFN-γ antibody interference; only the
  pre-antibody value (491 pg/mL) is analysed.
- One control's second sample, drawn during a disease flare — excluded from the
  control group, used only in the descriptive comparison of Figure S4.

## Data conventions

Subjects appear only as pseudonymous UPN labels; time points are relative labels
(T1, T2, T3), never calendar dates. Enforced on every data file by
`tests/test_data_contracts.py`.

Raw sequence data, the full NPX plate export and the flow-cytometry FCS files
are not in this repository; they are deposited under the accessions given in the
manuscript's Data Availability statement.

## Citation

See `CITATION.cff`. Archived snapshots of this repository are deposited on
Zenodo. Cite the concept DOI, which resolves to the latest version:

    10.5281/zenodo.22256790

The snapshot released alongside the manuscript is v1.0.0, DOI 10.5281/zenodo.22256791.

## License

MIT — see `LICENSE`.
