# inquistr_paper

Workflow and code for the inquiSTR publication.

This repository contains the [Snakemake](https://snakemake.readthedocs.io/) workflow used to
produce the benchmarks, tool comparisons, and analyses reported in the inquiSTR paper. It is
published for transparency and reproducibility — it is **not** intended as a turn-key pipeline:
many paths to tools, references, and input data are hard-coded for the authors' compute
environment (see [Environment assumptions](#environment-assumptions)).

## Workflow files

| File | Contents |
| --- | --- |
| `analysis.smk` | Main workflow: genotyping, benchmarking, tool comparison, and aggregation rules. `include`s `plotting.smk`. |
| `plotting.smk` | Plotting rules that turn the aggregated `*.tsv` results into interactive HTML plots. |

Run the workflow with `analysis.smk` as the Snakefile:

```bash
snakemake -s analysis.smk --use-conda --cores <N>
```

The default target (`rule all`) builds the complete set of paper outputs. Several smaller
**convenience target rules** are provided to build only part of the workflow (see below).

## Environment assumptions

The workflow expects the following to exist at hard-coded locations (defined at the top of
`analysis.smk`). Adapt these before running elsewhere:

- **inquiSTR** binary — `inquiSTR` variable.
- **TRGT**, **LongTR**, **STRaglr** — `TRGT`, `LongTR`, `STRAGLR` variables (STRaglr runs inside
  the `envs/straglr.yml` conda environment).
- **Reference genomes** — `GRCh38.fa` (`reference`) and the 1KG ONT Vienna reference used for the
  polymorphic analysis.
- **Cohort table** — `1000G_cohort.tsv`, a TSV with `sample`, `hg38_path`, `source`, and
  `Superpopulation code` columns. Sample cohorts (PCA, relatedness, genotyping) are derived from
  it deterministically (`.head(...)` rather than random sampling, so reruns don't change the
  selection and re-trigger downstream work).
- **Input alignments** — `pacbio.cram`, `ont.cram` (the downsampled ONT alignment), the adotto STR
  catalog, the HG002 adotto truth BED, the polymorphic-repeats BED, and the PureTarget BAMs.
- **conda environments** — `envs/straglr.yml`, `envs/minimap2.yml` (not committed; required for the
  STRaglr and HG01514 realignment rules).

The `.vscode/sftp.json` config mirrors the repository to the authors' compute host on save; it is
git-ignored and not needed to run the workflow.

## Target rules

`rule all` aggregates everything below. To build only a portion, ask Snakemake for one of these
convenience targets, e.g. `snakemake -s analysis.smk polymorphic --cores 8`:

| Target rule | Builds |
| --- | --- |
| `all` | Everything (default). |
| `straglr_chr21` | The chr21-only STRaglr vs inquiSTR runtime/accuracy comparison. |
| `polymorphic` | The polymorphic-repeat relatedness table and PCA plots. |
| `pathogenic` | inquiSTR genotyping of the PureTarget pathogenic-expansion BAMs + combined table and heatmap. |
| `versions` | Captured tool-version files. |

## Analysis components and their dependencies

The workflow is organised into independent analyses, each with its own genotype → aggregate →
plot chain. Below is the dependency flow for each.

### 1. inquiSTR scaling benchmark (`benchmarking/`)

Measures inquiSTR runtime/memory across 1–12 threads × 5 replicates for ONT and PacBio. Benchmark
rules are serialised with `resources: benchmark_slot=1` so only one runs at a time.

```
benchmark_call (per technology × threads × replicate)
  → aggregate_benchmark_results → benchmarking/results.tsv
      → plot_benchmark_results → benchmarking/runtime_plot.html
      → plot_memory_usage      → benchmarking/memory_plot.html
inquiSTR_accuracy (vs adotto HG002 truth) → benchmarking/accuracy_{technology}.tsv + .html
```

### 2. Tool comparison: inquiSTR vs TRGT vs LongTR (`tool_comparison/`)

Runtime, memory, and genotype-concordance comparison on the full adotto catalog (5 replicates each).

```
download_adotto → adotto_TRGT.bed.gz → decompress_catalog → adotto_LongTR.bed

# per-tool genotyping (timed)
TRGT_adotto / inquiSTR_adotto / LongTR_pacbio          (PacBio)
inquiSTR_ont / LongTR_ont                              (ONT)

# "filtered" pipelines: inquiSTR filters the catalog to variable loci,
# then TRGT/LongTR genotype only those loci
inquiSTR_adotto(_longtr) → filter_inquiSTR_* → {TRGT,LongTR}_*_filtered

# convert TRGT/LongTR VCFs to inquiSTR format (rep1) for concordance
convert_trgt_pacbio / convert_longtr_pacbio / convert_longtr_ont

# runtime/memory aggregation + plots
aggregate_tool_comparison → tool_comparison/results.tsv
  → plot_tool_comparison_time   → runtime_plot.html
  → plot_tool_comparison_memory → memory_plot.html

# genotype concordance (rep1) and --require-spanning variants
benchmark_inquistr_vs_{trgt,longtr}_{pacbio,ont} → *.tsv + .html + _discrepancies.tsv
inquiSTR_adotto_requirespanning / inquiSTR_ont_requirespanning
  → benchmark_inquistr_requirespanning_vs_{trgt,longtr}_* → *.tsv + .html + _discrepancies.tsv
```

### 3. chr21 STRaglr comparison (`tool_comparison/…chr21…`)

STRaglr is slow, so it is benchmarked on a chr21-only catalog against inquiSTR.

```
create_chr21_catalog_for_straglr → straglr_chr21_catalog.bed
  → run_straglr_chr21 / run_inquistr_chr21 (per technology × replicate, timed)
      → aggregate_straglr_chr21_results → straglr_chr21_results.tsv
          → plot_straglr_chr21_runtime / plot_straglr_chr21_memory
  → convert_straglr_to_inquistr_format
      → benchmark_straglr_chr21_accuracy / benchmark_inquistr_chr21_accuracy
        (vs adotto truth) → *_accuracy_{technology}.tsv + .html
```

### 4. Polymorphic repeats: relatedness & PCA (`polymorphic/`)

Genotypes ~100 1000G samples at Illumina polymorphic-repeat loci, then runs inquiSTR `relate` and
`pca`. CRAMs are downloaded per sample, genotyped, and removed. HG01514 is realigned from ONT reads
with minimap2 to demonstrate the father (HG01512)/child relationship.

```
genotype_polymorphic_sample (per sample, downloads CRAM)
realign_HG01514 → genotype_HG01514
  → combine_pca         → repeats_pca_combined.tsv → polymorphic_pca  → repeats_pca.html + scores
                                                       → plot_polymorphic_pca → repeats_pca_colored.html
  → combine_relatedness → repeats_relatedness_combined.tsv → polymorphic_relate → repeats_relate.tsv
```

### 5. Pathogenic expansions / PureTarget (`puretarget-calls/`)

```
genotype_puretarget (per BAM, inquiSTR --preset pathogenic)
  → combine_puretarget → combined.tsv → plot_puretarget_heatmap → heatmap.html
```

### 6. Tool versions and summary

```
straglr_version + inquiSTR_version → capture_tool_versions → tool_versions.txt
metrics_summary → tool_comparison/metrics_summary.txt
  (speed fold-changes from results.tsv + accuracy %s parsed from the benchmark TSVs)
```

## Outputs

Key paper artefacts are the HTML plots (`*_plot.html`, `*.html`), the aggregated `results.tsv` /
`accuracy_*.tsv` tables, the relatedness/PCA tables, and `tool_comparison/metrics_summary.txt`
(plain-text speed and accuracy summary). All timing files (`*.time`) are produced by
`/usr/bin/time -v` and parsed for wall-clock seconds and peak RSS.

## Histogram GOLGA8A lengths

```bash
cd ~/inquistr_paper/groupplot
./inquiSTR plot --region chr15:34419394-34419476 combined.tsv.gz --output groupplot.svg --sampleinfo cohort.tsv --condition aFTLD_U:PAT,CON --min 30
```
