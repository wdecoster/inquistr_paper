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
- **samtools** on `PATH` — used to subsample the deep ONT alignment (`--subsample` requires
  samtools 1.14 or newer).
- **medaka 2.2.1 or newer** — used only by the ONT medaka tandem comparison, and run from the
  conda environment named in the `run_medaka_tandem` rule. 2.2.1 is the release that added both
  CRAM input and gzipped BED regions to the tandem subcommand; older versions cannot pass a
  reference to htslib and fail every region with `[E::cram_next_slice] Failure to decode slice`.
  Note that medaka 2.2.0 is also the first release supporting Python 3.13, so an older Python can
  pin the installed medaka below the usable version. The environment must also provide `pyabpoa`,
  which the medaka wheel does not pull in and which the default `hybrid` phasing needs for its
  abPOA fallback.
- **Deep HG002 ONT alignment** — `ONT_FULL_CRAM`, roughly 50x, deeper than the `ont.cram` used
  by the rest of the workflow. `ONT_FULL_COVERAGE` records that estimate and sets the
  subsampling fractions.
- **conda environments** — `envs/straglr.yml`, `envs/minimap2.yml`, `envs/medaka.yml` (not
  committed; required for the STRaglr, HG01514 realignment and medaka tandem rules).

The `.vscode/sftp.json` config mirrors the repository to the authors' compute host on save; it is
git-ignored and not needed to run the workflow.

## Target rules

`rule all` aggregates everything below. To build only a portion, ask Snakemake for one of these
convenience targets, e.g. `snakemake -s analysis.smk polymorphic --cores 8`:

| Target rule | Builds |
| --- | --- |
| `all` | Everything (default). |
| `straglr_chr21` | The chr21-only STRaglr vs inquiSTR runtime/accuracy comparison. |
| `stratified_accuracy` | Genotype accuracy versus truth allele length, for every tool and technology. |
| `coverage_titration` | inquiSTR concordance versus ONT sequencing depth, 5x to 50x. |
| `medaka` | medaka tandem on ONT: runtime, memory, accuracy vs truth, concordance with inquiSTR. |
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

# accuracy against the adotto/GIAB truth, stratified by truth allele length
stratify_accuracy_by_length (per technology × tool)
  → combine_stratified_accuracy → tool_comparison/stratified_accuracy.tsv
  → combine_stratified_errors   → tool_comparison/stratified_errors.tsv.gz
      → plot_stratified_accuracy → tool_comparison/stratified_accuracy.html
```

Every tool is scored against the adotto/GIAB HG002 Tier1 truth through the same code path
(TRGT and LongTR via their VCFs converted to inquiSTR format), reproducing the filtering that
`inquiSTR benchmark --tier1 --max-locus` applies. Loci are binned by the **length of the truth
allele** — the reference span plus the length difference of the longer haplotype — using the
edges in `LENGTH_BIN_EDGES` (0, 100, 200, 500, 1000, 2000, 5000 bp, last bin open-ended).
`inquiSTR call --require-spanning` is included as its own call set, which isolates how much of
the behaviour at long alleles comes from soft-clip-derived genotypes; it tracks plain inquiSTR
closely enough that it is written to the table but left off the plot. Which call set represents
each tool is defined in `STRATIFIED_CALLSETS`.

Two details matter for reading the table:

- **The call-rate denominator is the calling catalog, not the truth set.** The truth holds
  ~1.64 M Tier1 loci while the adotto calling catalog holds ~0.94 M, and only ~11% share a start
  coordinate, which is what `inquiSTR benchmark` matches on. Truth loci with no catalog entry
  were never presented to any tool, so counting them as "not genotyped" puts every tool at ~11%.
  `n_targeted` is the number of truth loci in the catalog, `call_rate` is `n_called / n_targeted`,
  and `catalog_coverage` (`n_targeted / n_truth`) records how much of the truth set that leaves.
- **Two tolerances are scored.** `within_tolerance` uses the fixed 3 bp applied elsewhere in the
  paper; `within_scaled_tolerance` uses 3 bp or 1% of the truth allele length, whichever is
  larger. A fixed 3 bp window demands 0.06% accuracy on a 5 kb allele but 6% on a 50 bp one, so
  the two are identical below 300 bp and diverge above it only insofar as the fixed tolerance is
  what drives the decline rather than the genotypes themselves.

- **Per-locus errors are written out but not plotted.** `stratified_errors.tsv.gz` holds the
  signed error (call minus truth) for every genotyped locus, ~1.3M rows across the seven call
  sets, and the summary table carries `p90_abs_error` / `p99_abs_error` alongside the median and
  mean. Error magnitude grows sharply with allele length (mean |error| rises 0.23 → 2204 bp
  across the bins), but that is a different question from the one the figure answers, so the
  distributions are deliberately left out of it. The data is retained so the panel can be added
  back without redoing the analysis. The export button in the HTML is configured for a 3× PNG,
  since the browser default matches the on-screen size and is too coarse for print.

Per bin the summary table also reports the exact / within-1 bp fractions, the median and mean
absolute error, and Pearson *r*.

#### medaka tandem (ONT only, held separate)

Added at a reviewer's request as a stronger ONT alternative to LongTR. It is kept out of the
headline runtime, memory and accuracy figures for now and given its own outputs, the way STRaglr
is, so the metrics can be judged before deciding whether it joins the main comparison.

```
run_medaka_tandem → convert_medaka_to_inquistr_format → ont-medaka-adotto_rep{n}.inq.gz
  → medaka_accuracy_ont              (vs the adotto/GIAB truth)
  → benchmark_inquistr_vs_medaka_ont (inquiSTR vs medaka concordance)
aggregate_medaka_results → tool_comparison/medaka_results.tsv
  → plot_medaka_benchmark → tool_comparison/medaka_benchmark_plot.html
```

`aggregate_medaka_results` reuses the existing ONT inquiSTR and LongTR timings rather than
re-running them, so only medaka is newly timed.

The CLI is `medaka tandem <bam> <ref> <regions.bed> <sex> <outdir>`, where the sample sex is a
required positional (HG002 is male) and `--workers` is the parallelism setting. medaka reads the
bgzipped adotto catalog directly, using only its first three columns, and writes its genotypes to
`<outdir>/medaka_to_ref.TR.vcf`.

`convert_medaka_to_inquistr_format` turns that VCF into the inquiSTR individual-call layout via
`medaka_vcf_to_inquistr`. `inquiSTR convert` knows TRGT and LongTR VCFs but not medaka, and
converting in the workflow avoids touching the Rust binary, which would refresh
`inquiSTR_version.txt` and retrigger every rule depending on it. medaka reports the whole
haplotype-specific repeat sequence as the alternate allele, so the inquiSTR value is the length
of the sequence the genotype points at minus the reference span. Note `<DEL>`, medaka's symbolic
allele for a haplotype where the repeat is entirely absent: its length is 0, giving `-len(REF)`.
`POS - 1` and `len(REF)` reproduce the catalog interval exactly, which is what the benchmark
matches on.

Two settings affect how favourably medaka is judged and are worth revisiting:

- `MEDAKA_MODEL` is empty, so medaka uses its default consensus model. If the ONT reads were
  basecalled with a different model, medaka's consensus and therefore its accuracy are
  understated.
- medaka skips repeats whose estimated allele length exceeds 10 kb unless
  `--process_large_regions` is given, listing them in `<outdir>/skipped_large.bed`. inquiSTR's
  `--max-locus` filters on *reference span* instead, so medaka drops some loci inquiSTR calls.
  The default is kept because the alternative costs 14–23 GB of RAM and would distort the memory
  benchmark; the skipped count is recoverable from that file.

### 3. ONT coverage titration (`coverage/`)

Answers how far concordance holds up as sequencing depth falls. The deep HG002 ONT alignment is
subsampled to 5x-50x in 5x steps, genotyped, and scored against the adotto/GIAB truth.

```
downsample_ont_coverage (per level, temp() CRAM) → call_ont_coverage → benchmark_ont_coverage
  → aggregate_coverage_accuracy → coverage/accuracy_by_coverage.tsv
      → plot_coverage_accuracy  → coverage/accuracy_by_coverage.html
```

Levels are nominal depths derived from `ONT_FULL_COVERAGE` rather than measured ones, and the
subsampling seed is fixed so reruns select the same reads. The downsampled CRAMs are `temp()`:
the ten levels together hold ~5.5x the reads of the source file, and only the genotypes are
needed afterwards. The figure plots the loci assessed alongside the concordance metrics, since
at low depth a locus can fail the read-support threshold outright and concordance is only
computed over loci that produced a genotype.

### 4. chr21 STRaglr comparison (`tool_comparison/…chr21…`)

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

### 5. Polymorphic repeats: relatedness & PCA (`polymorphic/`)

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

### 6. Pathogenic expansions / PureTarget (`puretarget-calls/`)

```
genotype_puretarget (per BAM, inquiSTR --preset pathogenic)
  → combine_puretarget → combined.tsv → plot_puretarget_heatmap → heatmap.html
```

### 7. Tool versions and summary

```
straglr_version + inquiSTR_version → capture_tool_versions → tool_versions.txt
metrics_summary → tool_comparison/metrics_summary.txt
  (speed fold-changes from results.tsv + accuracy %s parsed from the benchmark TSVs)
metrics_accuracy → tool_comparison/metrics_accuracy.txt
  (accuracy %s parsed from the benchmark TSVs only)
```

## Outputs

Key paper artefacts are the HTML plots (`*_plot.html`, `*.html`), the aggregated `results.tsv` /
`accuracy_*.tsv` tables, `tool_comparison/stratified_accuracy.tsv` (accuracy per truth allele
length bin), the relatedness/PCA tables, `tool_comparison/metrics_summary.txt`
(plain-text speed and accuracy summary), and `tool_comparison/metrics_accuracy.txt`
(plain-text accuracy-only summary). All timing files (`*.time`) are produced by
`/usr/bin/time -v` and parsed for wall-clock seconds and peak RSS.

## Histogram GOLGA8A lengths

```bash
cd ~/inquistr_paper/groupplot
./inquiSTR plot --region chr15:34419394-34419476 combined.tsv.gz --output groupplot.svg --sampleinfo cohort.tsv --condition aFTLD_U:PAT,CON --min 30
```
