# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Snakemake workflow producing the benchmarks, tool comparisons and figures for the inquiSTR
publication. Published for transparency, **not** as a turn-key pipeline: tool paths, references and
input data are hard-coded for the authors' compute host.

## Commands

```bash
snakemake -s analysis.smk --use-conda --cores <N>          # default target: rule all
snakemake -s analysis.smk <target> --cores <N>             # build one analysis
snakemake -s analysis.smk <target> -n                      # dry run; the only cheap validation
snakemake -s analysis.smk <target> --rerun-triggers mtime  # see "Rerun triggers" below
```

Convenience target rules (input-only rules with no output): `all`, `stratified_accuracy`,
`coverage_titration`, `medaka`, `straglr_chr21`, `polymorphic`, `pathogenic`, `versions`.

There is no test suite. A dry run on the compute host is the validation step.

## The data is not on this machine

Inputs (CRAMs, catalogs, truth sets) live on the authors' compute host; `.vscode/sftp.json` mirrors
the repo there on save. **Do not try to run the workflow locally to check a change** — it will fail
on missing inputs. Verify work instead by:

- exercising the module-level helper functions standalone against small synthetic fixtures (they are
  plain Python and importable by slicing them out of `analysis.smk`);
- extracting a `run:` block body and executing it against a synthetic table;
- asking the user to dry-run or run on the host.

## Architecture

### Two files, one namespace

`analysis.smk` is the Snakefile: constants, ~9 module-level helper functions, and ~77 rules.
`plotting.smk` holds ~11 plotly rules and is pulled in by `include:` **at line 8, before the
constants are defined**. Consequently `plotting.smk` must not reference `analysis.smk` constants
anywhere evaluated at parse time (`input`, `output`, `expand`); it currently references none, and
uses literal paths instead. Bodies of `run:` blocks resolve globals at execution time, so helpers
and constants defined later in `analysis.smk` are available there — that is how the plotting and
aggregation rules reach `parse_accuracy_percentages`, `load_truth_alleles`, `medaka_vcf_to_inquistr`
and friends.

### Every caller is scored through `inquiSTR benchmark`

TRGT, LongTR, STRaglr and medaka tandem outputs are converted to the inquiSTR individual-call
format first, so all tools go through one comparison code path rather than per-tool scoring logic.
Conversions live in the workflow (an `awk` one-liner for STRaglr, `medaka_vcf_to_inquistr` for
medaka) rather than in `inquiSTR convert` — see "Rerun triggers".

The `.inq` format is `chromosome begin end info <sample>_H1 <sample>_H2`, preceded by `#` metadata
lines. Two conventions bite if missed:

- Allele values are the **length difference from the reference in bp** (positive = expansion), not
  absolute allele lengths and not repeat-unit counts. `NaN` means uncalled.
- `# file_type=individual_call` is how `benchmark` distinguishes a call file from an adotto-style
  truth BED. It is optional when the file is passed as `--test`, but **required when passed as
  `--truth`** — without it the file falls through to the BED parser and fails with
  "Malformed BED line 1 (expected 9 fields, got 6)".

### Rerun triggers

Twenty rules take `inquiSTR_version.txt` as an input, so rebuilding the inquiSTR binary retriggers
essentially the whole workflow. Prefer solving format and conversion problems in the workflow over
changing the Rust tool.

Snakemake also reruns on rule *code* changes by default, so editing a helper used by a `run:` block
can retrigger expensive genotyping that did not actually change. Pass `--rerun-triggers mtime` when
that is not wanted.

### Timed benchmark rules

Fourteen rules carry `resources: benchmark_slot=1` so only one timed run executes at a time.
`/usr/bin/time -v` writes `.time` files that the `aggregate_*` rules parse for wall-clock seconds
and peak RSS. Runtime and memory comparisons use 5 replicates (`REPLICATES`); accuracy comparisons
use rep1 only, since accuracy does not depend on the replicate or the thread count.

### New tools start as special cases

STRaglr and medaka tandem are deliberately kept out of the headline runtime, memory and accuracy
figures, each with its own aggregation table, plots and target rule. This lets their metrics be
judged before deciding whether they join the main comparison. Follow that pattern when adding a
caller.

### Denominators in the accuracy analyses

The adotto/GIAB truth holds ~1.64 M Tier1 loci while the calling catalog holds ~0.94 M, and
`inquiSTR benchmark` matches on `chromosome:begin` — only ~11% of truth loci share a start
coordinate with the catalog. Loci absent from the catalog were never presented to any caller, so
call rates must use truth ∩ catalog (`n_targeted`) as the denominator, never the full truth set.
Accuracy percentages are computed only over loci that produced a genotype, so they are conditional
on the call rate and the two must be read together; a caller with heavy dropout is otherwise
flattered by being scored on the subset it managed to call.

### Plotting conventions

Plotly figures built inside `run:` blocks and written with `write_html`, styled with
`plot_bgcolor='white'`, explicit font sizes and black axis lines. Newer figures pass
`config={'toImageButtonOptions': {..., 'scale': 3}}` so the browser's PNG export is print-resolution
rather than screen-resolution. Panel letters are added *after* the loop that restyles
`fig['layout']['annotations']`, otherwise that loop overwrites their font.

## Repository conventions

`revisions/` holds manuscript and reviewer-response drafts. It is intentionally **not committed**,
and is untracked rather than gitignored — avoid `git add -A` / `git commit -a`.
