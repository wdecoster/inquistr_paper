import pandas as pd
import glob
import os
import random 
import re

# Include plotting rules from separate file
include: "plotting.smk"

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"
# Load a TSV file into a DataFrame, containing URLs of 1000 Genomes cram files
df = pd.read_csv(file_path, sep='\t', usecols=['sample', 'hg38_path', 'source', 'Superpopulation code'])
df = df[df["source"] == "Noyvert/Schloissnig"] # these are cram files

# PCA cohort: 100 mixed samples, 20 per superpopulation
# we can't take a random sample, because that would make the result of the workflow change every iteration, triggering a rerun
samples_per_superpopulation = 20
selected_samples = pd.concat([
    df[df["Superpopulation code"] == "AFR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "AMR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "EAS"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "EUR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "SAS"].head(samples_per_superpopulation),
]).drop_duplicates(subset="sample")

# Relatedness cohort: 100 EUR samples, must include HG01512 (father of HG01514) and realigned HG01514
# The first 20 EUR overlap with the PCA cohort and don't need to be genotyped twice
relatedness_samples = pd.concat([
    df[df["sample"] == "HG01512"],
    df[df["Superpopulation code"] == "EUR"],
]).drop_duplicates(subset="sample").head(100)
# HG01514 is added via the realign rule and included in the relatedness combine

# All unique samples needing genotyping (union of PCA and relatedness cohorts)
# HG01514 has special handling via the realign rule
all_genotype_samples = pd.concat([selected_samples, relatedness_samples]).drop_duplicates(subset="sample")

# Build lookup dicts for per-sample polymorphic genotyping
POLYMORPHIC_SAMPLE_URLS = {row['sample']: row['hg38_path'] for _, row in all_genotype_samples.iterrows()}
POLYMORPHIC_SAMPLES = list(POLYMORPHIC_SAMPLE_URLS.keys())

# Sample lists for separate combine steps
PCA_SAMPLES = list(selected_samples['sample'])
RELATEDNESS_SAMPLES = list(relatedness_samples['sample'])


reference = "/home/AD/wdecoster/database/GRCh38.fa"
inquiSTR = "/home/AD/wdecoster/repositories/inquiSTR/target/x86_64-unknown-linux-musl/release/inquiSTR"
TRGT = "/home/AD/wdecoster/bin/trgt",
LongTR = "/home/AD/wdecoster/anaconda3/envs/longtr/bin/LongTR"
STRAGLR = "/home/AD/wdecoster/repositories/straglr/straglr.py" # requires environment with straglr dependencies
# medaka tandem, ONT only. Assumed to be on PATH; point this at an absolute path, or add a
# conda directive to the medaka rules the way the STRaglr rules use envs/straglr.yml, if it
# lives in its own environment. medaka tandem additionally needs pyabpoa installed, which the
# medaka wheel does not pull in: without it the default `hybrid` phasing cannot fall back to
# abPOA clustering for loci where haplotags are insufficient.
MEDAKA = "medaka"
# Leave empty to use medaka's default consensus model, or set the model matching the basecaller
# used for the ONT reads (e.g. "dna_r10.4.1_e8.2_400bps_hac@v4.1.0:consensus").
MEDAKA_MODEL = ""

# Benchmark parameters
TECHNOLOGIES = ["ont", "pacbio"]
THREAD_COUNTS = list(range(1, 13))  # 1 to 12 threads
REPLICATES = list(range(1,6))  # 5 replicates
MAX_LOCUS = 10000  # limit to loci shorter than 10kb for genotyping

# Truth set used for all accuracy benchmarking: the adotto/GIAB HG002 tandem repeat benchmark
ADOTTO_TRUTH = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz"

# Length-stratified accuracy analysis: the call set representing each tool per technology.
# rep1 is representative, since accuracy does not depend on the replicate. The TRGT/LongTR
# call sets are their VCFs converted to inquiSTR format, so every tool is evaluated through
# exactly the same code. All of these were genotyped from the same adotto catalog
# coordinates, which keeps the per-length-bin locus sets identical across tools.
STRATIFIED_CALLSETS = {
    ("pacbio", "inquiSTR"): "tool_comparison/pacbio-inquistr-adotto_rep1.inq.gz",
    ("ont", "inquiSTR"): "tool_comparison/ont-inquistr-adotto_rep1.inq.gz",
    ("pacbio", "inquiSTR-spanning"): "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz",
    ("ont", "inquiSTR-spanning"): "tool_comparison/ont-inquistr-adotto-requirespanning.inq.gz",
    ("pacbio", "TRGT"): "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
    ("pacbio", "LongTR"): "tool_comparison/pacbio-longtr-adotto_rep1.inq.gz",
    ("ont", "LongTR"): "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
}

# Lower edges (bp) of the truth allele length bins; the last bin is open-ended
LENGTH_BIN_EDGES = [0, 100, 200, 500, 1000, 2000, 5000]

# Coverage titration. This HG002 ONT alignment is deeper than the downsampled `ont.cram` used
# elsewhere in the paper; ONT_FULL_COVERAGE is its approximate depth and sets the subsampling
# fractions, so correct it here if the estimate changes.
ONT_FULL_CRAM = "/home/AD/wdecoster/optimize_inquiSTR/benchmark/hg002_ont/ont.cram"
ONT_FULL_COVERAGE = 50
# 7.5x is included because the loci-genotyped curve does nearly all of its rising between
# 5x and 10x, so that interval is the one worth resolving.
COVERAGE_LEVELS = [5, 7.5] + list(range(10, ONT_FULL_COVERAGE + 1, 5))


def coverage_label(value):
    """Filename-safe label for a coverage level: 5 -> "5", 7.5 -> "7.5"."""
    return f"{value:g}"

# Randomize order of thread counts for benchmarking
THREAD_ORDER = THREAD_COUNTS * len(REPLICATES)
random.shuffle(THREAD_ORDER)

# if the "tool_versions.txt" file exists, remove it to ensure fresh capture of tool versions
if os.path.exists("tool_versions.txt"):
    os.remove("tool_versions.txt")

# pathogenic expansions to genotype are in /home/AD/wdecoster/inquiSTR_paper/puretarget-data and have to be globbed to a list and genotyped with inquiSTR --pathogenic
puretarget_files = [os.path.basename(f).replace(".bam", "") for f in glob.glob("/home/AD/wdecoster/inquiSTR_paper/puretarget-data/*.bam")]


def parse_accuracy_percentages(path):
    with open(path, "r") as handle:
        text = handle.read()

    patterns = {
        "within_3bp": r"Within\s+3\s+bp\s+tolerance:\s*\d+\s*\(([\d.]+)%\)",
        "within_1bp": r"Maximally\s+off\s+by\s+one:\s*\d+\s*\(([\d.]+)%\)",
        "exact": r"Exact\s+matches:\s*\d+\s*\(([\d.]+)%\)",
    }

    parsed = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if not match:
            raise ValueError(f"Could not parse '{key}' percentage from {path}")
        parsed[key] = float(match.group(1))
    return parsed


def length_bin_labels(edges=None):
    """Labels for the truth allele length bins defined by `edges` (lower edges, bp).

    The final bin is open-ended, e.g. [0, 100, 200] -> ["0-100", "100-200", ">200"].
    """
    edges = LENGTH_BIN_EDGES if edges is None else edges
    labels = [f"{low}-{high}" for low, high in zip(edges, edges[1:])]
    labels.append(f">{edges[-1]}")
    return labels


def load_truth_alleles(bed_path, max_locus, mode="MAX"):
    """Load the adotto/GIAB HG002 truth BED the way `inquiSTR benchmark` does.

    The BED has 9 columns: column 4 holds the tier and the last two columns hold the
    per-haplotype length difference relative to the reference, with the opposite sign
    convention from inquiSTR (hence the negation). Only Tier1 loci are kept and loci
    with a reference span above `max_locus` are dropped, mirroring `--tier1 --max-locus`.

    Returns a DataFrame with the selected truth allele (MAX or MIN of both haplotypes,
    in bp relative to the reference) and its absolute length, i.e. the reference span
    plus that difference. The absolute length is what the stratification bins on.
    """
    import pandas as pd

    truth = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3, 7, 8])
    truth.columns = ["chromosome", "begin", "end", "tier", "h1", "h2"]
    truth[["h1", "h2"]] = -truth[["h1", "h2"]]
    truth = truth[(truth["tier"] == "Tier1") & (truth["end"] - truth["begin"] <= max_locus)]

    alleles = truth[["h1", "h2"]]
    truth["truth_allele"] = alleles.max(axis=1) if mode == "MAX" else alleles.min(axis=1)
    truth = truth.dropna(subset=["truth_allele"])
    truth["truth_length"] = ((truth["end"] - truth["begin"]) + truth["truth_allele"]).clip(lower=0)
    return truth[["chromosome", "begin", "end", "truth_allele", "truth_length"]]


def load_catalog_keys(catalog_path):
    """Load the `chromosome:begin` keys of the catalog the tools were asked to genotype.

    The adotto truth set covers roughly 1.6M Tier1 loci while the calling catalog holds
    ~0.9M, and only ~11% share a start coordinate. Truth loci with no catalog entry were
    never presented to any tool, so they have to be excluded from the call-rate denominator
    or every tool looks like it drops ~90% of loci. The two ONT catalogs
    (`adotto_TRGT.bed.gz`, `adotto_LongTR.bed`) carry identical coordinates and differ only
    in the fourth column, so one catalog serves every call set.
    """
    import pandas as pd

    catalog = pd.read_csv(catalog_path, sep="\t", header=None, usecols=[0, 1], comment="#")
    catalog.columns = ["chromosome", "begin"]
    return pd.MultiIndex.from_arrays([catalog["chromosome"], catalog["begin"]])


def load_call_alleles(call_path, mode="MAX"):
    """Load an inquiSTR individual call file (`.inq`, optionally gzipped) and reduce each
    locus to a single allele (MAX or MIN of both haplotypes), as `inquiSTR benchmark` does.

    TRGT and LongTR calls are converted to this format first, so the same loader serves
    every tool. Loci without a genotype keep their NaN and are counted as not called.
    """
    import gzip
    import pandas as pd

    opener = gzip.open if str(call_path).endswith(".gz") else open
    metadata_lines = 0
    with opener(call_path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            metadata_lines += 1

    calls = pd.read_csv(call_path, sep="\t", skiprows=metadata_lines, usecols=[0, 1, 2, 4, 5])
    calls.columns = ["chromosome", "begin", "end", "h1", "h2"]
    alleles = calls[["h1", "h2"]]
    calls["call"] = alleles.max(axis=1) if mode == "MAX" else alleles.min(axis=1)
    # keep="last" mirrors the HashMap insertion order of `inquiSTR benchmark`
    calls = calls.drop_duplicates(subset=["chromosome", "begin"], keep="last")
    return calls[["chromosome", "begin", "call"]]


def medaka_vcf_to_inquistr(vcf_path, out_path, sample_name="medaka"):
    """Convert a medaka tandem VCF (`medaka_to_ref.TR.vcf`) to inquiSTR individual-call format.

    medaka reports the whole haplotype-specific repeat sequence as the alternate allele, so an
    allele length is the length of the sequence its genotype points at and the inquiSTR value is
    that minus the reference span (positive = expansion). Three cases need care:

    * `<DEL>` is medaka's symbolic allele for a haplotype in which the repeat is entirely absent,
      so its length is 0 and the value is -len(REF), not len("<DEL>"). It appears both as the
      only alternate and as one of several.
    * A missing genotype allele (`.`) becomes NaN, which inquiSTR reads as uncalled.
    * A haploid genotype (one allele, e.g. chrX in a male sample) fills H1 and leaves H2 NaN,
      matching what `inquiSTR call` writes for haploid chromosomes.

    Loci medaka omits, including the >10 kb repeats it lists in `skipped_large.bed`, simply do
    not appear and count as uncalled in the benchmark. `POS - 1` and `len(REF)` reproduce the
    catalog interval exactly, which is what `inquiSTR benchmark` matches on.
    """
    import gzip

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    out_opener = gzip.open if str(out_path).endswith(".gz") else open

    with opener(vcf_path, "rt") as handle, out_opener(out_path, "wt") as out:
        # The `# file_type=` line is what `inquiSTR benchmark` uses to tell an individual call
        # file from an adotto-style truth BED. Without it the file is still usable as --test,
        # but passing it as --truth falls through to the BED parser and fails with
        # "Malformed BED line 1 (expected 9 fields, got 6)".
        out.write("# file_type=individual_call\n")
        out.write("# source=medaka tandem, converted by the inquiSTR paper workflow\n")
        out.write(f"# sample={sample_name}\n")
        out.write("chromosome\tbegin\tend\tinfo\t{0}_H1\t{0}_H2\n".format(sample_name))
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chromosome, pos, ref, alt, fmt, sample = (
                fields[0], int(fields[1]), fields[3], fields[4], fields[8], fields[9]
            )

            genotype = dict(zip(fmt.split(":"), sample.split(":"))).get("GT", ".")
            alleles = [ref] + ([] if alt == "." else alt.split(","))

            values = []
            for call in genotype.replace("|", "/").split("/"):
                if not call.isdigit() or int(call) >= len(alleles):
                    values.append("NaN")
                    continue
                sequence = alleles[int(call)]
                length = 0 if sequence == "<DEL>" else len(sequence)
                values.append(f"{length - len(ref):.1f}")

            while len(values) < 2:
                values.append("NaN")

            begin = pos - 1
            out.write(
                f"{chromosome}\t{begin}\t{begin + len(ref)}\t.\t{values[0]}\t{values[1]}\n"
            )


def parse_loci_assessed(path):
    """Number of loci matched between the call set and the truth, from a benchmark output."""
    with open(path, "r") as handle:
        match = re.search(r"^LOCI_ASSESSED:\s*(\d+)", handle.read(), flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Could not parse LOCI_ASSESSED from {path}")
    return int(match.group(1))


def render_accuracy_table(rows):
    """Render the genotype-accuracy table as text lines.

    `rows` is a list of (technology, test_set, truth_set, metrics) tuples, where metrics is
    the dict from parse_accuracy_percentages. The test set is the caller being evaluated; the
    truth set is either the adotto/GIAB HG002 truth or another caller.
    """
    headers = [
        "Technology", "Test set", "Truth set", "Within 3 bp", "Within 1 bp", "Exact matches",
    ]
    table = [
        [
            technology,
            test_set,
            truth_set,
            f"{metrics['within_3bp']:.2f}%",
            f"{metrics['within_1bp']:.2f}%",
            f"{metrics['exact']:.2f}%",
        ]
        for technology, test_set, truth_set, metrics in rows
    ]
    widths = [len(h) for h in headers]
    for row in table:
        for i, value in enumerate(row):
            widths[i] = max(widths[i], len(value))

    def fmt(values):
        return " | ".join(value.ljust(widths[i]) for i, value in enumerate(values))

    out = [
        "Accuracy table",
        "--------------",
        fmt(headers),
        "-+-".join("-" * width for width in widths),
    ]
    out.extend(fmt(row) for row in table)
    return out


rule all:
    input:
        "benchmarking/results.tsv",
        "benchmarking/runtime_plot.html",
        "benchmarking/memory_plot.html",
        #"genotyping/adotto_combined_selected_samples.tsv",
        # PacBio TRGT
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR filter catalog (TRGT format)
        expand("tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-pacbio_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR + filter catalog (LongTR format)
        expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time", replicate=REPLICATES),
        # PacBio TRGT filtered
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # PacBio LongTR
        expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio LongTR filtered
        expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR
        expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR filter catalog (TRGT format, unused by LongTR)
        expand("tool_comparison/adotto-variable-catalog-ont_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-ont_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR + filter catalog (LongTR format)
        expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time", replicate=REPLICATES),
        # ONT LongTR
        expand("tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # ONT LongTR filtered
        expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # inquiSTR convert (rep1 only - representative conversion)
        "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
        "tool_comparison/pacbio-longtr-adotto_rep1.inq.gz",
        "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
        # inquiSTR benchmark: compare inquiSTR genotypes against TRGT/LongTR genotypes
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.html",
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio_discrepancies.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.html",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio_discrepancies.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont.html",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont_discrepancies.tsv",
        # inquiSTR --require-spanning vs TRGT benchmark
        "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.html",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio_discrepancies.tsv",
        # inquiSTR --require-spanning vs LongTR benchmark (ONT)
        "tool_comparison/ont-inquistr-adotto-requirespanning.inq.gz",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont.tsv",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont.html",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont_discrepancies.tsv",
        expand("benchmarking/accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        # TRGT and LongTR genotype accuracy vs the adotto/GIAB HG002 truth
        "tool_comparison/trgt_accuracy_pacbio.tsv",
        "tool_comparison/trgt_accuracy_pacbio.html",
        expand("tool_comparison/longtr_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/longtr_accuracy_{technology}.html", technology=TECHNOLOGIES),
        # genotype accuracy stratified by truth allele length
        "tool_comparison/stratified_accuracy.tsv",
        "tool_comparison/stratified_accuracy.html",
        "tool_comparison/stratified_errors.tsv.gz",
        # genotype accuracy versus ONT sequencing depth
        "coverage/accuracy_by_coverage.tsv",
        "coverage/accuracy_by_coverage.html",
        "tool_comparison/results.tsv",
        "tool_comparison/runtime_plot.html",
        "tool_comparison/memory_plot.html",
        "tool_comparison/straglr_chr21_catalog.bed",
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        "tool_comparison/straglr_chr21_results.tsv",
        "tool_comparison/straglr_chr21_runtime_plot.html",
        "tool_comparison/straglr_chr21_memory_plot.html",
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        "tool_versions.txt",
        "straglr_version.txt",
        "inquiSTR_version.txt",
        "tool_comparison/metrics_summary.txt"

rule straglr_chr21:
    """Target rule for the chr21 STRaglr vs inquiSTR runtime comparison subset."""
    input:
        "tool_comparison/straglr_chr21_catalog.bed",
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.tsv", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.inq.gz", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.inq.gz", technology=TECHNOLOGIES, replicate=["1"]),
        "tool_comparison/straglr_chr21_results.tsv",
        "tool_comparison/straglr_chr21_runtime_plot.html",
        "tool_comparison/straglr_chr21_memory_plot.html",
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES)

rule versions:
    input:
        "tool_versions.txt",
        "straglr_version.txt"

rule benchmark_callers:
    input:
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        "tool_comparison/metrics_accuracy.txt"

rule straglr_version:
    """Capture STRaglr version using the STRaglr conda environment."""
    input:
        STRAGLR
    output:
        "straglr_version.txt"
    params:
        straglr = STRAGLR
    conda:
        "envs/straglr.yml"
    log:
        "logs/straglr_version.log"
    shell:
        """
        {params.straglr} --version > {output} 2> {log} || {params.straglr} -h > {output} 2>> {log}
        """

rule inquiSTR_version:
    """Capture inquiSTR version to trigger reruns when version changes"""
    input:
        inquiSTR
    output:
        "inquiSTR_version.txt"
    params:
        inquiSTR = inquiSTR
    log:
        "logs/inquiSTR_version.log"
    shell:
        """
        {params.inquiSTR} --version > {output} 2>&1
        """

rule capture_tool_versions:
    input:
        straglr_version = "straglr_version.txt"
    output:
        "tool_versions.txt"
    params:
        inquiSTR = inquiSTR,
        TRGT = TRGT,
        LongTR = LongTR
    log:
        "logs/capture_tool_versions.log"
    shell:
        """
        {{
            echo "Tool Versions Summary"
            echo "====================="
            echo ""
            echo "inquiSTR:"
            {params.inquiSTR} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "TRGT:"
            {params.TRGT} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "LongTR:"
            {params.LongTR} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "STRaglr:"
            cat {input.straglr_version}
            echo ""
            echo "Generated on: $(date)"
        }} > {output} 2> {log}
        """

# using polymorphic repeats from illumina https://zenodo.org/records/8329210/files/polymorphic_repeats.hg38.bed?download=1
rule polymorphic:
    input:
        "polymorphic/repeats_relate.tsv",
        "polymorphic/repeats_pca.html",
        "polymorphic/repeats_pca_colored.html",

rule create_manifest_selected_samples:
    output:
        manifest = "selected_samples_manifest.tsv"
    log:
        "logs/create_manifest.log"
    run:
        with open(output.manifest, 'w') as f:
            f.write("bam_path\tsample_name\n")
            for idx, row in selected_samples.iterrows():
                f.write(f"{row['hg38_path']}\t{row['sample']}\n")

rule genotype_polymorphic_sample:
    input:
        bed = "polymorphic_repeats.hg38.bed",
    output:
        "polymorphic/individual/{sample}.tsv"
    params:
        url = lambda wildcards: POLYMORPHIC_SAMPLE_URLS[wildcards.sample],
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        cram = lambda wildcards: f"polymorphic/tmp/{wildcards.sample}.cram"
    threads: 4
    log:
        "logs/genotype_polymorphic_{sample}.log"
    shell:
        """
        mkdir -p polymorphic/tmp 2> {log}
        wget -q -O {params.cram} {params.url} 2>> {log}
        wget -q -O {params.cram}.crai {params.url}.crai 2>> {log}
        {params.inquiSTR} call {params.cram} \
            --region-file {input.bed} \
            --reference {params.reference} \
            --threads {threads} \
            --unphased \
            --max-locus {params.max_locus} > {output} 2>> {log}
        rm -f {params.cram} {params.cram}.crai 2>> {log}
        """

rule realign_HG01514:
    input:
        "/home/AD/wdecoster/study322-ONT_genomes/all_files/rr_HG01514/LCYT/209418/v7.3.11/cram/rr_HG01514§LCYT.cram"
    output:
        cram = "polymorphic/tmp/HG01514_realigned.cram",
        crai = "polymorphic/tmp/HG01514_realigned.cram.crai"
    params:
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
    threads:
        20
    conda:
        "envs/minimap2.yml" # also includes samtools
    log:
        "logs/realign_HG01514.log"
    shell:
        """
        mkdir -p polymorphic/tmp
        samtools fastq -@ {threads} {input} 2>> {log} \
        | minimap2 -ax map-ont -t {threads} {params.reference} - 2>> {log} \
        | samtools sort --write-index -o {output.cram} - 2>> {log}
        """

rule genotype_HG01514:
    # run inquiSTR call on the realigned HG01514
    input:
        bed = "polymorphic_repeats.hg38.bed",
        cram = "polymorphic/tmp/HG01514_realigned.cram",
        crai = "polymorphic/tmp/HG01514_realigned.cram.crai"
    output:
        "polymorphic/individual/HG01514.tsv"
    params:
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    threads: 4
    log:
        "logs/genotype_HG01514.log"
    shell:
        """
        {params.inquiSTR} call {input.cram} \
            --region-file {input.bed} \
            --reference {params.reference} \
            --threads {threads} \
            --unphased \
            --max-locus {params.max_locus} > {output} 2> {log}
        """


rule combine_pca:
    input:
        expand("polymorphic/individual/{sample}.tsv", sample=PCA_SAMPLES),
    output:
        "polymorphic/repeats_pca_combined.tsv"
    log:
        "logs/combine_pca.log"
    threads:
        4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input} > {output} 2> {log}
        """

rule combine_relatedness:
    input:
        samples = expand("polymorphic/individual/{sample}.tsv", sample=RELATEDNESS_SAMPLES),
        hg01514 = "polymorphic/individual/HG01514.tsv"
    output:
        "polymorphic/repeats_relatedness_combined.tsv"
    log:
        "logs/combine_relatedness.log"
    threads:
        4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input.samples} {input.hg01514} > {output} 2> {log}
        """



rule polymorphic_relate:
    input:
        combined = "polymorphic/repeats_relatedness_combined.tsv",
    output:
        "polymorphic/repeats_relate.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR,
        min_spacing = 100000, # also the default, but set explicitly for transparency
        tolerance = 1, # also the default, but set explicitly for transparency
    log:
        "logs/polymorphic_relate.log"
    shell:
        """
        {params.inquiSTR} relate \
            --output {output} \
            --threads {threads} \
            --min-spacing {params.min_spacing} \
            --tolerance {params.tolerance} \
            {input.combined} \
            > {log} 2>&1
        """


rule polymorphic_pca:
    input:
        combined = "polymorphic/repeats_pca_combined.tsv",
    output:
        plot = "polymorphic/repeats_pca.html",
        scores = "polymorphic/repeats_pca_scores.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/polymorphic_pca.log"
    shell:
        """
        {params.inquiSTR} pca \
            --output {output.plot} \
            --threads {threads} \
            --scores {output.scores} \
            {input.combined} \
            > {log} 2>&1
        """


rule genotype_adotto:
    input:
        manifest = "selected_samples_manifest.tsv",
        bed = "adotto.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        "genotyping/adotto_combined_selected_samples.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads: 16
    log:
        "logs/genotype_adotto.log"
    shell:
        """
        {params.inquiSTR} batch {input.manifest} \
            --region-file {input.bed} \
            --output {output} \
            --noextend \
            --threads {threads} \
            --parallel-samples 4 \
            --max-locus {params.max_locus} \
            --reference {params.reference} \
            --resume \
            --save-individual genotyping/adotto_individual \
            > {log} 2>&1
        """


# Benchmark rules
rule benchmark_call:
    input:
        cram = "{technology}.cram", # for ont, this refers to the `ont_downsampled.cram` file
        bed = "adotto.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        result = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.tsv",
        timing = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.time"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads: lambda wildcards: int(wildcards.threads)
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    log:
        "logs/benchmark_{technology}_threads{threads}_rep{replicate}.log"
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.cram} \
            --region-file {input.bed} \
            --threads {threads} \
            --noextend \
            --max-locus {params.max_locus} \
            --reference {params.reference} \
            > {output.result} 2> {log}
        """


rule aggregate_benchmark_results:
    input:
        expand("benchmarking/data/{technology}_threads{threads}_rep{replicate}.time",
               technology=TECHNOLOGIES,
               threads=THREAD_COUNTS,
               replicate=REPLICATES)
    output:
        "benchmarking/results.tsv"
    run:
        import re
        results = []
        
        for timing_file in input:
            # Parse filename to get metadata
            match = re.search(r'benchmarking/data/(\w+)_threads(\d+)_rep(\d+)\.time', timing_file)
            if match:
                technology, threads, replicate = match.groups()
                
                # Parse timing output from /usr/bin/time -v
                elapsed_time = None
                max_memory_kb = None
                with open(timing_file, 'r') as f:
                    for line in f:
                        if 'Elapsed (wall clock) time' in line:
                            # Format is "Elapsed (wall clock) time (h:mm:ss or m:ss): 6:15.40"
                            # Split at "): " to get the time value
                            time_str = line.split('): ', 1)[1].strip()
                            # Convert to seconds
                            parts = time_str.split(':')
                            if len(parts) == 3:  # h:mm:ss
                                h, m, s = parts
                                elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                            elif len(parts) == 2:  # mm:ss
                                m, s = parts
                                elapsed_time = int(m) * 60 + float(s)
                            elif len(parts) == 1:  # just ss
                                elapsed_time = float(parts[0])
                        elif 'Maximum resident set size' in line:
                            # Format is "Maximum resident set size (kbytes): 123456"
                            max_memory_kb = int(line.split(':')[1].strip())
                
                if elapsed_time is not None and max_memory_kb is not None:
                    results.append({
                        'technology': technology,
                        'threads': int(threads),
                        'replicate': int(replicate),
                        'elapsed_seconds': elapsed_time,
                        'max_memory_gb': max_memory_kb / (1024 * 1024)  # Convert KB to GB
                    })
        
        # Write results
        import pandas as pd
        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'threads', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)


rule download_adotto:
    output:
        catalog = "adotto_TRGT.bed.gz",
    log:
        "logs/download_adotto.log"
    params:
        url = "https://zenodo.org/records/8329210/files/adotto_repeats.hg38.bed.gz",
    shell:
        """
        wget -O {output.catalog} {params.url} &> {log}
        """

rule decompress_catalog:
    input:
        catalog = "adotto_TRGT.bed.gz"
    output:
        catalog = "adotto_LongTR.bed"
    shell:
        "zcat {input.catalog} | awk 'BEGIN{{OFS=\"\\t\"}} {{match($4, /MOTIFS=([^;]+)/, m); $4=m[1]; print}}' > {output.catalog}"

rule TRGT_adotto:
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.time"
    log:
        "logs/TRGT_adotto_rep{replicate}.log"
    params:
        TRGT = "/home/AD/wdecoster/bin/trgt",
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.TRGT} genotype --genome {params.reference} \
            --repeats {input.catalog} \
            --reads {input.pacbio} \
            --threads {threads} \
            --output-prefix tool_comparison/pacbio-trgt-adotto_rep{wildcards.replicate} &> {log}
        """

rule inquiSTR_adotto:
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
        timing = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time"
    log:
        "logs/inquiSTR_adotto_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads:
        4
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """

rule filter_inquiSTR_pacbio:
    input:
        pacbio_inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-pacbio_rep{replicate}.time"
    log:
        "logs/filter_inquistr_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.pacbio_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' > {output.catalog} 2>> {log}"""

rule inquiSTR_adotto_longtr:
    """Run inquiSTR call on PacBio using the LongTR-format catalog, so the filtered
    output catalog has LongTR-compatible 4th-field motifs."""
    input:
        catalog = "adotto_LongTR.bed",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz",
        timing = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time"
    log:
        "logs/inquiSTR_adotto_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule filter_inquiSTR_pacbio_longtr:
    """Filter the LongTR-catalog inquiSTR PacBio calls to variable loci,
    producing a LongTR-compatible catalog (motif-only 4th field)."""
    input:
        pacbio_inq = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time"
    log:
        "logs/filter_inquistr_pacbio_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.pacbio_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' | grep -v '^chromosome' > {output.catalog} 2>> {log}"""


rule TRGT_adotto_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time"
    log:
        "logs/TRGT_adotto_filtered_rep{replicate}.log"
    params:
        TRGT = "/home/AD/wdecoster/bin/trgt",
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.TRGT} genotype --genome {params.reference} \
            --repeats {input.catalog} \
            --reads {input.pacbio} \
            --threads {threads} \
            --output-prefix tool_comparison/pacbio-trgt-adotto-filtered_rep{wildcards.replicate} &> {log}
        """

rule LongTR_pacbio:
    input:
        catalog = "adotto_LongTR.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.time"
    log:
        "logs/LongTR_pacbio_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.pacbio} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule LongTR_pacbio_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time"
    log:
        "logs/LongTR_pacbio_filtered_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.pacbio} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule inquiSTR_ont:
    input:
        catalog = "adotto_TRGT.bed.gz",
        ont = "ont.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz",
        timing = "tool_comparison/ont-inquistr-adotto_rep{replicate}.time"
    log:
        "logs/inquiSTR_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.ont} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """

rule filter_inquiSTR_ont:
    input:
        ont_inq = "tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-ont_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-ont_rep{replicate}.time"
    log:
        "logs/filter_inquistr_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.ont_inq} --minchange 20 2> {log} | cut -f1-4 > {output.catalog} 2>> {log}
        """

rule inquiSTR_ont_longtr:
    """Run inquiSTR call on ONT using the LongTR-format catalog, so the filtered
    output catalog has LongTR-compatible 4th-field motifs."""
    input:
        catalog = "adotto_LongTR.bed",
        ont = "ont.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz",
        timing = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time"
    log:
        "logs/inquiSTR_ont_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.ont} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule filter_inquiSTR_ont_longtr:
    """Filter the LongTR-catalog inquiSTR ONT calls to variable loci,
    producing a LongTR-compatible catalog (motif-only 4th field)."""
    input:
        ont_inq = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz"
    output:
        catalog = "tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time"
    log:
        "logs/filter_inquistr_ont_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \\
        {params.inquiSTR} filter {input.ont_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' | grep -v '^chromosome' > {output.catalog} 2>> {log}
        """

rule LongTR_ont:
    input:
        catalog = "adotto_LongTR.bed",
        ont = "ont.cram",
    output:
        vcf = "tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/ont-longtr-adotto_rep{replicate}.time"
    log:
        "logs/LongTR_ont_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.ont} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule LongTR_ont_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed",
        ont = "ont.cram",
    output:
        vcf = "tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time"
    log:
        "logs/LongTR_ont_filtered_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.ont} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule convert_trgt_pacbio:
    """
    Convert TRGT VCF to inquiSTR format. The genotyper is auto-detected from the VCF
    header; TRGT coordinates (POS = catalog start) and allele lengths (AL field) are
    handled accordingly.
    """
    input:
        vcf = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz"
    output:
        inq = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_trgt_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule convert_longtr_pacbio:
    input:
        vcf = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz",
    output:
        inq = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_longtr_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule convert_longtr_ont:
    input:
        vcf = "tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz",
    output:
        inq = "tool_comparison/ont-longtr-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_longtr_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule aggregate_tool_comparison:
    input:
        # PacBio timing files
        pacbio_inquistr_time = expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_trgt_time = expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_time = expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_filter_time = expand("tool_comparison/filter-inquistr-pacbio_rep{replicate}.time", replicate=REPLICATES),
        pacbio_trgt_filtered_time = expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_filtered_time = expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_inquistr_time = expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_filter_time = expand("tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time", replicate=REPLICATES),
        # ONT timing files
        ont_inquistr_time = expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_time = expand("tool_comparison/ont-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        ont_filter_time = expand("tool_comparison/filter-inquistr-ont_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_filtered_time = expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_inquistr_time = expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_filter_time = expand("tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time", replicate=REPLICATES)
    output:
        "tool_comparison/results.tsv"
    run:
        import re
        import pandas as pd

        def parse_timing_file(filepath):
            """Parse /usr/bin/time -v output file."""
            elapsed_time = None
            max_memory_kb = None
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Elapsed (wall clock) time' in line:
                        time_str = line.split('): ', 1)[1].strip()
                        parts = time_str.split(':')
                        if len(parts) == 3:  # h:mm:ss
                            h, m, s = parts
                            elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                        elif len(parts) == 2:  # mm:ss
                            m, s = parts
                            elapsed_time = int(m) * 60 + float(s)
                        elif len(parts) == 1:  # just ss
                            elapsed_time = float(parts[0])
                    elif 'Maximum resident set size' in line:
                        max_memory_kb = int(line.split(':')[1].strip())
            return elapsed_time, max_memory_kb

        def get_replicate(filepath):
            match = re.search(r'rep(\d+)', filepath)
            return int(match.group(1)) if match else None

        def add_single_tool(file_list, tool, technology, results):
            """Parse individual tool timing files and append to results."""
            for filepath in file_list:
                replicate = get_replicate(filepath)
                elapsed, memory = parse_timing_file(filepath)
                if elapsed and memory and replicate:
                    results.append({
                        'tool': tool,
                        'technology': technology,
                        'replicate': replicate,
                        'elapsed_seconds': elapsed,
                        'max_memory_gb': memory / (1024 * 1024)
                    })

        def add_combined_tool(inq_files, filter_files, tool_files, tool_label, technology, results):
            """Sum timing across three steps (inquiSTR + filter + tool) and append to results."""
            for inq_file, filt_file, tool_file in zip(inq_files, filter_files, tool_files):
                replicate = get_replicate(inq_file)
                inq_elapsed, inq_memory = parse_timing_file(inq_file)
                filt_elapsed, filt_memory = parse_timing_file(filt_file)
                tool_elapsed, tool_memory = parse_timing_file(tool_file)
                if all([inq_elapsed, filt_elapsed, tool_elapsed, replicate]):
                    results.append({
                        'tool': tool_label,
                        'technology': technology,
                        'replicate': replicate,
                        'elapsed_seconds': inq_elapsed + filt_elapsed + tool_elapsed,
                        'max_memory_gb': max(inq_memory, filt_memory, tool_memory) / (1024 * 1024)
                    })

        results = []

        # PacBio single tools
        add_single_tool(input.pacbio_inquistr_time, 'inquiSTR', 'pacbio', results)
        add_single_tool(input.pacbio_trgt_time, 'TRGT', 'pacbio', results)
        add_single_tool(input.pacbio_longtr_time, 'LongTR', 'pacbio', results)

        # PacBio combined pipelines
        add_combined_tool(
            input.pacbio_inquistr_time, input.pacbio_filter_time, input.pacbio_trgt_filtered_time,
            'inquiSTR+TRGT', 'pacbio', results
        )
        add_combined_tool(
            input.pacbio_longtr_inquistr_time, input.pacbio_longtr_filter_time, input.pacbio_longtr_filtered_time,
            'inquiSTR+LongTR', 'pacbio', results
        )

        # ONT single tools
        add_single_tool(input.ont_inquistr_time, 'inquiSTR', 'ont', results)
        add_single_tool(input.ont_longtr_time, 'LongTR', 'ont', results)

        # ONT combined pipeline
        add_combined_tool(
            input.ont_longtr_inquistr_time, input.ont_longtr_filter_time, input.ont_longtr_filtered_time,
            'inquiSTR+LongTR', 'ont', results
        )

        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)

rule benchmark_inquistr_vs_trgt_pacbio:
    """
    Compare inquiSTR PacBio genotypes (--test) against TRGT PacBio genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/pacbio-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_vs_trgt_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule inquiSTR_adotto_requirespanning:
    """Run inquiSTR call on PacBio with --require-spanning to restrict genotypes to
    loci fully covered by spanning reads, for direct comparison against TRGT."""
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz"
    log:
        "logs/inquiSTR_adotto_requirespanning.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    shell:
        """
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend \
            --require-spanning 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule benchmark_inquistr_requirespanning_vs_trgt_pacbio:
    """Compare inquiSTR PacBio --require-spanning genotypes against TRGT genotypes.
    This isolates spanning-read-only calls to assess how much accuracy improves
    when soft-clipped genotypes are excluded."""
    input:
        test = "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz",
        truth = "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_requirespanning_vs_trgt_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule inquiSTR_ont_requirespanning:
    """Run inquiSTR call on ONT with --require-spanning to restrict genotypes to
    loci fully covered by spanning reads, for direct comparison against LongTR."""
    input:
        catalog = "adotto_LongTR.bed",
        ont = "ont.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/ont-inquistr-adotto-requirespanning.inq.gz"
    log:
        "logs/inquiSTR_ont_requirespanning.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    shell:
        """
        {params.inquiSTR} call {input.ont} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend \
            --require-spanning 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule benchmark_inquistr_requirespanning_vs_longtr_ont:
    """Compare inquiSTR ONT --require-spanning genotypes against LongTR genotypes.
    This isolates spanning-read-only calls to assess how much accuracy improves
    when soft-clipped genotypes are excluded."""
    input:
        test = "tool_comparison/ont-inquistr-adotto-requirespanning.inq.gz",
        truth = "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont.tsv",
        plot = "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont.html",
        diff_out = "tool_comparison/benchmark_inquistr_requirespanning_vs_longtr_ont_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_requirespanning_vs_longtr_ont.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_vs_longtr_pacbio:
    """
    Compare inquiSTR PacBio genotypes (--test) against LongTR PacBio genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/pacbio-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/pacbio-longtr-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_vs_longtr_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \\
            --test {input.test} \\
            --truth {input.truth} \\
            --plot {output.plot} \\
            --diff-out {output.diff_out} \\
            --max-locus {params.max_locus} \\
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_vs_longtr_ont:
    """
    Compare inquiSTR ONT genotypes (--test) against LongTR ONT genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/ont-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_longtr_ont.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_longtr_ont_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_vs_longtr_ont.log"
    shell:
        """
        {params.inquiSTR} benchmark \\
            --test {input.test} \\
            --truth {input.truth} \\
            --plot {output.plot} \\
            --diff-out {output.diff_out} \\
            --max-locus {params.max_locus} \\
            > {output.txt} 2> {log}
        """


rule inquiSTR_accuracy:
    """
    execute the inquiSTR benchmark analysis (named here accuracy because we already have a benchmark rule which is on time and memory)
    this command takes the inquiSTR call output and compares it to the adotto truth genotypes (in bed format)
    we will do this fo both ont and pacbio technologies, for a specific thread count and replicate (as that shouldn't matter for accuracy)
    """
    input:
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        genotypes = "benchmarking/data/{technology}_threads4_rep1.tsv",
        version = "inquiSTR_version.txt"
    output:
        txt = "benchmarking/accuracy_{technology}.tsv",
        plot = "benchmarking/accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS # limit to same length as used in genotyping
    log:
        "logs/inquiSTR_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """


rule trgt_accuracy:
    """
    Compare TRGT (PacBio) genotypes against the adotto/GIAB HG002 truth, mirroring
    inquiSTR_accuracy. Uses the rep1 TRGT calls converted to inquiSTR format as --test.
    """
    input:
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        genotypes = "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/trgt_accuracy_pacbio.tsv",
        plot = "tool_comparison/trgt_accuracy_pacbio.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS # limit to same length as used in genotyping
    log:
        "logs/trgt_accuracy_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """


rule longtr_accuracy:
    """
    Compare LongTR genotypes against the adotto/GIAB HG002 truth (both ONT and PacBio),
    mirroring inquiSTR_accuracy. Uses the rep1 LongTR calls converted to inquiSTR format as --test.
    """
    input:
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        genotypes = "tool_comparison/{technology}-longtr-adotto_rep1.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/longtr_accuracy_{technology}.tsv",
        plot = "tool_comparison/longtr_accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS # limit to same length as used in genotyping
    log:
        "logs/longtr_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """

rule stratify_accuracy_by_length:
    """Accuracy against the adotto/GIAB HG002 truth, stratified by truth allele length.

    Answers whether genotyping accuracy degrades for long repeat alleles, and reports the
    call rate per length bin so that dropout (loci the tool leaves uncalled) is separated
    from miscalling. The `inquiSTR-spanning` call sets are `inquiSTR call --require-spanning`,
    which excludes genotypes derived from soft-clipped reads, isolating that contribution.

    Two tolerances are scored: the fixed 3 bp used elsewhere in the paper, and one that
    scales with the allele (3 bp or 1% of the truth allele length, whichever is larger).
    A fixed 3 bp window demands 0.06% accuracy on a 5 kb allele but 6% on a 50 bp one, so
    the pair separates a genuine length effect from the tolerance simply getting stricter.
    """
    input:
        truth_bed = ADOTTO_TRUTH,
        catalog = "adotto_TRGT.bed.gz",
        genotypes = lambda wildcards: STRATIFIED_CALLSETS[(wildcards.technology, wildcards.tool)]
    output:
        tsv = "tool_comparison/stratified_accuracy/{technology}-{tool}.tsv",
        errors = "tool_comparison/stratified_errors/{technology}-{tool}.tsv.gz"
    wildcard_constraints:
        technology = "|".join(TECHNOLOGIES),
        # longest first so that e.g. "inquiSTR-spanning" is not truncated to "inquiSTR"
        tool = "|".join(sorted({tool for _, tool in STRATIFIED_CALLSETS}, key=lambda t: (-len(t), t)))
    params:
        tolerance = 3,
        relative_tolerance = 0.01,
        max_locus = MAX_LOCUS,
        edges = LENGTH_BIN_EDGES
    run:
        import numpy as np
        import pandas as pd

        truth = load_truth_alleles(input.truth_bed, params.max_locus)
        calls = load_call_alleles(input.genotypes)

        # Left join: truth loci absent from the call set (or called NaN) count as not called
        merged = truth.merge(calls, on=["chromosome", "begin"], how="left")
        merged["abs_error"] = (merged["call"] - merged["truth_allele"]).abs()

        # Only truth loci that are in the calling catalog were ever presented to a tool;
        # the rest cannot be called by anyone and must stay out of the call-rate denominator
        catalog_keys = load_catalog_keys(input.catalog)
        merged["targeted"] = pd.MultiIndex.from_arrays(
            [merged["chromosome"], merged["begin"]]
        ).isin(catalog_keys)

        # Tolerance that scales with the allele, floored at the fixed tolerance
        scaled_tolerance = np.maximum(
            params.tolerance, params.relative_tolerance * merged["truth_length"]
        )

        labels = length_bin_labels(params.edges)
        merged["length_bin"] = pd.cut(
            merged["truth_length"],
            bins=params.edges + [np.inf],
            right=False,
            labels=labels,
        )
        merged["within_scaled"] = merged["abs_error"] <= scaled_tolerance

        rows = []
        for label in labels:
            binned = merged[merged["length_bin"] == label]
            targeted = binned[binned["targeted"]]
            called = targeted.dropna(subset=["call"])
            n_truth = len(binned)
            n_targeted = len(targeted)
            n_called = len(called)
            errors = called["abs_error"]

            exact = int((errors == 0).sum())
            within_1bp = int((errors <= 1).sum())
            within_tolerance = int((errors <= params.tolerance).sum())
            within_scaled = int(called["within_scaled"].sum())

            def percent(count):
                return 100 * count / n_called if n_called else float("nan")

            rows.append({
                "technology": wildcards.technology,
                "tool": wildcards.tool,
                "length_bin": label,
                "n_truth": n_truth,
                "n_targeted": n_targeted,
                "catalog_coverage": n_targeted / n_truth if n_truth else float("nan"),
                "n_called": n_called,
                "call_rate": n_called / n_targeted if n_targeted else float("nan"),
                "exact": exact,
                "exact_percent": percent(exact),
                "within_1bp": within_1bp,
                "within_1bp_percent": percent(within_1bp),
                "within_tolerance": within_tolerance,
                "within_tolerance_percent": percent(within_tolerance),
                "within_scaled_tolerance": within_scaled,
                "within_scaled_tolerance_percent": percent(within_scaled),
                "tolerance_bp": params.tolerance,
                "relative_tolerance": params.relative_tolerance,
                "median_abs_error": errors.median() if n_called else float("nan"),
                "mean_abs_error": errors.mean() if n_called else float("nan"),
                # Quantiles of the absolute error. Unlike the mean these are not dragged by a
                # handful of extreme misses, and unlike the range of the error distribution
                # they do not grow simply because a bin holds more loci - the bins span
                # 142,330 down to 27 loci, so only an n-robust statistic can be compared
                # across them.
                "p90_abs_error": errors.quantile(0.90) if n_called else float("nan"),
                "p99_abs_error": errors.quantile(0.99) if n_called else float("nan"),
                "pearson_r": called["call"].corr(called["truth_allele"]) if n_called > 2 else float("nan"),
            })

        pd.DataFrame(rows).to_csv(output.tsv, sep="\t", index=False, float_format="%.4f")

        # Per-locus signed errors, for the error-distribution panel and as supplementary data.
        # Signed rather than absolute: under-calling a long allele and over-calling one are
        # different failures, and only the signed value distinguishes them. Every genotyped
        # locus is written; the figure subsamples the dense bins, this file does not.
        errors = merged[merged["targeted"] & merged["call"].notna()].copy()
        errors["error"] = errors["call"] - errors["truth_allele"]
        errors["technology"] = wildcards.technology
        errors["tool"] = wildcards.tool
        errors[[
            "technology", "tool", "chromosome", "begin", "length_bin",
            "truth_length", "truth_allele", "call", "error",
        ]].to_csv(output.errors, sep="\t", index=False)


rule combine_stratified_accuracy:
    """Concatenate the per-tool length-stratified accuracy tables into one long-format table."""
    input:
        expand("tool_comparison/stratified_accuracy/{technology}-{tool}.tsv",
               zip,
               technology=[technology for technology, _ in STRATIFIED_CALLSETS],
               tool=[tool for _, tool in STRATIFIED_CALLSETS])
    output:
        "tool_comparison/stratified_accuracy.tsv"
    run:
        import pandas as pd

        combined = pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True)
        combined.to_csv(output[0], sep="\t", index=False)


rule combine_stratified_errors:
    """Concatenate the per-tool, per-locus genotype errors into one long-format table."""
    input:
        expand("tool_comparison/stratified_errors/{technology}-{tool}.tsv.gz",
               zip,
               technology=[technology for technology, _ in STRATIFIED_CALLSETS],
               tool=[tool for _, tool in STRATIFIED_CALLSETS])
    output:
        "tool_comparison/stratified_errors.tsv.gz"
    run:
        import pandas as pd

        combined = pd.concat([pd.read_csv(path, sep="\t") for path in input], ignore_index=True)
        combined.to_csv(output[0], sep="\t", index=False)


rule stratified_accuracy:
    """Target rule for the accuracy-versus-truth-allele-length analysis."""
    input:
        "tool_comparison/stratified_accuracy.tsv",
        "tool_comparison/stratified_accuracy.html",
        "tool_comparison/stratified_errors.tsv.gz"


rule downsample_ont_coverage:
    """Subsample the deep HG002 ONT alignment to a target coverage.

    The fraction is the target divided by ONT_FULL_COVERAGE, so the levels are nominal depths
    based on that estimate rather than measured ones. The seed is fixed so reruns reproduce the
    same read subset. Outputs are temporary: the ten levels together hold ~5.5x the reads of the
    source file, and only the genotypes are needed afterwards.
    """
    input:
        cram = ONT_FULL_CRAM
    output:
        cram = temp("coverage/ont_{coverage}x.cram"),
        crai = temp("coverage/ont_{coverage}x.cram.crai")
    params:
        reference = reference,
        fraction = lambda wildcards: float(wildcards.coverage) / ONT_FULL_COVERAGE,
        seed = 42
    wildcard_constraints:
        coverage = r"\d+(\.\d+)?"
    threads: 4
    conda:
        "envs/minimap2.yml" # also includes samtools
    log:
        "logs/downsample_ont_{coverage}x.log"
    shell:
        """
        samtools view -C -T {params.reference} \
            --subsample {params.fraction} --subsample-seed {params.seed} \
            --threads {threads} -o {output.cram} {input.cram} 2> {log}
        samtools index {output.cram} 2>> {log}
        """


rule call_ont_coverage:
    """Genotype one downsampled ONT alignment, with the same settings as the main ONT run."""
    input:
        cram = "coverage/ont_{coverage}x.cram",
        crai = "coverage/ont_{coverage}x.cram.crai",
        catalog = "adotto_TRGT.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        inq = "coverage/ont_{coverage}x.inq.gz"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    log:
        "logs/call_ont_{coverage}x.log"
    shell:
        """
        {params.inquiSTR} call {input.cram} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule benchmark_ont_coverage:
    """Score one coverage level against the adotto/GIAB HG002 truth."""
    input:
        truth_bed = ADOTTO_TRUTH,
        genotypes = "coverage/ont_{coverage}x.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "coverage/accuracy_{coverage}x.tsv"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_ont_{coverage}x.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """


rule aggregate_coverage_accuracy:
    """Collect the per-coverage benchmark outputs into one table."""
    input:
        expand("coverage/accuracy_{coverage}x.tsv",
               coverage=[coverage_label(c) for c in COVERAGE_LEVELS])
    output:
        "coverage/accuracy_by_coverage.tsv"
    run:
        import pandas as pd

        rows = []
        for coverage, path in zip(COVERAGE_LEVELS, input):
            metrics = parse_accuracy_percentages(path)
            rows.append({
                "coverage": coverage,
                "loci_assessed": parse_loci_assessed(path),
                "exact_percent": metrics["exact"],
                "within_1bp_percent": metrics["within_1bp"],
                "within_3bp_percent": metrics["within_3bp"],
            })

        pd.DataFrame(rows).to_csv(output[0], sep="\t", index=False)


rule coverage_titration:
    """Target rule for the concordance-versus-coverage analysis."""
    input:
        "coverage/accuracy_by_coverage.tsv",
        "coverage/accuracy_by_coverage.html"


rule pathogenic:
    input:
        inquistr = expand("puretarget-calls/{sample}.inq", sample=puretarget_files),
        combined = "puretarget-calls/combined.tsv",
        heatmap = "puretarget-calls/heatmap.html",
        heatmap_png = "puretarget-calls/heatmap.png",

rule combine_puretarget:
    input:
        expand("puretarget-calls/{sample}.inq", sample=puretarget_files),
    output:
        "puretarget-calls/combined.tsv"
    log:
        "logs/combine_puretarget.log"
    threads: 4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input} > {output} 2> {log}
        """

rule genotype_puretarget:
    input:
        bam = "/home/AD/wdecoster/inquiSTR_paper/puretarget-data/{sample}.bam",
    output:
        inq = "puretarget-calls/{sample}.inq",
    log:
        "logs/genotype_puretarget_{sample}.log"
    params:
        reference = reference,
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} call --preset pathogenic --imbalance 0.1 --unphased {input.bam} > {output.inq} 2> {log}
        """


rule create_chr21_catalog_for_straglr:
    """Create a minimal STR catalog containing only chr21 loci."""
    input:
        catalog = "adotto_TRGT.bed.gz"
    output:
        catalog = "tool_comparison/straglr_chr21_catalog.bed"
    log:
        "logs/create_chr21_catalog_for_straglr.log"
    shell:
        """
        zcat {input.catalog} | awk 'BEGIN{{OFS="\t"}} $1=="chr21" {{
            n=split($4, fields, ";");
            motif="";
            for (i=1; i<=n; i++) {{
                if (fields[i] ~ /^STRUC=/) {{
                    motif=fields[i];
                    sub(/^STRUC=\\(/, "", motif);
                    sub(/\\)n.*$/, "", motif);
                    break;
                }}
            }}
            if (motif != "" && length(motif) >= 2 && length(motif) <= 50) print $1, $2, $3, motif;
        }}' > {output.catalog} 2> {log}
        """


rule run_straglr_chr21:
    """Run STRaglr on a chr21-only catalog for both PacBio and ONT."""
    input:
        catalog = "tool_comparison/straglr_chr21_catalog.bed",
        cram = "{technology}.cram"
    output:
        tsv = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.tsv",
        bed = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed",
        timing = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.time"
    log:
        "logs/straglr_chr21_{technology}_rep{replicate}.log"
    params:
        straglr = STRAGLR,
        reference = reference
    resources:
        benchmark_slot=1
    threads:
        4
    conda:
        "envs/straglr.yml"
    shell:
        """
        prefix=$(echo {output.tsv} | sed 's/\\.tsv$//')
        /usr/bin/time -v -o {output.timing} \
        {params.straglr} {input.cram} {params.reference} $prefix \
            --loci {input.catalog} \
            --genotype_in_size \
            --nprocs {threads} &> {log}
        """


rule run_inquistr_chr21:
    """Run inquiSTR on the same chr21-only catalog for both PacBio and ONT."""
    input:
        catalog = "tool_comparison/straglr_chr21_catalog.bed",
        cram = "{technology}.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/{technology}-inquistr-chr21_rep{replicate}.inq.gz",
        timing = "tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time"
    log:
        "logs/inquistr_chr21_{technology}_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    resources:
        benchmark_slot=1
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.cram} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule aggregate_straglr_chr21_results:
    input:
        straglr_time = expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        inquistr_time = expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES)
    output:
        "tool_comparison/straglr_chr21_results.tsv"
    run:
        import re
        import pandas as pd

        def parse_timing_file(filepath):
            elapsed_time = None
            max_memory_kb = None
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Elapsed (wall clock) time' in line:
                        time_str = line.split('): ', 1)[1].strip()
                        parts = time_str.split(':')
                        if len(parts) == 3:
                            h, m, s = parts
                            elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                        elif len(parts) == 2:
                            m, s = parts
                            elapsed_time = int(m) * 60 + float(s)
                        elif len(parts) == 1:
                            elapsed_time = float(parts[0])
                    elif 'Maximum resident set size' in line:
                        max_memory_kb = int(line.split(':')[1].strip())
            return elapsed_time, max_memory_kb

        def parse_metadata(filepath):
            match = re.search(r'tool_comparison/(ont|pacbio)-(straglr|inquistr)-chr21_rep(\d+)\.time$', filepath)
            if not match:
                return None, None, None
            technology, tool, replicate = match.groups()
            tool_name = 'Straglr' if tool == 'straglr' else 'inquiSTR'
            return technology, tool_name, int(replicate)

        results = []
        for timing_file in input:
            technology, tool_name, replicate = parse_metadata(str(timing_file))
            if technology is None:
                continue
            elapsed, memory = parse_timing_file(timing_file)
            if elapsed is None or memory is None:
                continue
            results.append({
                'technology': technology,
                'tool': tool_name,
                'replicate': replicate,
                'elapsed_seconds': elapsed,
                'max_memory_gb': memory / (1024 * 1024)
            })

        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)

rule convert_straglr_to_inquistr_format:
    """Convert STRaglr BED output to inquiSTR format for benchmarking.
    Removes header lines (starting with '#') and keeps columns 1-5 and 8.
    Adds inquiSTR-compatible header.
    """
    input:
        bed = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed"
    output:
        inq = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.inq.gz"
    log:
        "logs/convert_straglr_to_inquistr_{technology}_rep{replicate}.log"
    shell:
        """
           awk 'BEGIN {{OFS="\t"; print "chromosome", "begin", "end", "info", "straglr_H1", "straglr_H2"}} \
               $0 !~ /^#/ {{print $1, $2, $3, $4, $5, $8}}' {input.bed} \
             | gzip > {output.inq} 2> {log}
        """


rule benchmark_straglr_chr21_accuracy:
    """Benchmark STRaglr rep1 genotypes (rep1 for both ont and pacbio) against adotto truth."""
    input:
        test_inq = "tool_comparison/{technology}-straglr-chr21_rep1.inq.gz",
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/straglr_chr21_accuracy_{technology}.tsv",
        plot = "tool_comparison/straglr_chr21_accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_straglr_chr21_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test_inq} \
            --truth {input.truth_bed} \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_chr21_accuracy:
    """Benchmark inquiSTR rep1 genotypes (rep1 for both ont and pacbio) against adotto truth."""
    input:
        test_inq = "tool_comparison/{technology}-inquistr-chr21_rep1.inq.gz",
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/inquistr_chr21_accuracy_{technology}.tsv",
        plot = "tool_comparison/inquistr_chr21_accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_chr21_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test_inq} \
            --truth {input.truth_bed} \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


# ---------------------------------------------------------------------------------------------
# medaka tandem (ONT only), kept separate from the headline comparison for now
#
# Requested by a reviewer as a stronger ONT alternative to LongTR. Handled the same way STRaglr
# was: its own run, its own aggregation and its own plots, so the metrics can be inspected before
# deciding whether medaka tandem joins the main figures or stays a separate comparison.
#
# Two decisions worth revisiting, both flagged because they affect how favourably medaka is
# judged in a comparison a reviewer asked for:
#
#   1. MEDAKA_MODEL is empty, so medaka uses its default consensus model. If the HG002 ONT reads
#      were basecalled with a different model, medaka's consensus - and therefore its accuracy -
#      is understated. Set MEDAKA_MODEL to the matching model, or add `--auto_model consensus`
#      pointed at the CRAM, once the basecaller version is known.
#   2. medaka skips repeats whose estimated allele length exceeds 10 kb unless
#      `--process_large_regions` is given, listing them in `skipped_large.bed`. inquiSTR's
#      --max-locus instead filters on *reference span*, so medaka will drop some loci inquiSTR
#      calls. The default is kept because the alternative costs 14-23 GB of RAM, which would
#      distort the memory benchmark; the skipped count is recoverable from skipped_large.bed.
# ---------------------------------------------------------------------------------------------

rule run_medaka_tandem:
    """Genotype the ONT alignment with medaka tandem, timed like the other callers.

    CLI is `medaka tandem <bam> <ref> <regions.bed> <sex> <outdir>`; the sample sex is a required
    positional and HG002 is male. medaka reads the bgzipped adotto catalog directly, using only
    its first three columns. Results land in `<outdir>/medaka_to_ref.TR.vcf`.

    Requires medaka 2.2.1 or newer, the release that added both CRAM input and gzipped BED
    regions to the tandem subcommand. Older versions cannot pass a reference to htslib and fail
    every region with `[E::cram_next_slice] Failure to decode slice` followed by
    `Retrieved too few reads (0 < 3)`, surfacing at the end as
    `Medaka failed to generate a consensus sequence for the input regions`.
    """
    input:
        catalog = "adotto_TRGT.bed.gz",
        ont = "ont.cram",
    output:
        vcf = "tool_comparison/medaka_rep{replicate}/medaka_to_ref.TR.vcf",
        timing = "tool_comparison/ont-medaka-adotto_rep{replicate}.time"
    log:
        "logs/medaka_tandem_rep{replicate}.log"
    params:
        medaka = MEDAKA,
        reference = reference,
        sex = "male",  # HG002
        model = f"--model {MEDAKA_MODEL}" if MEDAKA_MODEL else ""
    threads:
        4
    conda:
        "/home/AD/wdecoster/inquiSTR_paper/envs/medaka.yml"
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    shell:
        """
        outdir=$(dirname {output.vcf})
        rm -rf "$outdir"
        /usr/bin/time -v -o {output.timing} \
        {params.medaka} tandem {params.model} \
            --workers {threads} \
            --sample_name medaka \
            {input.ont} \
            {params.reference} \
            {input.catalog} \
            {params.sex} \
            "$outdir" &> {log}
        """


rule convert_medaka_to_inquistr_format:
    """Convert the medaka tandem VCF to inquiSTR format for benchmarking."""
    input:
        vcf = "tool_comparison/medaka_rep{replicate}/medaka_to_ref.TR.vcf"
    output:
        inq = "tool_comparison/ont-medaka-adotto_rep{replicate}.inq.gz"
    run:
        medaka_vcf_to_inquistr(input.vcf, output.inq)


rule medaka_version:
    """Capture the medaka version alongside the other tool versions."""
    output:
        "medaka_version.txt"
    params:
        medaka = MEDAKA
    log:
        "logs/medaka_version.log"
    conda:
        "/home/AD/wdecoster/inquiSTR_paper/envs/medaka.yml"
    shell:
        """
        {params.medaka} --version > {output} 2> {log}
        """


rule medaka_accuracy_ont:
    """Compare medaka tandem genotypes against the adotto/GIAB HG002 truth.

    Mirrors longtr_accuracy so the two ONT callers are scored identically.
    """
    input:
        truth_bed = ADOTTO_TRUTH,
        genotypes = "tool_comparison/ont-medaka-adotto_rep1.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/medaka_accuracy_ont.tsv",
        plot = "tool_comparison/medaka_accuracy_ont.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/medaka_accuracy_ont.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_vs_medaka_ont:
    """Concordance between inquiSTR and medaka tandem genotypes on ONT.

    Mirrors benchmark_inquistr_vs_longtr_ont: medaka is passed as the truth set, so the
    percentages describe agreement between the two callers rather than accuracy.
    """
    input:
        test = "tool_comparison/ont-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/ont-medaka-adotto_rep1.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_medaka_ont.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_medaka_ont.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_medaka_ont_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_vs_medaka_ont.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule aggregate_medaka_results:
    """Runtime and peak memory for the three ONT callers on the full adotto catalog.

    inquiSTR and LongTR timings already exist from the main tool comparison, so this only adds
    medaka to the same table rather than re-running anything.
    """
    input:
        medaka = expand("tool_comparison/ont-medaka-adotto_rep{replicate}.time", replicate=REPLICATES),
        inquistr = expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        longtr = expand("tool_comparison/ont-longtr-adotto_rep{replicate}.time", replicate=REPLICATES)
    output:
        "tool_comparison/medaka_results.tsv"
    run:
        import re
        import pandas as pd

        def parse_timing_file(filepath):
            elapsed_time = None
            max_memory_kb = None
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Elapsed (wall clock) time' in line:
                        time_str = line.split('): ', 1)[1].strip()
                        parts = time_str.split(':')
                        if len(parts) == 3:
                            h, m, s = parts
                            elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                        elif len(parts) == 2:
                            m, s = parts
                            elapsed_time = int(m) * 60 + float(s)
                        elif len(parts) == 1:
                            elapsed_time = float(parts[0])
                    elif 'Maximum resident set size' in line:
                        max_memory_kb = int(line.split(':')[1].strip())
            return elapsed_time, max_memory_kb

        tool_names = {'medaka': 'medaka tandem', 'inquistr': 'inquiSTR', 'longtr': 'LongTR'}

        def parse_metadata(filepath):
            match = re.search(r'tool_comparison/ont-(medaka|inquistr|longtr)-adotto_rep(\d+)\.time$',
                              filepath)
            if not match:
                return None, None
            tool, replicate = match.groups()
            return tool_names[tool], int(replicate)

        results = []
        for timing_file in input:
            tool_name, replicate = parse_metadata(str(timing_file))
            if tool_name is None:
                continue
            elapsed, memory = parse_timing_file(timing_file)
            if elapsed is None or memory is None:
                continue
            results.append({
                'technology': 'ont',
                'tool': tool_name,
                'replicate': replicate,
                'elapsed_seconds': elapsed,
                'max_memory_gb': memory / (1024 * 1024),
            })

        df = pd.DataFrame(results).sort_values(['tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)


rule medaka:
    """Target rule for the medaka tandem comparison (ONT only)."""
    input:
        "tool_comparison/medaka_results.tsv",
        "tool_comparison/medaka_benchmark_plot.html",
        "tool_comparison/medaka_accuracy_ont.tsv",
        "tool_comparison/benchmark_inquistr_vs_medaka_ont.tsv",
        "medaka_version.txt"


rule metrics_summary:
    """Create a plain-text summary of speed differences and accuracy metrics."""
    input:
        runtimes = "tool_comparison/results.tsv",
        adotto_pacbio = "benchmarking/accuracy_pacbio.tsv",
        adotto_ont = "benchmarking/accuracy_ont.tsv",
        trgt_pacbio = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        longtr_pacbio = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        longtr_ont = "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        trgt_truth_pacbio = "tool_comparison/trgt_accuracy_pacbio.tsv",
        longtr_truth_pacbio = "tool_comparison/longtr_accuracy_pacbio.tsv",
        longtr_truth_ont = "tool_comparison/longtr_accuracy_ont.tsv"
    output:
        "tool_comparison/metrics_summary.txt"
    run:
        import re
        import pandas as pd

        def fold_summary(mean_df, faster_tool, slower_tool):
            techs = sorted(set(mean_df[mean_df["tool"] == faster_tool]["technology"]) &
                           set(mean_df[mean_df["tool"] == slower_tool]["technology"]))

            if not techs:
                raise ValueError(f"No overlapping technologies for {faster_tool} and {slower_tool}")

            folds = []
            for tech in techs:
                fast = mean_df[(mean_df["tool"] == faster_tool) & (mean_df["technology"] == tech)]["elapsed_seconds"].iloc[0]
                slow = mean_df[(mean_df["tool"] == slower_tool) & (mean_df["technology"] == tech)]["elapsed_seconds"].iloc[0]
                folds.append((tech, slow / fast))

            mean_fold = sum(fold for _, fold in folds) / len(folds)
            return mean_fold, folds

        runtimes = pd.read_csv(input.runtimes, sep="\t")
        mean_runtime = runtimes.groupby(["tool", "technology"], as_index=False)["elapsed_seconds"].mean()

        mean_trgt_fold, trgt_folds = fold_summary(mean_runtime, "inquiSTR", "TRGT")
        mean_longtr_fold, longtr_folds = fold_summary(mean_runtime, "inquiSTR", "LongTR")

        # Test set = the caller being evaluated; truth set = adotto/GIAB or another caller.
        accuracy_rows = [
            ("PacBio", "inquiSTR", "adotto/GIAB", parse_accuracy_percentages(input.adotto_pacbio)),
            ("ONT", "inquiSTR", "adotto/GIAB", parse_accuracy_percentages(input.adotto_ont)),
            ("PacBio", "TRGT", "adotto/GIAB", parse_accuracy_percentages(input.trgt_truth_pacbio)),
            ("PacBio", "LongTR", "adotto/GIAB", parse_accuracy_percentages(input.longtr_truth_pacbio)),
            ("ONT", "LongTR", "adotto/GIAB", parse_accuracy_percentages(input.longtr_truth_ont)),
            ("PacBio", "inquiSTR", "TRGT", parse_accuracy_percentages(input.trgt_pacbio)),
            ("PacBio", "inquiSTR", "LongTR", parse_accuracy_percentages(input.longtr_pacbio)),
            ("ONT", "inquiSTR", "LongTR", parse_accuracy_percentages(input.longtr_ont)),
        ]

        lines = []
        lines.append("Metrics Summary")
        lines.append("===============")
        lines.append("")
        lines.append("Mean speed differences")
        lines.append("----------------------")
        lines.append(f"inquiSTR vs TRGT: {mean_trgt_fold:.2f}x faster")
        for tech, fold in trgt_folds:
            lines.append(f"  {tech}: {fold:.2f}x")
        lines.append(f"inquiSTR vs LongTR: {mean_longtr_fold:.2f}x faster")
        for tech, fold in longtr_folds:
            lines.append(f"  {tech}: {fold:.2f}x")
        lines.append("")
        lines.extend(render_accuracy_table(accuracy_rows))

        with open(output[0], "w") as out_handle:
            out_handle.write("\n".join(lines) + "\n")


rule metrics_accuracy:
    """Create a plain-text summary of accuracy metrics without runtime aggregation."""
    input:
        adotto_pacbio = "benchmarking/accuracy_pacbio.tsv",
        adotto_ont = "benchmarking/accuracy_ont.tsv",
        trgt_pacbio = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        longtr_pacbio = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        longtr_ont = "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        trgt_truth_pacbio = "tool_comparison/trgt_accuracy_pacbio.tsv",
        longtr_truth_pacbio = "tool_comparison/longtr_accuracy_pacbio.tsv",
        longtr_truth_ont = "tool_comparison/longtr_accuracy_ont.tsv"
    output:
        "tool_comparison/metrics_accuracy.txt"
    run:
        # Test set = the caller being evaluated; truth set = adotto/GIAB or another caller.
        accuracy_rows = [
            ("PacBio", "inquiSTR", "adotto/GIAB", parse_accuracy_percentages(input.adotto_pacbio)),
            ("ONT", "inquiSTR", "adotto/GIAB", parse_accuracy_percentages(input.adotto_ont)),
            ("PacBio", "TRGT", "adotto/GIAB", parse_accuracy_percentages(input.trgt_truth_pacbio)),
            ("PacBio", "LongTR", "adotto/GIAB", parse_accuracy_percentages(input.longtr_truth_pacbio)),
            ("ONT", "LongTR", "adotto/GIAB", parse_accuracy_percentages(input.longtr_truth_ont)),
            ("PacBio", "inquiSTR", "TRGT", parse_accuracy_percentages(input.trgt_pacbio)),
            ("PacBio", "inquiSTR", "LongTR", parse_accuracy_percentages(input.longtr_pacbio)),
            ("ONT", "inquiSTR", "LongTR", parse_accuracy_percentages(input.longtr_ont)),
        ]

        lines = ["Accuracy Summary", "================", ""]
        lines.extend(render_accuracy_table(accuracy_rows))

        with open(output[0], "w") as out_handle:
            out_handle.write("\n".join(lines) + "\n")
