import pandas as pd

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"
# Load a TSV file into a DataFrame, containing URLs of 1000 Genomes cram files
df = pd.read_csv(file_path, sep='\t', usecols=['sample', 'hg38_path', 'source'])
df = df[df["source"] == "Gustafson"]
reference = "/home/AD/wdecoster/database/GRCh38s.fa"
inquiSTR = "/home/AD/wdecoster/repositories/inquiSTR/target/x86_64-unknown-linux-musl/release/inquiSTR"

rule all:
    input:
        "codis_relate.tsv",
        "codis_pca.html"

rule create_codis_manifest:
    output:
        manifest = "codis_manifest.tsv"
    log:
        "logs/create_codis_manifest.log"
    run:
        with open(output.manifest, 'w') as f:
            f.write("bam_path\tsample_name\n")
            for idx, row in df.iterrows():
                f.write(f"{row['hg38_path']}\t{row['sample']}\n")


rule genotype_codis:
    input:
        "codis_manifest.tsv"
    output:
        "codis_combined.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR
    threads: 8
    log:
        "logs/genotype_codis.log"
    shell:
        """
        {params.inquiSTR} batch {input} \
            --output {output} \
            --preset codis \
            --unphased \
            --threads {threads} \
            --reference {params.reference} \
            --resume \
            > {log} 2>&1
        """


rule codis_relate:
    input:
        "codis_combined.tsv"
    output:
        "codis_relate.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/codis_relate.log"
    shell:
        """
        {params.inquiSTR} relate \
            --combined {input} \
            --output {output} \
            --threads {threads} \
            > {log} 2>&1
        """


rule codis_pca:
    input:
        "codis_combined.tsv"
    output:
        "codis_pca.html"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/codis_pca.log"
    shell:
        """
         {params.inquiSTR} pca \
            --combined {input} \
            --output {output} \
            --threads {threads} \
            > {log} 2>&1
        """


