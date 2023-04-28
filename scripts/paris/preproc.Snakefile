import datetime

configfile: "config_preproc.yaml"

FRICTIONLESSER=config["executable"]

rule all:
    input:
        "data/inter/ranks.tsv",
        "cache/trans.cache.dat",
        expand("cache/size_{size}.cache.dat", size=config["sizes"]),
        "data/inter/paris+ranks.h5an.gz"

rule preproc_mara:
    input:
        task="preproc_mara.py",
        counts="data/input/counts.npz",
        features="data/input/features.csv",
        meta="data/input/meta.csv"
    output:
        "data/inter/paris.h5an"
    shell:
        "python3 {input.task} {input.counts} {input.features} {input.meta} {output}"

rule counts_csv:
    input:
        task="counts_to_csv.py",
        data="data/inter/paris.h5an"
    output:
        protected("data/inter/counts.csv")
    shell:
        "python3 {input.task} {input.data} > {output}"

rule ranks_tsv:
    input:
        "data/inter/counts.csv"
    output:
        protected("data/inter/ranks.tsv")
    shell:
        "{FRICTIONLESSER} --exprs={input} > {output}"

rule ranks:
    input:
        task="ranks_tsv_to_h5an.py",
        csv="data/inter/ranks.tsv",
        h5an="data/inter/paris.h5an"
    output:
        "data/inter/paris+ranks.h5an.gz"
    shell:
        "python3 {input.task} {input.csv} {input.h5an} {output}"

rule save_cache_transcriptome:
    input:
        "data/inter/ranks.tsv"
    output:
        protected("cache/trans.cache.dat")
    log: "data/inter/logs/save_cache_transcriptome.log"
    benchmark: "data/inter/logs/save_cache_transcriptome.bench"
    shell:
        "{FRICTIONLESSER}"
        "  --ranks={input}"
        "  --cache-transcriptome={output}"
        "  --cache-only"
        "  2> {log}"

rule save_cache_size:
    input:
        ranks="data/inter/ranks.tsv",
        transcache="cache/trans.cache.dat"
    output:
        protected("cache/size_{size}.cache.dat")
    wildcard_constraints:
        # Wildcard {size} should be numeric.
        size="\d+"
    log: "data/inter/logs/save_cache_size-{size}.log"
    benchmark: "data/inter/logs/save_cache_size-{size}.bench"
    resources:
        job_name=lambda wildcards: f"cache_z{wildcards.size}"
    shell:
        "{FRICTIONLESSER}"
        "  --ranks={input.ranks}"
        "  --cache-transcriptome={input.transcache}"
        "  --cache-size={output}"
        "  --cache-only"
        "  2> {log}"

