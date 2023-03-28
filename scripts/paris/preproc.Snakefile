import datetime

configfile: "config_preproc.yaml"

FRICTIONLESSER=config["executable"]

rule all:
    input:
        "data/inter/ranks.tsv",
        "cache/trans.cache.dat",
        expand("cache/size_{size}.cache.dat", size=config["sizes"])

rule preprocessing:
    input:
        counts="data/input/counts.npz",
        features="data/input/features.csv",
        meta="data/input/eta.csv"
    output:
        "data/inter/counts.mara.hdf5"
    shell:
        "python3 preproc-mara__npz-to-hdf5.py {input.counts} {input.features} {input.meta}"

rule counts:
    input:
        "data/inter/counts.mara.hdf5"
    output:
        "data/inter/counts.csv"
    shell:
        "python3 counts__hdf5-to-csv.py {input} > {output}"

rule ranks:
    input:
        "data/inter/counts.csv"
    output:
        "data/inter/ranks.tsv"
    shell:
        "{FRICTIONLESSER} --exprs={input} > {output}"

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

