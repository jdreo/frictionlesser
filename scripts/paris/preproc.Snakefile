import datetime

configfile: "config.yaml"

# NOW=datetime.date.today().isoformat()
# workdir: "expe_{name}_{date}".format(name=config["expe-name"], date=NOW)

SEEDS=list(range(0,config["runs"]))

FRICTIONLESSER=config["executable"]

rule all:
    input:
        "data/inter/ranks.tsv",
        "cache/trans.cache.dat",
        expand("cache/size_{size}.cache.dat", size=config["sizes"])

rule preprocessing:
    input:
        counts="data/input/2022_02_18_version_2_EOC_counts.npz",
        features="data/input/2022_02_18_version_2_EOC_features.csv",
        meta="data/input/2022_02_18_version_2_EOC_meta.csv"
    output:
        "data/inter/2022_02_18_version_2_EOC_counts.mara.hdf5"
    shell:
        "python3 preproc-mara__npz-to-hdf5.py {input.counts} {input.features} {input.meta}"

rule counts:
    input:
        "data/inter/2022_02_18_version_2_EOC_counts.mara.hdf5"
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
    log: "logs/save_cache_transcriptome.log"
    benchmark: "logs/save_cache_transcriptome.bench"
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
    log: "logs/save_cache_size-{size}.log"
    benchmark: "logs/save_cache_size-{size}.bench"
    resources:
        job_name=lambda wildcards: f"cache_z{wildcards.size}"
    shell:
        "{FRICTIONLESSER}"
        "  --ranks={input.ranks}"
        "  --cache-transcriptome={input.transcache}"
        "  --cache-size={output}"
        "  --cache-only"
        "  2> {log}"

