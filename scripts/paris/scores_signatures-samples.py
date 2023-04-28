import sys
import pandas as pd
import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 4)
    size = int(sys.argv[1])
    if size == 0:
        size = None
    fout_ssscores = sys.argv[2]
    fsignatures = sys.argv[3:]

    print("Load signatures from: ", fsignatures, file=sys.stderr, flush=True)
    genesets,genome,breaks,genesets_scores,genesets_samplescores = signatures.load(fsignatures, filter_size = size, with_scores=True)
    nsignatures = len(genesets_samplescores)
    assert(len(genesets_samplescores) > 0)
    assert(len(list(genesets_samplescores.values())[0]) > 0)
    nsamples = len(list(genesets_samplescores.values())[0])
    print("Found",nsignatures,"unique signatures", file=sys.stderr, flush=True)

    joined = {}
    for sign in genesets_samplescores:
        j = " ".join(sorted([g for g in sign]))
        print(j, file=sys.stderr, flush=True)
        joined[j] = genesets_samplescores[sign]

    dfscores = pd.DataFrame.from_dict(joined, orient='index', columns = list(range(nsamples)))
    dfscores.sort_index(ascending=True, inplace=True)
    dfscores.to_csv(fout_ssscores)

    print("Done", file=sys.stderr, flush=True)
