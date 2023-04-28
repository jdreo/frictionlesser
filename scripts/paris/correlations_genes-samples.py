import sys
import scipy
import numpy as np
import anndata as ad
import pandas as pd
from sortedcollections import OrderedSet

import signatures

if __name__ == "__main__":

    assert(len(sys.argv) >= 5)
    franks = sys.argv[1]
    size = int(sys.argv[2])
    if size == 0:
        size = None
    fout_gccorr = sys.argv[3]
    # fout_sccorr = sys.argv[4]
    # fout_sgcorr = sys.argv[5]
    # fplot_gcorr = sys.argv[4]
    # fplot_selfcorr = sys.argv[7]
    fsignatures = sys.argv[4:]

    # Precision for floating-point computation asserts.
    epsilon = 1e-6
    # epsilon = 1e-1

    ###########################################################################
    # LOAD
    ###########################################################################

    print("Load annotated ranks data from: ", franks, file=sys.stderr, flush=True)
    cells_allgenes = ad.read(franks)
    print(cells_allgenes, file=sys.stderr, flush=True)
    ncells = cells_allgenes.shape[0]
    ngenes_all = cells_allgenes.shape[1]

    # cells_allgenes structure:
    #
    #   ←  all genes  →
    # ┌────────────────┐
    # │ var:"id",[…]   │
    # └────────────────┘
    # ┏━━━━━━━━━━━━━━━━┓  ┌────────────────┐
    # ┃ X:             ┃  │ obs:           │  ↑
    # ┃ counts         ┃┓ │ "sample",[…]   │ cells
    # ┃                ┃┃ │                │  ↓
    # ┗━━━━━━━━━━━━━━━━┛┃┓└────────────────┘
    #  ┃layers["ranks"] ┃┃
    #  ┗━━━━━━━━━━━━━━━━┛┃
    #   ┃layers["zscore"]┃
    #   ┗━━━━━━━━━━━━━━━━┛
    #    ←  all genes  →

    print("Loaded", ncells,"cells and", ngenes_all,"genes", file=sys.stderr, flush=True)

    print("Load signatures from: ", fsignatures, file=sys.stderr, flush=True)
    unique_signatures, unique_scores, origins = signatures.load_unique_signatures(fsignatures, size)
    nsignatures = len(unique_signatures)
    genome = signatures.load_genome(fsignatures, filter_size = size)
    print("Found",nsignatures,"unique signatures", file=sys.stderr, flush=True)

    print("Gather samples...", file=sys.stderr, flush=True)
    samples = OrderedSet()
    for s in cells_allgenes.obs["sample"]:
        samples.add(s)
    print("Found",len(samples),"samples", file=sys.stderr, flush=True)
    nsamples = len(samples)
    assert(len(samples) > 0)

    ###########################################################################
    # GENES CORRELATIONS
    # > All genes
    ###########################################################################

    print("Compute cell-allgenes z-scores sample fractions...", file=sys.stderr, flush=True)

    # cells_allgenes structure:
    #
    #   ←  all genes  →
    # ┌────────────────┐
    # │ var:"id",[…]   │
    # └────────────────┘
    # ┏━━━━━━━━━━━━━━━━┓  ┌────────────────┐
    # ┃ X:             ┃  │ obs:           │  ↑
    # ┃ counts         ┃┓ │ "sample",[…]   │ cells
    # ┃                ┃┃ │                │  ↓
    # ┗━━━━━━━━━━━━━━━━┛┃┓└────────────────┘
    #  ┃layers["ranks"] ┃┃
    #  ┗━━━━━━━━━━━━━━━━┛┃
    #   ┃layers["zscore"]┃
    #   ┗━━━━━━━━━━━━━━━━┛
    #    ←  all genes  →
    # ┌────────────────┐
    # │ varp:          │┐
    # │ "correlations" ││  ↑
    # │                ││ all genes
    # │                ││  ↓
    # │                ││
    # └────────────────┘│
    #  │ "p-values"     │
    #  └────────────────┘
    #  ←  all genes  →

    ranks = cells_allgenes.layers["ranks"]
    cells_allgenes.layers["zscore"] = np.zeros(ranks.shape)
    for i,s in enumerate(samples):
        sample_cells = (cells_allgenes.obs["sample"] == s)
        print("\t", "Sample", i, ":", file=sys.stderr, flush=True)
        print("\t\t", np.sum(sample_cells), "sample_cells", file=sys.stderr, flush=True)
        assert(np.sum(sample_cells) > 1)

        # Z-score computation avoiding null standard deviations:
        #  1. compute stdev
        #  2. find index for when stdev != 0
        #  3. compute z-score
        #  4. other sample_cells are already zero since instantiation
        #
        #                              ←  all genes  →
        #         ╭─╮                  ┏━━━━━━━━━━━━━━━━┓layers["ranks"]
        #         │ │                ↑ ┃                ┃
        #         │ │                  ┃┌──────────────┐┃
        #         │ │                  ┃│              │┃ ╮
        #         │ │            cells ┃│   sranks     │┃ ├ sample
        #         │ │                  ┃│              │┃ ╯
        #         │ │                ↓ ┃└──────────────┘┃  ⮧
        #         │ │                  ┗━━━━━━━━━━━━━━━━┛
        #         │ │                    ⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣    ┌──────────────┐
        #         │ │             smean ┌──────────────┐   │ Ⓝ  Ⓝ ⓃⓃ Ⓝ   Ⓝ│
        #         │ │                   └──────────────┘ ⭢ │ Ⓝ z-scores  Ⓝ│
        #           │              sdev ┌──────────────┐   │ Ⓝ  Ⓝ ⓃⓃ Ⓝ   Ⓝ│
        #      loop │                   └──────────────┘ ⭢ └──────────────┘
        #      over │                    ⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣
        #   samples │           sdev!=0  ◌●◌◌●◌●●◌●◌◌◌●           ⭣
        #           │                    ⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣⭣     ┌──────────────┐ ┌──────────────┐
        #         ▲ │           sample_  ◌●◌◌●◌●●◌●◌◌◌●    │ Ⓝ  Ⓝ ⓃⓃ Ⓝ   Ⓝ│ │              │
        #         │ │           uplifted ◌●◌◌●◌●●◌●◌◌◌● ⭠  │ Ⓝ z-scores  Ⓝ│∨│     zeros    │
        #         │ │             _genes ◌●◌◌●◌●●◌●◌◌◌●    │ Ⓝ  Ⓝ ⓃⓃ Ⓝ   Ⓝ│ │              │
        #         │ │                          ⭣           └──────────────┘ └──────────────┘
        #         │ │                   ┌──────────────┐
        #         │ │                   │ 0  0 00 0   0│
        #         │ │                   │ 0 szscores  0│
        #         │ │                   │ 0  0 00 0   0│
        #         │ │                   └──────────────┘
        #         ╰─╯                          ┊
        #                              ┏━━━━━━━┊━━━━━━━━┓layers["zscore"]
        #                            ↑ ┃       ┊        ┃┓layers["ranks"]
        #                              ┃┌──────▽───────┐┃┃
        #                              ┃│              │┃┃
        #                        cells ┃│   zscores    │┃┃
        #                              ┃│              │┃┃
        #                            ↓ ┃└──────────────┘┃┃
        #                   ⮦          ┗━━━━━━━━━━━━━━━━┛┃
        #                               ┗━━━━━━━━━━━━━━━━┛
        #                                ←  all genes  →
        #   layers["zscore"].T        ×
        #          ┏━━━━━━━━━━━━━━━━┓  ┏━━━━━━━━━━━━━━━━┓varp["correlations"]
        #        ↑ ┃                ┃  ┃                ┃
        #          ┃                ┃  ┃                ┃
        #    all   ┃                ┃  ┃                ┃
        #    genes ┃     zscores    ┃  ┃  correlations  ┃
        #          ┃                ┃  ┃                ┃
        #        ↓ ┃                ┃  ┃                ┃
        #          ┗━━━━━━━━━━━━━━━━┛  ┗━━━━━━━━━━━━━━━━┛

        sranks = cells_allgenes.layers["ranks"][sample_cells,:]

        smean = np.mean(sranks, axis=0)
        assert(smean.shape[0] == ngenes_all)

        # Derivate standard-deviation from the computed mean,
        # to ensure we do not mess with pop VS sample statistics.
        sdev  = np.sqrt(np.sum(np.power(sranks - smean, 2), axis=0))
        assert(sdev.shape[0] == ngenes_all)

        uplifted_genes = (sdev != 0)
        print("\t\t", np.sum(uplifted_genes), "uplifted genes (with std-dev ≠ 0)", file=sys.stderr, flush=True)
        assert(np.sum(uplifted_genes) > 0)

        # Following sections assigned the exploded subarray of sample's cells × uplifted genes to the final array.
        # I would have like to do something like:
        # cells_allgenes.layers["zscore"][sample_cells,:][:,uplifted_genes] = (sranks[:,uplifted_genes] - smean[:,uplifted_genes]) / sdev[:,uplifted_genes] * 1/np.sqrt(nsamples)
        # but such a compound assignement does not occur (nor does it raise an error…).

        # The std dev may be null, so we must allow for division by zero.
        np_err_settings = np.geterr() # Save the previous config.
        # Allow for division by zero and invalid floating-point operation.
        np.seterr(divide='ignore', invalid='ignore')
        # Actually do (some) divisions by zero.
        zs = (sranks - smean) / sdev * 1/np.sqrt(nsamples)
        # print( np.mean(zs[:,uplifted_genes] * np.sqrt(nsamples), axis=1).T, file=sys.stderr, flush=True)
        # print(1-np.std(zs[:,uplifted_genes] * np.sqrt(nsamples), axis=1).T, file=sys.stderr, flush=True)
        # assert(all(  np.mean(zs[:,uplifted_genes] * np.sqrt(nsamples), axis=1) <= epsilon ))
        # assert(all( 1-np.std(zs[:,uplifted_genes] * np.sqrt(nsamples), axis=1) <= epsilon ))

        # Double check the old way:
        # g_zscores = np.zeros(sranks.shape)
        # print("g_zscores.shape:", g_zscores.shape, file=sys.stderr, flush=True)
        # for j,g in enumerate(cells_allgenes.var["id"]):
        #     # print("gene",j,g, file=sys.stderr, flush=True)
        #     granks = cells_allgenes.layers["ranks"][sample_cells, (cells_allgenes.var["id"] == g) ]
        #     # granks = sranks[:, (cells_allgenes.var["id"] == g) ]
        #     # print("cells ranks:",granks.shape,granks, file=sys.stderr, flush=True)
        #     # assert(len(granks.shape) == 1) # only one gene should match.
        #     gncells = granks.shape[0]
        #     # print("gncells:",gncells, file=sys.stderr, flush=True)
        #     assert(gncells > 0)
        #     gmean = np.mean(granks) # Numpy uses the population mean.
        #     assert(gmean == (gncells+1)/2)
        #     # Derivate standard-deviation from the computed mean,
        #     # to ensure we do not mess with pop VS sample statistics.
        #     gsdev = np.sqrt(np.sum(np.power(granks-gmean,2)))
        #     # print("\tgsdev:",gsdev, file=sys.stderr, flush=True)
        #     if gsdev != 0:
        #         # Standardized z-score across cells,
        #         # and pre-division for preparing the incoming dot product.
        #         g_zscores[:,j] = (granks-gmean)/gsdev * 1/np.sqrt(nsamples)

        #         # print("\tz-scores:",g_zscores[:,j], file=sys.stderr, flush=True)

        #         ssp2 = np.sum(np.power(g_zscores[:,j], 2)) * nsamples
        #         # print("\tssp2:", ssp2, file=sys.stderr, flush=True)
        #         assert(np.abs(1-ssp2) <= epsilon)
        #     else:
        #         g_zscores[:,j] = np.zeros(gncells)
        #         assert(np.sum(np.power(g_zscores[:,j], 2)) * nsamples == 0)
        #
        # np.testing.assert_array_equal(g_zscores, np.nan_to_num(zs,nan=0.0))

        assert(zs.shape == sranks.shape)
        # Assert divisions by zeros occured where expected.
        assert(np.isfinite(zs[:,uplifted_genes]).all())
        assert(np.logical_not(np.isfinite(zs[:,np.logical_not(uplifted_genes)])).all())
        # Get back the standard way of handling division by zero.
        np.seterr(**np_err_settings)

        # Duplicate the uplifted_genes index down to fit the sample's cells array shape.
        sample_uplifted_genes = np.tile(uplifted_genes, (zs.shape[0],1))
        assert(sample_uplifted_genes.shape == zs.shape)
        # If a cell is from an uplifted gene column, get its z-score,
        # else set it to zero (non uplifted columns should hold NaNs).
        # FIXME test if an in-place np.nan_to_num may be faster.
        szscores = np.where(sample_uplifted_genes, zs, np.zeros(zs.shape))
        assert(szscores.shape == zs.shape)
        # Assign the result back to the full array.
        cells_allgenes.layers["zscore"][sample_cells,:] = szscores

        if __debug__:
            # print(szscores, file=sys.stderr, flush=True)
            np.testing.assert_array_equal(szscores, cells_allgenes.layers["zscore"][sample_cells,:])

            sznull  = cells_allgenes.layers["zscore"][sample_cells,:][:,np.logical_not(uplifted_genes)]
            assert(all(np.sum(sznull, axis=1) == 0))

            # print(np.mean(szscores * np.sqrt(nsamples), axis=1), file=sys.stderr, flush=True)
            # assert(all(  np.mean(szscores * np.sqrt(nsamples), axis=1) <= epsilon ))
            # assert(all( 1-np.std(szscores * np.sqrt(nsamples), axis=1) <= epsilon ))

            # print(np.sqrt(np.sum(np.power(szscores[:,uplifted_genes], 2) * nsamples, axis=0)))
            # FIXME this should pass:
            assert(all( 1-np.sqrt(np.sum(np.power(szscores[:,uplifted_genes], 2) * nsamples, axis=0)) <= epsilon ))

    # Assert sum of scores for each gene.
    if __debug__:
        for g in genome:
            cgene = cells_allgenes.layers["zscore"][:,(cells_allgenes.var["id"] == g)]
            # print(np.sum(cgene, axis=1), file=sys.stderr, flush=True)
            assert(all(0 <= s <= 1 for s in np.sum(np.power(cgene,2), axis=1)))

    print("Compute pairwise average correlation matrix for all genes over samples...", file=sys.stderr, flush=True)
    # Correlation is the average of the sum of the product of z-scores.
    # We already prepared for the division by the number of samples, hence only the sum of products remains.
    # Clipping is necessary because of rounding errors.
    cells_allgenes.varp["correlations"] = np.clip(cells_allgenes.layers["zscore"].T @ cells_allgenes.layers["zscore"], -1,1)
    # if __debug__:
    #     for i,c in np.ndenumerate(cells_allgenes.varp["correlations"]):
    #         if not -1 <= c <= 1:
    #             print(i,c, file=sys.stderr, flush=True)
    #             assert(1-np.abs(c) <= epsilon)
    #

    print("Compute p-values for all genes over samples...", file=sys.stderr, flush=True)
    dist = scipy.stats.beta(ncells/2 -1, ncells/2 - 1, loc=-1, scale=2)
    cells_allgenes.varp["p-values"] = 2 * dist.cdf(-np.abs(cells_allgenes.varp["correlations"]))

    print("Update data into",franks, file=sys.stderr, flush=True)
    cells_allgenes.write(franks, compression="gzip")
    print(cells_allgenes, file=sys.stderr, flush=True)

    # FIXME genes frequency distribution


    ###########################################################################
    # GENES CORRELATIONS
    # > Signatures' genes
    ###########################################################################

    print("Compute cell-gene z-scores sample fractions...", file=sys.stderr, flush=True)

    # No genes selected.
    genome_ids = np.zeros(len(cells_allgenes.var["id"]), dtype=bool)
    for g in genome:
        # Add this gene.
        genome_ids |= (cells_allgenes.var["id"] == g) # FIXME try with "in genome"?
    # cells_sgenes.X = cells_allgenes.layers["zscore"][ :, genome_ids ]
    assert(np.sum(genome_ids) == len(genome))

    cells_sgenes = cells_allgenes[:,genome_ids]
    assert("zscore" in cells_allgenes.layers)
    assert("correlations" in cells_allgenes.varp)

    # cells_sgenes = ad.AnnData(
    #     np.zeros((ncells,ngenes)),
    #     cells_allgenes.obs, {"gene":np.asarray(genome)},
    #     dtype=cells_allgenes.layers["ranks"].dtype )

    # cells_sgenes structure:
    #
    #  ← genes in sign →
    # ┌────────────────┐
    # │ var:"gene"     │
    # └────────────────┘
    # ┏━━━━━━━━━━━━━━━━┓  ┌────────────────┐
    # ┃ X: z-score     ┃  │ obs:           │  ↑
    # ┃                ┃  │ "sample"       │ cells
    # ┃                ┃  │                │  ↓
    # ┗━━━━━━━━━━━━━━━━┛  └────────────────┘
    # ┌────────────────┐
    # │ varp:          │
    # │ "correlations" │  ↑
    # │                │ genes in sign
    # │                │  ↓
    # │                │
    # └────────────────┘
    #  ← genes in sign →
    #

    # # FIXME just extract this from the matrix with all genes.
    # genome_ids = np.zeros(len(cells_allgenes.var["id"]), dtype=bool)
    # for g in genome:
    #     genome_ids |= (cells_allgenes.var["id"] == g)
    # # cells_sgenes = cells_allgenes[ :, genome_ids ]
    # cells_sgenes.X = cells_allgenes.layers["zscore"][ :, genome_ids ]
    # assert(np.sum(genome_ids) == len(genome))

    # print("Compute pairwise average correlation matrix for signatures' genes over samples...", file=sys.stderr, flush=True)
    # # Correlation is the average of the sum of the product of z-scores.
    # # We already prepared for the division by the number of samples, hence only the sum of products remains.
    # # cells_sgenes.varp["correlations"] = np.clip(cells_sgenes.X.T @ cells_sgenes.X, -1,1)
    # cells_sgenes.varp["correlations"] = cells_sgenes.X.T @ cells_sgenes.X
    # if __debug__:
    #     for i in range(cells_sgenes.varp["correlations"].shape[0]):
    #         for j in range(cells_sgenes.varp["correlations"].shape[1]):
    #             # if not (-1 <= cells_sgenes.varp["correlations"][i,j] <= 1):
    #                 # print(i,j,cells_sgenes.varp["correlations"][i,j], file=sys.stderr, flush=True)
    #             assert(np.abs(cells_sgenes.varp["correlations"][i,j]) <= 1+epsilon)


    print("Save cell-gene z-score and gene correlations...", file=sys.stderr, flush=True)
    print(cells_sgenes, file=sys.stderr, flush=True)
    cells_sgenes.write(fout_gccorr, compression = "gzip")

    ###########################################################################
    # PLOT
    ###########################################################################

    # FIXME clustermap plot

    print("Done", file=sys.stderr, flush=True)
