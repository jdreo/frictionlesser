
Frictionlesser being a small project, any contribution are welcomed in any form.
Just be bold and either post an issue, a pull request or just send a message to <johann.dreo@pasteur.fr>.


Architecture
============

Frictionlesser implements:

- a local search algorithm,
- manipulating swaps over a binary partition of a subset of human genes,
- evaluated by an objective function, using RNA single-cell sequencing data.

The software comes in two main parts: the application and its library.
The application implements a command line executable, running the search algorithm.
Its entry point is the `app/frictionlesser.cpp` file.

The `datatester.cpp` binary basically re-implements the data checks
and can be used to double-check if input data file are correct, without running any algorithm.

The library holds the data structures, the objective function and some common functions.
Its entry point is the `include/frictionless` headers directory, along with the
`src/` implementation directory.


Terms
-----

Paradiseo and Frictionlesser sometimes differ on how they call this or that.
To help modularizing the code, all those terms may be used for different objects.

Nonetheless, they more or less points to similar concepts:
- solution, individual ≈ signature,
- objective function ≈ score, quality,
- search algorithm ≈ [meta]-heuristic, optimization algorithm, evolutionary algorithm,


Paradiseo
---------

The source code heavily rely on the [Paradiseo](https://github.com/nojhan/paradiseo)
framework for everything related to search algorithmics:

- the local search is implemented with the Paradiseo-MO module,
  which allows for easily modifying and extending the algorithm by just combining operators.
- the binary partition data structure and the corresponding "swap" neighborhood
  follows the (very light) Paradiseo conventions.
  (They are actually designed to be ultimately a part of Paradiseo)
- The objective function inherits from the (quite light) Paradiseo interface,
  which allows to be easily plugged into other software, thanks to Paradiseo's tooling.

The idea behind Paradiseo is to modularize optimization/search algorithms.
As such, it may be difficult to follow, as each module has its own set of interfaces ("operators"),
and various implementations are available.

The main design pattern may be hard to graps at first if you are not fluent in object programming.
See the ["20 years" preprint](https://arxiv.org/pdf/2105.00420.pdf) for a high-level view on it.
You can also look at the
[algopattern project](https://github.com/nojhan/algopattern/blob/master/cpp/strategy.cpp)
which is a gentle introduction to this kind of design pattern (albeit for another kind of algorithm).


Search Algorithm
----------------

The algorithm is only an assembling of Paradiseo-MO components.
It is thus completely implemented in a few lines, near the end of the `app/frictionlesser.cpp` file.
Most of the code is actually managing various way to log its execution.

If you want to have a look at the algorithm itself, you need to browse Paradiseo's code.

- The code of the [`moRandomBestHC` class](https://github.com/nojhan/paradiseo/blob/master/mo/src/algo/moRandomBestHC.h)
  is just a wrapper, actually pre-assembling an "explorer" for you.
- The [`moLocalSearch` class](https://github.com/nojhan/paradiseo/blob/master/mo/src/algo/moLocalSearch.h)
  contains some actual code from which you can follow the important operators.


Objective Function
------------------

The objective function is the high-level interface that computes a signature's quality.

The entry point for the objective function is the file `include/frictionless/eval.h`.

The objective function follows Paradiseo-MO's architucture for partial evaluations.
This allows to drastically reduce the amount of computations when evaluating a solution
that is just a gene swap away from another.

The entry points are:

- The `frictionless::EvalFull` implements the full evaluation of a completely new solution.
- The `frictionless::EvalSwap` implements a partial evaluation for swap neighborhoods.

These two classes heavily rely on the `FriedmanScore` class (`score.h`),
which computes the main statistic, and computes the data cache that allows
the partial evaluation.

The `FriedmanScore` itself relies on a `Transcriptome` (`transcriptome.h`),
which holds the input RNA expression data, along with various accessors onto it.

The score is computed for a given binary partition of the genes space,
which is held by the `Signature` class (see below).

Note that the name of members in the `FriedmanScore` follows the notation used
in the [Frictionlesser technical report](https://www.overleaf.com/project/6166fe78f282a9f39c869372).


Signature
---------

A solution to the problem is called a `Signature` (`signature.h`),
which is actually a `moBinaryPartition` (`moBinaryPartition.h`).
The binary partition is just a set of "selected" genes
(and its counterpart, a set of "rejected" genes).

It is coupled with a "Fitness", which is the slang term in Paradiseo for
"objective function value of a solution to the problem".

In Frictionlesser, the Fitness of a Signature is a `Score` (`signature.h`).
This `Score` essentially holds the cache allowing the partial evaluation.
It also hold the score value (a scalar), and the atomic score values by samples;
see the `ScoreDetails` class (`signature.h`).


Cache
-----

The cache system are the low-level data structure that are to be updated
when the score of a signature is updated after some change.
It is structured in three layers, depending on what is changing
when encountering a rew signature.

In the current setup, only the swap cache is supposed to be used during the search.
The two other caches are involved during data load,
and are managed by the high-level application (see `frictionlesser.cpp`).

All the details related to the cache system are in `cache.h`:
- `CacheTranscriptome`, for Friedman score's intermediate results that are tied to a given *transcriptome*,
- `CacheSize`, for results that are tied to a given *signature size*,
- `CacheSwap`, for results involved in *swaping two genes*.


Neighborhood
------------

The neighborhood describes how to "move" from one signature to another.

In Paradiseo-MO, this concept is at the core of the modularization,
and may be difficult to fully grasp at first.
You may first read the [Paradiseo-MO preprint](https://inria.hal.science/hal-00665421/)
to get an introduction.

A Paradiseo-MO "Neighbor" is not just another `Signature`,
but it implements *how* to move from one signature to another.
In `moBinaryPartitionSwapNeighbor` (`moBinaryPartitionSwapNeighbor.h`),
it stores a couple of genes: one to be selected, the other to be rejected,
hence modelling a swap that can be applied on a `Signature`.

The `moBinaryPartitionSwapNeighborhood` class (`moBinaryPartitionSwapNeighborhood.h`)
implements a way to *enumerate* all the possible neighbors of a given signature.
It actually generates *neighbors* and not *solutions*.


Other
-----

Frictionlesser uses the [clutchlog project](https://nojhan.github.io/clutchlog/)
for having nice, colored, logs that shows the log location.
Its configuration is set in `log.cpp`.

Frictionlesser also uses the [exceptions project](https://github.com/nojhan/exceptions)
for having clean exception classes declarations, holding the errors location.

The `frictionless.h` file holds some convenience functions.

The `src/pgamma.cpp` file is borrowed from the [R project](https://www.r-project.org).



Licensing
=========

TL;DR: *Frictionlesser is available under the AGPL v3*.

Frictionlesser itself is distributed under the GNU Affero Public General License v3.0 license (AGPL).
It's source code is (so far) fully copyrighted to the Institut Pasteur,
except for the code of the `pgamma` function, which is borrowed from the R project (under GPL).

Frictionlesser compiles against the [Paradiseo](https://github.com/nojhan/paradiseo)
project code, which is distributed under the LGPL v2.0 (for its core)
and the CeCILL license v2.1 (for the MO module).

The CeCILL license is fully compatible with the GPL, and the AGPL is basically a GPL
with added clauses on using the software as a service over a network.
Hence, the most restrictive license applies, *which is the **AGPL v3**.*

This means that any derivative work should be licensed under the same term,
which basically guarantee that you will always be able to get access to
the source code, whatever the setting in which you use this software.
