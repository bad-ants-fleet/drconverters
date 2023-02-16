# DrConverters: Example input for DrForna visualization

This repository contains examples to convert the output of different
cotranscriptional folding simulators into the \*.drf visualization file format.

## Installation
We provide two scripts:
    - `DrKinfold`: A Python wrapper for Kinfold.
    - `DrKinefold`: A Python wrapper for Kinefold.
Both scripts rely on 

We provide a few example sequences:


## Cites
The program wrappers provided in this repository use published software. Use them
at your own risk, and do not forget to cite the original software:
 - `DrKinfold` depends on [Kinfold], publised in [Flamm et al. (2001)]
 - `DrKinefold` depends on [Kinefold], published in [Isambert (2001)]

The \*.drf output format was primarily developed for [DrForna] visualization.
If you find this visualization helpful during your analysis, consider citing
[Tanasie et al. (2023)]

The program [DrPlotter] provides additional Python visualization options, it is
part of the [DrTransformer] Python package, which should be cited as [Badelt et
al. (2023)].


## Installation
All conversion scripts require the Python bindings from the [ViennaRNA] package,
which are most easily obtained when installing the ViennaRNA package via bioconda.

```sh
> conda config --add channels defaults
> conda config --add channels bioconda
> conda config --add channels conda-forge
> conda config --set channel_priority strict
> conda install viennarna
```

## Kinfold -- A stochastic model using single base-pair transitions (no pseudoknots)
[Kinefold]

## Kinefold -- A stochastic model using helix transitions (including pseudoknots)

## DrTransformer -- A heuristic model based on base-pair transitions (no pseudoknots)

[//]: References
[ViennaRNA]: <http://www.tbi.univie.ac.at/RNA>
[DrForna]: <https://github.com/ViennaRNA/drforna>
[Kinfold]: <https://github.com/ViennaRNA/drforna>
[Kinefold]: <https://github.com/ViennaRNA/drforna>
[DrTransformer]: <https://github.com/ViennaRNA/drforna>
