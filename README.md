# Example input for DrForna visualization

This repository contains examples to produce the [DrForna] \*.drf visualization
file format from different cotranscriptional folding simulators.

### Installation
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
