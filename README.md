# DrConverters: unified cotranscriptional folding output.

This repository contains scripts to convert the output of different
cotranscriptional folding simulators into the \*.drf file format.
This format is used for visualization by the JavaScript application [DrForna],
and it can also be plotted in various flavors by the Python script `DrPlotter`
from the [DrTransformer] package for cotranscriptional folding.

## Installation
**Note:** This repository provides **wrapper scripts** for other programs, 
the respective dependencies must be installed separately (see below)!

The sripts provided in this repository can be installed using:
```sh
pip install .
```

### Dependencies

#### Kinfold
[Kinfold] is part of the [ViennaRNA package], which can be installed via
bioconda. First, make sure bioconda is set up properly with:
```sh
  ~$ conda config --add channels defaults
  ~$ conda config --add channels bioconda
  ~$ conda config --add channels conda-forge
  ~$ conda config --set channel_priority strict
```
Second, install or update your ViennaRNA installation.
```sh
  ~$ conda install 'viennarna>=2.5.1'
```

#### Kinefold
[Kinefold] must be downloaded (follow the link [Kinefold]). It is important
to have the executable `kinefold_long_static` placed in the working directory
where `DrKinefold` is used.

### Testing
Test the functionality via:

```sh
DrKinfold --help
DrKinefold --help
```

## Cites
The program wrappers provided in this repository use published software. Use
them at your own risk, and do not forget to cite the original software:
- `DrKinfold` depends on [Kinfold], publised in [Flamm et al. (2000)]
- `DrKinefold` depends on [Kinefold], published in [Xayaphoummine et al. (2005)]

The \*.drf output format was primarily developed for [DrForna] visualization.
If you find this visualization helpful during your analysis, consider citing
[Tanasie et al. (2023)]

The program [DrPlotter] provides additional Python visualization options, it is
part of the [DrTransformer] Python package, which should be cited as [Badelt et
al. (2023)].


[//]: References
[ViennaRNA package]: <http://www.tbi.univie.ac.at/RNA>
[ViennaRNA github]: <https://www.github.com/ViennaRNA/ViennaRNA>
[DrForna]: <https://github.com/ViennaRNA/drforna>
[Kinfold]: <https://www.tbi.univie.ac.at/RNA/Kinfold.1.html>
[Kinefold]: <http://kinefold.curie.fr/download.html>
[DrTransformer]: <https://github.com/ViennaRNA/drtransformer>
[Flamm et al. (2000)]: <https://doi.org/10.1017/s1355838200992161>
[Xayaphoummine et al. (2005)]: <doi.org/10.1093/nar/gki447>
[Tanasie et al. (2023)]: <https://>
[Badelt et al. (2023)]: <https://doi.org/10.1093/bioinformatics/btad034>
