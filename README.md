# DrConverters: unified cotranscriptional folding output.

This repository contains scripts to convert the output of different
cotranscriptional folding simulators into the \*.drf file format.
This format is used for visualization by the JavaScript application [DrForna].
Note that the [DrTransformer] package for cotranscriptional folding
natively supports the drf file format, and so this format can also be
plotted in various flavors by the Python script `DrPlotter`
from the [DrTransformer] package.

## Installation
The scripts provided in this repository can be installed using:
```sh
pip install .
```
**Note:** This repository only provides functions to convert output from other
programs, as well as **wrapper scripts** to call those programs and convert
output automatically. The respective dependencies must be installed separately (see below).

### Dependencies

#### Kinfold [[Flamm et al. (2000)]] & ViennaRNA [[Lorenz et al. (2011)]]
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

#### Kinefold [[Xayaphoummine et al. (2005)]]
[Kinefold] must be downloaded (follow the link [Kinefold]). It is important
to have the executable `kinefold_long_static` placed in the working directory
where `DrKinefold` is used.
DrKinefold **also depends on ViennaRNA (>=2.5.1)**, due to library functions that 
help with the conversion of Kinefold output to pseudoknotted dot-bracket strings.

### Testing
Test the functionality of wrapper scripts via:

```sh
DrKinfold --help
DrKinefold --help
```

## Contributing
Did you find a bug? Or do you want to provide support for a different
cotranscriptional folding software? Please fork the repository and submit
a pull request. 


## Cites
The program wrappers provided in this repository use published software. Use
them at your own risk -- ask if you have questions -- and do not forget to cite the original software:
- `DrKinfold` depends on [Kinfold], publised in [Flamm et al. (2000)]
- `DrKinefold` depends on [Kinefold], published in [Xayaphoummine et al. (2005)]

The \*.drf output format was primarily developed for [DrForna] visualization.
If you find this visualization helpful during your analysis, consider citing
[Tanasie et al. (2023)].

The program `DrPlotter` provides additional Python visualization options, it is
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
[Lorenz et al. (2011)]: <https://doi.org/10.1186/1748-7188-6-26>

