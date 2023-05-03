# About 
The files here are examples to help users generate the drf file format with
different cotranscriptional folding software. The particluar simulation results
are not necessarily useful, as parameters have not been chosen to match
experimental data.

# Reproducing results:

The \*.drf output files in this directory were produced using the following commandline calls:

- SRPn_kinfold.drf: 
    ```sh
    cat sequences/SRPn.fa | DrKinfold -p 20 -n 5 --name SRPn_kinfold
    ```
- SRPn_kinefold.drf: 
    ```sh
    cat sequences/SRPn.fa | DrKinefold -p 100 --name SRPn_kinefold
    ```
- SRPn_drtransformer.drf: 
    ```sh
    cat sequences/SRPn.fa | DrTransformer --t-lin 5 --t-log 30 --t-end 1e5 --o-prune 0.001 --name SRPn_drtransformer
    ```
- grow_kinfold.drf: 
    ```sh
    cat sequences/grow.fa | DrKinfold -p 20 -n 5 --name grow_kinfold
    ```
- grow_kinefold.drf: 
    ```sh
    cat sequences/grow.fa | DrKinefold -p 100 --name grow_kinefold
    ```
- grow_kinefold.drf: 
    ```sh
    cat sequences/grow.fa | DrTransformer --t-lin 5 --t-log 30 --t-end 1e5 --o-prune 0.001 --name grow_drtransformer
    ```
    
## Dependencies
Remember, `kinefold_long_static` must be in this directory to use `DrKinefold`.

```sh
wget http://kinefold.curie.fr/download/kinefold_long_static.tgz
tar -xzf kinefold_long_static.tgz
cp TEST/kinefold_long_static .
```
