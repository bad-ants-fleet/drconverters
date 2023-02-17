# Reproducing results:

The \*.drf output files in this directory were produced using the following commandline calls:

- SRPn_kinfold.drf: 
    ```sh
    cat sequences/SRPn.fa | DrKinfold -p 5 -n 10 --name SRPn_kinfold
    ```
- SRPn_kinefold.drf: 
    ```sh
    cat sequences/SRPn.fa | DrKinefold -p 10 --name SRPn_kinefold
    ```
    
## Dependencies
Remember, `kinefold_long_static` must be in this directory to use `DrKinefold`.

```sh
wget http://kinefold.curie.fr/download/kinefold_long_static.tgz
tar -xzf kinefold_long_static.tgz
cp TEST/kinefold_long_static .
```
