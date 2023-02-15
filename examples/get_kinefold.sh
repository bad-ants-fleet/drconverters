#!/bin/bash
wget http://kinefold.curie.fr/download/kinefold_long_static.tgz
tar -xzf kinefold_long_static.tgz
cp TEST/kinefold_long_static .
rm -rf TEST
rm -rf kinefold_long_static.tgz
