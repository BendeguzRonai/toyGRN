# toyGRN
Extension of [toyLIFE](https://github.com/PabloCatalan/toylife) by Pablo Catalán.\
Using code from the [simulator](https://github.com/tamames/toyLIFE_simulator) by Pablo González.

This is an extended version of the framework used to study the evolution of [spatiotemporal gene expression patterns](https://royalsocietypublishing.org/rsif/article/17/167/20190843/36040).
Specifically, we have added:
  1. An extended genotype-phenotype map based on toyLIFE rules, that results in a richer phenotype space
  2. An individual-based modeling framework that allows the evolution of diverse populations in a 2-dimensional environment
  3. Large-scale mutations (deletion, cut-paste, copy-paste of genes)

# Usage
**Edit the output folder & file on lines 22-23 of main.cc**\
**Unzip 'PatternsData.zip'**

To print usage:
```
./main -help
```

# Modes
## Standard
Evolution experiment towards a specific target pattern with ancestor tracing and population snapshots.

## Multiple experiments
Fast mode for 20 independent evolution experiments towards a specific target pattern.

## Sample genomes
Creates N (user-defined) random genomes and outputs a file with 
1. Every pattern ID that was discovered & number of times it was discovered
2. Number of unique patterns discovered vs genomes generated (new phenotype discovery rate)

## Print genome
Input: genome(s) & target pattern\
Output: phenotype and fitness

## Compare genomes
Input: 2 genomes\
Output: If genomes have the same phenotype: "Same pattern", otherwise: both phenotypes




