# Supplementary data for a paper on delta learning CCSD(T) water

This repository contains supplementary data supporting the findings of the paper:

Towards Routine Condensed Phase Simulations with Delta-Learned Coupled Cluster
Accuracy: Application to Liquid Water:
[https://doi.org/10.1021/acs.jctc.5c01377](https://doi.org/10.1021/acs.jctc.5c01377)

by Niamh O'Neill, Benjamin Xu Shi, William J.Baldwin, William C. Witt, G'abor Cs'anyi, Julian D. Gale, Angelos Michaelides and Christoph Schran

## Contents
* `models/`:
The final r^2SCAN baseline (`models/baseline/model`) and delta [r^2SCAN -> CCSD(T)] (`models/baseline/model`) models used in the main paper.
3 models formats are included, LAMMPS compatible GPU and CPU models and ASE compatible GPU model.

* `models/*/dataset`:
  The datasets for training the above models, including the MACE input training files. Atomic energies are also provided.

* `cluster-cutting-code/cut_water.py`:
  A small python script to cleave out water clusters for the delta model from periodic configurations.

* `inputs/`:
LAMMPS forcefield inputs for both Symmetrix and GPU LAMMPS




