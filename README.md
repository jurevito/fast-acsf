# 🧪 Fast ACSF Featurizer
High performance implementation of atom-centered symmetry functions featurizer. It featurizer ligand - receptor complex into a fixed size vector. It is heavly inspired by [AQDnet featurizer](https://github.com/koji11235/AQDnet).

## Dependencies
Only dependencies are `numpy` and `rdkit`. You also need to have C compiler installed to get shared library.

## Usage
1. Compile C source code using example command `make compile`.
2. Import and define featurizer in Python.
3. Pass it ligand atom coordinates and protein atom coordinates along with their respective atomic numbers.

## ⚡ Performance
Testing was done on `i7-9700K` processor with `32GB` of `3200MHz` RAM. I featurized `aa2ar` from [PDBBind](https://www.pdbbind-plus.org.cn/) dataset. Receptor has 4570 atoms, while ligand has 40. On average it takes 75ms with default settings to generate feature vector with size 11583.
