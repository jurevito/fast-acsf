# 🧪 Fast ACSF Featurizer
High performance implementation of atom-centered symmetry functions featurizer. It featurizes ligand-receptor complex into a fixed size vector. It is heavly inspired by [AQDnet featurizer](https://github.com/koji11235/AQDnet).

## Dependencies
Only dependencies are `numpy` and `rdkit`. You need to have C compiler installed to compile source code into shared library.

## Usage
1. Compile C source code using example command `make compile`.
2. Import and define featurizer in Python.
3. Pass it ligand atom coordinates and protein atom coordinates along with their respective atomic numbers.

```python
from acsf.featurizer import Featurizer

featurizer = Featurizer(
    elements=['H', 'C', 'O', 'N'],
    angular_cutoff=8
)
```

## ⚡ Performance
Testing was done on `i7-9700K` processor with `32GB` of `3200MHz` RAM. Featurization of `aa2ar` complex from [PDBBind](https://www.pdbbind-plus.org.cn/) dataset took on average 75ms per pose. Receptor has 4570 atoms, while ligand has 40. Default settings were used to generate feature vector with size 11583.
