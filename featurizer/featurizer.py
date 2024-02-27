import ctypes
import numpy as np
import os

# Load shared _library and set argument types.
_lib = ctypes.CDLL('./featurizer/featurizer.dll')

class Coord(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_double),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double),
        ("atom_index", ctypes.c_int),
    ]

_lib.square.argtypes = [ctypes.c_int]
_lib.square.restype = ctypes.c_int

_lib.featurize.argtypes = [ctypes.POINTER(Coord), ctypes.c_int]
_lib.featurize.restype = None

def square(x: int) -> int:
    global _lib
    result = _lib.square(x)
    return int(result)

"""
Num: 1, Sym: H
Num: 6, Sym: C
Num: 7, Sym: N
Num: 8, Sym: O
Num: 15, Sym: P
Num: 16, Sym: S
Num: 17, Sym: Cl
Num: 30, Sym: Zn
Other Du
"""
atom_to_index = dict(zip([1, 6, 7, 8, 15, 16, 17, 30], range(10)))

def convert_to_c_atom(coords: np.array, atom_nums: list[int]):
    coords_array = (Coord * len(coords))()
    for i, ((x, y, z), atom_num) in enumerate(zip(coords, atom_nums)):
        coords_array[i].x = x
        coords_array[i].y = y
        coords_array[i].z = z
        coords_array[i].atom_index = atom_to_index.get(atom_num, len(atom_to_index))

    return coords_array

def featurize(mol_coords: np.array, mol_atom_nums: list[int], protein_coords: np.array, protein_atom_nums: list[int]):
    mol_array = convert_to_c_atom(mol_coords, mol_atom_nums)
    protein_array = convert_to_c_atom(protein_coords, protein_atom_nums)

    _lib.featurize(mol_array, len(mol_coords), protein_array, len(protein_coords))