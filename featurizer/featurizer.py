import ctypes
import numpy as np
import os

# Load shared _library and set argument types.
_lib = ctypes.CDLL('./featurizer/featurizer.dll')

class Atom(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_double),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double),
        ("atom_index", ctypes.c_int),
    ]

_lib.featurize.argtypes = [ctypes.POINTER(Atom), ctypes.c_int, ctypes.POINTER(Atom), ctypes.c_int]
_lib.featurize.restype = ctypes.POINTER(ctypes.c_double)

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
    atom_array = (Atom * len(coords))()
    for i, ((x, y, z), atom_num) in enumerate(zip(coords, atom_nums)):
        atom_array[i].x = x
        atom_array[i].y = y
        atom_array[i].z = z
        atom_array[i].atom_index = atom_to_index.get(atom_num, len(atom_to_index))

    return atom_array

def featurize(mol_coords: np.array, mol_atom_nums: list[int], protein_coords: np.array, protein_atom_nums: list[int]):
    mol_array = convert_to_c_atom(mol_coords, mol_atom_nums)
    protein_array = convert_to_c_atom(protein_coords, protein_atom_nums)

    result_ptr = _lib.featurize(mol_array, len(mol_coords), protein_array, len(protein_coords))

    # Convert the result pointer to a NumPy array
    result = np.ctypeslib.as_array(result_ptr, shape=(11583,)) # FIXME: Size is hardcoded.
    result = np.copy(result)

    # Free the memory allocated in C
    #_lib.free(result_ptr)

    return result