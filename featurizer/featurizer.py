import ctypes
import numpy as np
from rdkit import Chem

# Load shared _library and set argument types. # FIXME: load also on linux.
_lib = ctypes.CDLL('./featurizer/featurizer.dll')

class Atom(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_double),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double),
        ("atom_index", ctypes.c_int),
    ]

class Config(ctypes.Structure):
    _fields_ = [
        ("radial_cutoff", ctypes.c_double),
        ("angular_cutoff", ctypes.c_double),
        ("radial_step", ctypes.c_double),
        ("angular_step", ctypes.c_double),
        ("num_theta", ctypes.c_int),
        ("num_elems", ctypes.c_int),
    ]

class Result(ctypes.Structure):
    _fields_ = [
        ("features", ctypes.POINTER(ctypes.c_double)),
        ("size", ctypes.c_int)
    ]

_lib.featurize.argtypes = [ctypes.POINTER(Atom), ctypes.c_int, ctypes.POINTER(Atom), ctypes.c_int, Config]
_lib.featurize.restype = ctypes.POINTER(ctypes.c_double)

class Featurizer:
    def __init__(self, 
            radial_cutoff: float = 12,
            angular_cutoff: float = 6,
            radial_step: float = 0.5,
            angular_step: float = 2.0,
            num_theta: int = 8,
            elements: list[str] = ['H', 'C', 'N', 'O', 'P', 'S', 'Cl', 'Zn']
        ):
        self.radial_cutoff = radial_cutoff
        self.radial_step = radial_step
        self.angular_step = angular_step
        self.angular_cutoff = angular_cutoff
        self.num_theta = num_theta
        self.elements = elements

        # Convert element atomic numbers to index offsets.
        symbol_to_num = dict([(Chem.Atom(i).GetSymbol(), i) for i in range(1, 119)])
        self.atomic_num_to_index = dict([(symbol_to_num[sym], i) for i, sym in enumerate(elements)])

    def __convert_to_c_atom(self, coords: np.array, atom_nums: list[str]):
        atom_array = (Atom * len(coords))()
        for i, ((x, y, z), atom_num) in enumerate(zip(coords, atom_nums)):
            atom_array[i].x = x
            atom_array[i].y = y
            atom_array[i].z = z
            atom_array[i].atom_index = self.atomic_num_to_index.get(atom_num, len(self.elements))

        return atom_array
    
    def __setup_config(self):
        config = Config()
        config.radial_cutoff = self.radial_cutoff
        config.angular_cutoff = self.angular_cutoff
        config.radial_step = self.radial_step
        config.angular_step = self.angular_step
        config.num_theta = self.num_theta
        config.num_elems = len(self.elements) + 1

        return config

    def featurize(self, mol_coords: np.array, mol_atom_nums: list[int], protein_coords: np.array, protein_atom_nums: list[int]) -> np.array:
        mol_array = self.__convert_to_c_atom(mol_coords, mol_atom_nums)
        protein_array = self.__convert_to_c_atom(protein_coords, protein_atom_nums)
        config = self.__setup_config()

        result_ptr = _lib.featurize(mol_array, len(mol_coords), protein_array, len(protein_coords), config)

        result = np.ctypeslib.as_array(result_ptr, shape=(11583,)) # FIXME: Size is hardcoded.
        result = np.copy(result)

        # FIXME: free result vector memory in C.
        #_lib.free(result_ptr)

        return result
