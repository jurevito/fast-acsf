import ctypes
import numpy as np
import itertools
from rdkit import Chem
import os

# Load shared library.
current_dir = os.path.dirname(os.path.abspath(__file__))
lib_name = 'featurizer.dll' if os.name == 'nt' else 'featurizer.so'
lib_path = os.path.join(current_dir, lib_name)
_lib = ctypes.CDLL(lib_path)

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
        ("num_mols", ctypes.c_int),
    ]

class Result(ctypes.Structure):
    _fields_ = [
        ("features", ctypes.POINTER(ctypes.c_double)),
        ("num_rows", ctypes.c_int),
        ("num_cols", ctypes.c_int),
    ]

_lib.featurize.argtypes = [ctypes.POINTER(ctypes.POINTER(Atom)), ctypes.c_int, ctypes.POINTER(Atom), ctypes.c_int, Config]
_lib.featurize.restype = Result

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
        self.atomic_num_to_index = dict([(symbol_to_num[sym], i) for i, sym in enumerate(self.elements)])

        # Generate labels.
        elem_labels = self.elements + ['DU']
        rs_radial = np.arange(0.5, radial_cutoff, radial_step).tolist()
        rs_angular = np.arange(0.5, angular_cutoff, angular_step).tolist()
        theta_list = np.linspace(0, np.pi*2, num_theta, endpoint=False, dtype='float32').tolist()

        radial_pairs = ['_'.join(comb) for comb in itertools.product(elem_labels, repeat=2)]
        radial_labels = [f'{elems}_{rs}' for elems, rs in itertools.product(radial_pairs, range(len(rs_radial)))]

        elem_triplets = [f"{e_j}_{e_i}_{e_k}" for (e_i, (e_j, e_k)) in 
            list(itertools.product(
                elem_labels,
                list(itertools.combinations_with_replacement(elem_labels, 2))
            ))
        ]
        angular_labels = [f"{elem_triplet}_{theta}_{rs}" 
            for elem_triplet in elem_triplets
                for theta in range(len(theta_list))
                    for rs in range(len(rs_angular))]
        
        self.labels = radial_labels + angular_labels

    def __convert_to_c_atom(self, coords: np.array, atom_nums: list[str]):
        atom_array = (Atom * len(coords))()
        for i, ((x, y, z), atom_num) in enumerate(zip(coords, atom_nums)):
            atom_array[i].x = x
            atom_array[i].y = y
            atom_array[i].z = z
            atom_array[i].atom_index = self.atomic_num_to_index.get(atom_num, len(self.elements))

        return atom_array
    
    def __convert_to_c_atom_array(self, coords_matrix: np.array, atom_num_matrix: list[str]):
        atom_arrays = (ctypes.POINTER(Atom) * len(coords_matrix))()

        for i in range(len(coords_matrix)):
            atom_array = self.__convert_to_c_atom(coords_matrix[i], atom_num_matrix[i])
            atom_arrays[i] = ctypes.cast(atom_array, ctypes.POINTER(Atom))

        return atom_arrays
    
    def __setup_config(self, num_mols: int):
        config = Config()
        config.radial_cutoff = self.radial_cutoff
        config.angular_cutoff = self.angular_cutoff
        config.radial_step = self.radial_step
        config.angular_step = self.angular_step
        config.num_theta = self.num_theta
        config.num_elems = len(self.elements) + 1
        config.num_mols = num_mols

        return config

    def featurize(self, mol_coords_matrix: np.array, mol_atom_num_matrix: np.array, protein_coords: np.array, protein_atom_nums: np.array) -> np.array:
        if len(mol_coords_matrix) == 0:
            raise ValueError("Molecule coordinates matrix is empty.")

        mol_arrays = self.__convert_to_c_atom_array(mol_coords_matrix, mol_atom_num_matrix)
        protein_array = self.__convert_to_c_atom(protein_coords, protein_atom_nums)
        config = self.__setup_config(len(mol_coords_matrix))

        result = _lib.featurize(mol_arrays, len(mol_coords_matrix[0]), protein_array, len(protein_coords), config)
        res = np.ctypeslib.as_array(result.features, shape=(result.num_rows, result.num_cols))
        res = np.copy(res)

        _lib.free_features(result.features)
        return res
