from rdkit import Chem
from featurizer.featurizer import featurize
import numpy as np
import time

if __name__ == '__main__':
    POSES_PATH = './data/poses.sdf'
    PROTEIN_PATH = './data/protein.pdb'
    
    supp = Chem.SDMolSupplier(POSES_PATH)
    protein = Chem.MolFromPDBFile(PROTEIN_PATH)

    protein = Chem.AddHs(protein, addCoords=True, addResidueInfo=True)
    protein_coords = protein.GetConformer().GetPositions()
    protein_atom_nums = [atom.GetAtomicNum() for atom in protein.GetAtoms()]

    for mol in supp:
        mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
        mol_coords = mol.GetConformer().GetPositions()
        mol_atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

        start_time = time.perf_counter()
        featurize(mol_coords, mol_atom_nums, protein_coords, protein_atom_nums)
        stop_time = time.perf_counter()
        print(f'Execution time: {(stop_time - start_time)*1000:.2f}ms')
        break