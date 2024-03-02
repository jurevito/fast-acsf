from rdkit import Chem
from featurizer.featurizer import featurize
import numpy as np
import pandas as pd
import time

if __name__ == '__main__':
    POSES_PATH = './data/poses.sdf'
    PROTEIN_PATH = './data/protein.pdb'

    # Compare data to truth.
    df = pd.read_csv('./data/aa2ar_features.csv', index_col=0)
    
    supp = Chem.SDMolSupplier(POSES_PATH)
    protein = Chem.MolFromPDBFile(PROTEIN_PATH)

    protein = Chem.AddHs(protein, addCoords=True, addResidueInfo=True)
    protein_coords = protein.GetConformer().GetPositions()
    protein_atom_nums = [atom.GetAtomicNum() for atom in protein.GetAtoms()]

    row_example = df.iloc[0].values
    row_res = None

    for mol in supp:
        mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
        mol_coords = mol.GetConformer().GetPositions()
        mol_atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

        start_time = time.perf_counter()
        row_res = featurize(mol_coords, mol_atom_nums, protein_coords, protein_atom_nums)
        stop_time = time.perf_counter()
        print(f'Execution time: {(stop_time - start_time)*1000:.2f}ms')
        break
    
    abs_diff = np.abs(row_example - row_res)
    max_diff = np.max(abs_diff)
    max_diff_index = np.argmax(abs_diff)
    print(f'Max diff: {max_diff:.3f} at ({max_diff_index}) with {row_example[max_diff_index]} != {row_res[max_diff_index]}')