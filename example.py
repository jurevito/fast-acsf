from rdkit import Chem
from acsf.feat import Featurizer
import numpy as np
import pandas as pd

POSES_PATH = './data/poses.sdf'
PROTEIN_PATH = './data/protein.pdb'

# Read ligand and protein file.
supp = Chem.SDMolSupplier(POSES_PATH)
protein = Chem.MolFromPDBFile(PROTEIN_PATH)

# Protonate protein and get atom information.
protein = Chem.AddHs(protein, addCoords=True, addResidueInfo=True)
protein_coords = protein.GetConformer().GetPositions()
protein_atom_nums = [atom.GetAtomicNum() for atom in protein.GetAtoms()]

# Initialize featurizer.
featurizer = Featurizer(
    elements=['H', 'C', 'N', 'O', 'P', 'S', 'Cl', 'Zn']
)

mol_coords_list = []
mol_atoms_list = []

# Prepare ligand poses.
for i, mol in enumerate(supp):
    mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_coords = mol.GetConformer().GetPositions()
    mol_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    mol_coords_list.append(mol_coords)
    mol_atoms_list.append(mol_atoms)

mol_coords_matrix = np.stack(mol_coords_list, 0)
mol_atoms_matrix = np.stack(mol_atoms_list, 0)

# Featurize all ligand - protein poses. 
feats = featurizer.featurize(
    mol_coords_matrix, 
    mol_atoms_matrix, 
    protein_coords, 
    protein_atom_nums
)

# Create a pandas dataframe with the features.
df = pd.DataFrame(feats, columns=featurizer.labels)
print(df.head(10))
