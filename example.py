from rdkit import Chem
from acsf.featurizer import Featurizer

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

for i, mol in enumerate(supp):

    # Protonate ligand pose and get atom information.
    mol = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    mol_coords = mol.GetConformer().GetPositions()
    mol_atom_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

    # Featurize ligand - protein complex. 
    res = featurizer.featurize(
        mol_coords, 
        mol_atom_nums, 
        protein_coords, 
        protein_atom_nums
    )
    
    break
