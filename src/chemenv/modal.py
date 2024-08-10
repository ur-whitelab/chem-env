from modal import App, web_endpoint
from modal import Image

rdkit_image = Image.debian_slim(python_version="3.12").pip_install(["rdkit", "numpy"])
with rdkit_image.imports():
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    import numpy as np


app = App()


@app.function(image=rdkit_image)
@web_endpoint()
def tanimoto(s1: str, s2: str) -> float:
    """Calculate the Tanimoto similarity of two SMILES strings."""
    try:
        mol1 = Chem.MolFromSmiles(s1)
        mol2 = Chem.MolFromSmiles(s2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except (TypeError, ValueError, AttributeError):
        return "Error: Not a valid SMILES string"


app.function(image=rdkit_image)
def get_number_of_topologically_distinct_atoms(smiles: str, atomic_number: int = 1):
    """Return the number of unique `element` environments based on environmental topology.
    This corresponds to the number of peaks one could maximally observe in an NMR spectrum.

    Args:
        smiles (str): SMILES string
        atomic_number (int, optional): Atomic number. Defaults to 1.

    Returns:
        int: Number of unique environments.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)

        if atomic_number == 1:
            # add hydrogen
            mol = Chem.AddHs(molecule)
        else:
            mol = molecule

        # Get unique canonical atom rankings
        atom_ranks = list(Chem.rdmolfiles.CanonicalRankAtoms(mol, breakTies=False))

        # Select the unique element environments
        atom_ranks = np.array(atom_ranks)

        # Atom indices
        atom_indices = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == atomic_number
        ]
        # Count them
        return len(set(atom_ranks[atom_indices]))
    except (TypeError, ValueError, AttributeError):
        return "Error: Not a valid SMILES string"
