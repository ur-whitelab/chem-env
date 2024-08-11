def tanimoto(s1: str, s2: str) -> float:
    """
    Calculate the Tanimoto similarity of two SMILES strings.

    Args:
        s1 (str): The first SMILES string.
        s2 (str): The second SMILES string.

    Returns:
        float: The Tanimoto similarity between the two SMILES strings.

    Raises:
        ValueError: If either of the input SMILES strings is invalid.

    This function calculates the Tanimoto similarity between two SMILES strings.
    It uses the RDKit library to calculate the Morgan fingerprints of the molecules.
    The Tanimoto similarity is then calculated using the DataStructs module of the RDKit library.

    Example:
        >>> tanimoto("CCO", "CC")
        0.143
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    try:
        mol1 = Chem.MolFromSmiles(s1)
        mol2 = Chem.MolFromSmiles(s2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except (TypeError, ValueError, AttributeError) as e:
        raise ValueError("Invalid SMILES strings") from e


def get_number_of_topologically_distinct_atoms(smiles: str, atomic_number: int = 1):
    """Return the number of unique `element` environments based on environmental topology.
    This corresponds to the number of peaks one could maximally observe in an NMR spectrum.
    Args:
        smiles (str): SMILES string
        atomic_number (int, optional): Atomic number. Defaults to 1.

    Returns:
        int: Number of unique environments.

    Raises:
        ValueError: If not a valid SMILES string

    Example:
        >>> get_number_of_topologically_distinct_atoms("CCO", 1)
        2

        >>> get_number_of_topologically_distinct_atoms("CCO", 6)
        1
    """

    from rdkit import Chem
    import numpy as np

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
        atom_indices = np.array(
            [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == atomic_number
            ]
        )
        # Count them
        return len(set(atom_ranks[atom_indices]))
    except (TypeError, ValueError, AttributeError):
        return "Error: Not a valid SMILES string"
