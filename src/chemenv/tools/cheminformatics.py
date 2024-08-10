def tanimoto(s1: str, s2: str) -> float:
    """Calculate the Tanimoto similarity of two SMILES strings."""
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