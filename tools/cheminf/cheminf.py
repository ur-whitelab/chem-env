from modal import App
from modal import Image

rdkit_image = Image.debian_slim(python_version="3.12").pip_install("rdkit")
with rdkit_image.imports():
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem


app = App()


@app.function(image=rdkit_image)
def tanimoto(s1: str, s2: str) -> float:
    """Calculate the Tanimoto similarity of two SMILES strings."""
    try:
        mol1 = Chem.MolFromSmiles(s1)
        mol2 = Chem.MolFromSmiles(s2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except (TypeError, ValueError, AttributeError) as e:
        raise ValueError("Invalid SMILES strings") from e
