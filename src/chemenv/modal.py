from modal import App
from modal import Image
from chemenv.tools.cheminformatics import tanimoto

rdkit_image = Image.debian_slim(python_version="3.12").pip_install("rdkit")

app = App()

tanimoto = app.function(tanimoto,image=rdkit_image)