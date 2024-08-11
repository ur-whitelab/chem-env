from modal import App
from modal import Image
from chemenv.tools.cheminformatics import (
    tanimoto,
    get_number_of_topologically_distinct_atoms,
)

rdkit_image = (
    Image.debian_slim(python_version="3.12").pip_install("rdkit").pip_install("numpy")
)

app = App()

tanimoto = app.function(tanimoto, image=rdkit_image)
get_number_of_topologically_distinct_atoms = app.function(
    get_number_of_topologically_distinct_atoms, image=rdkit_image
)
