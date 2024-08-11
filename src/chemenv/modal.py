from modal import App
from modal import Image
from chemenv.tools.cheminformatics import (
    get_tanimoto_similarity,
    get_number_of_topologically_distinct_atoms,
    get_element_info,
)

rdkit_image = (
    Image.debian_slim(python_version="3.12").pip_install("rdkit").pip_install("numpy")
)

mendeleev_image = Image.debian_slim().pip_install("mendeleev")

app = App()

get_tanimoto_similarity = app.function(get_tanimoto_similarity, image=rdkit_image)
get_number_of_topologically_distinct_atoms = app.function(
    get_number_of_topologically_distinct_atoms, image=rdkit_image
)
get_element_info = app.function(get_element_info, image=mendeleev_image)
