from modal import App
import os

from cheminf.cheminf import app as cheminf_app

# set a post-fix for testing/staging
chemenv_name = os.getenv("CHEMENV_NAME", "")
if chemenv_name and not chemenv_name.startswith("-"):
    chemenv_name = f"-{chemenv_name}"
app = App(f"chemenv{chemenv_name}")
app.include(cheminf_app)


