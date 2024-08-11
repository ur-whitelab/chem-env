import os
import subprocess
import fire
from chemenv.modal_app import app


class ChemEnvCLI:
    def __init__(self):
        self.app = app

    def _get_modal_file(self):
        return os.path.join(os.path.dirname(__file__), "modal_app.py")

    def deploy(self):
        """Deploy the ChemEnv Modal app."""
        print("Deploying ChemEnv app")
        modal_file = self._get_modal_file()
        subprocess.run(["modal", "deploy", modal_file], check=True)

    def serve(self):
        """Serve the ChemEnv Modal app (for development)."""
        print("Serving ChemEnv app")
        modal_file = self._get_modal_file()
        subprocess.run(["modal", "serve", modal_file], check=True)

    # Add methods for specific functions if needed, for example:
    def tanimoto_similarity(self, smiles1, smiles2):
        """Calculate Tanimoto similarity between two SMILES strings."""
        result = self.app.get_tanimoto_similarity.remote(smiles1, smiles2)
        print(f"Tanimoto similarity: {result}")


def main():
    fire.Fire(ChemEnvCLI)


if __name__ == "__main__":
    main()
