import os
import pytest


@pytest.fixture(scope="function")
def app_name():
    chemenv_name = os.getenv("CHEMENV_NAME", "")
    return f"chemenv{chemenv_name}"
