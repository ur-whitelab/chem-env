from chemenv.tools.cheminformatics import (
    tanimoto,
    get_number_of_topologically_distinct_atoms,
)
import pytest


def test_tanimoto():
    assert pytest.approx(tanimoto("CCO", "CC"), abs=0.001) == 0.143
    assert tanimoto("CCO", "CCO") == 1


def test_get_number_of_topologically_distinct_atoms():
    assert get_number_of_topologically_distinct_atoms("CCO", 1) == 2
    assert get_number_of_topologically_distinct_atoms("CCO", 6) == 1
