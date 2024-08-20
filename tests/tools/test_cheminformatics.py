from chemenv.tools.cheminformatics import (
    get_tanimoto_similarity,
    get_number_of_topologically_distinct_atoms,
    get_element_info,
)
import pytest


def test_get_tanimoto_similarity():
    assert pytest.approx(get_tanimoto_similarity("CCO", "CC"), abs=0.001) == 0.143
    assert get_tanimoto_similarity("CCO", "CCO") == 1


def test_get_number_of_topologically_distinct_atoms():
    assert get_number_of_topologically_distinct_atoms("CCO", 1) == 3
    assert get_number_of_topologically_distinct_atoms("CCO", 6) == 2


def test_get_element_info():
    assert get_element_info("H")["name"] == "Hydrogen"
    assert get_element_info("H")["symbol"] == "H"
    assert get_element_info("H")["atomic_number"] == 1
    assert get_element_info("H")["mass"] == 1.008

    assert get_element_info(1)["name"] == "Hydrogen"
