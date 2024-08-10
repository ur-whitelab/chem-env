from chemenv.tools.cheminformatics import tanimoto
import pytest
def test_tanimoto():
    assert pytest.approx(tanimoto("CCO", "CC") , abs=0.001) == 0.143
    assert tanimoto("CCO", "CCO") == 1