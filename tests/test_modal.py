import pytest
import modal


@pytest.mark.asyncio
async def test_get_tanimoto(app_name):
    fxn = modal.Function.lookup(app_name, "get_tanimoto")
    result = await fxn.remote.aio("CCO", "CC")
    assert abs(result - 0.143) < 0.001


@pytest.mark.asyncio
async def test_get_number_of_topologically_distinct_atoms(app_name):
    fxn = modal.Function.lookup(app_name, "get_number_of_topologically_distinct_atoms")
    result = await fxn.remote.aio("CCO", 6)
    assert result == 2


@pytest.mark.asyncio
async def test_get_element_info(app_name):
    fxn = modal.Function.lookup(app_name, "get_element_info")
    result = await fxn.remote.aio("H")
    assert result["name"] == "Hydrogen"
    assert result["symbol"] == "H"
    assert result["atomic_number"] == 1
    assert result["mass"] == 1.008
