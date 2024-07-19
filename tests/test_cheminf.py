import pytest
import modal

@pytest.mark.asyncio
async def test_tanimoto(app_name):
    fxn = modal.Function.lookup(app_name, "tanimoto")
    result = await fxn.remote.aio("CCO", "CC")
    assert abs(result - 0.143) < 0.001
