import pytest
from your_module import get_element_info

def test_get_element_info():
    # Test with a known element (Oxygen)
    result = get_element_info("Oxygen")

    expected = {
        "name": "Oxygen",
        "symbol": "O",
        "atomic_number": 8,
        "mass": 15.999,
        "electron_configuration": "1s2 2s2 2p4",
        "electronegativity": 3.44,
        "group": "Chalcogen",
        "period": 2,
        "block": "p",
    }

    assert result == expected

def test_invalid_element():
    result = get_element_info("Invalid")
    assert result == "Error: 'Invalid' is not a valid element identifier."

# Optional: parametrized test for multiple elements
@pytest.mark.parametrize("input_element, expected_symbol", [
    ("Hydrogen", "H"),
    ("Carbon", "C"),
    ("Iron", "Fe"),
    ("Gold", "Au"),
])
def test_multiple_elements(input_element, expected_symbol):
    result = get_element_info(input_element)
    assert result["symbol"] == expected_symbol
