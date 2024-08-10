import modal

stub = modal.Stub("element-info")

@stub.function(
    image=modal.Image.debian_slim().pip_install("mendeleev"),
)
def get_element_info(identifier: str):
    from mendeleev import element

    try:
        # Try to get the element
        if identifier.isdigit():
            el = element(int(identifier))
        else:
            el = element(identifier)

        # Collect basic information
        info = {
            "name": el.name,
            "symbol": el.symbol,
            "atomic_number": el.atomic_number,
            "mass": el.mass,
            "electron_configuration": el.electronic_configuration,
            "electronegativity": el.electronegativity,
            "group": el.group.name if el.group else None,
            "period": el.period,
            "block": el.block,
        }

        return info

    except ValueError:
        return {"error": f"'{identifier}' is not a valid element identifier."}
