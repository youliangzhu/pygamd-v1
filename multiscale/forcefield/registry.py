"""Force-field registry for Multiscale."""
from __future__ import annotations
from typing import Any
from .opls import summarize_opls

def list_forcefields() -> list[dict[str, Any]]:
    opls = summarize_opls()
    data = opls.assets.to_dict()
    data.update({"status": "available" if opls.ok else "missing-or-incomplete",
                 "atom_type_count": opls.atom_type_count,
                 "unique_moltemplate_types": opls.unique_moltemplate_types,
                 "elements": list(opls.elements)})
    return [data]
