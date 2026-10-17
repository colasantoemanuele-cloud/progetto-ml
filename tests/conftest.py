"""Utilita' condivise dai test: caricamento degli script di `src/` per percorso.

I nomi degli script iniziano con una cifra (es. `05_baseline.py`) e non sono
importabili con `import`: vengono caricati con `importlib.util.spec_from_file_location`.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

_CACHE: dict[str, object] = {}


def load_script(name: str):
    """Carica `src/<name>.py` (o `src/sim/<name>.py`) come modulo, con cache."""
    if name in _CACHE:
        return _CACHE[name]
    path = SRC / f"{name}.py"
    if not path.exists():
        path = SRC / "sim" / f"{name}.py"
    mod_name = "s_" + name.replace("-", "_")
    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    _CACHE[name] = mod
    return mod
