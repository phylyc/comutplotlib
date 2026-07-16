"""Configuration loading for comutplotlib.

Institution/study-specific vocabulary (gene-symbol aliases, default metadata
rows, sample-classification rules, palette extensions) is shipped as editable
JSON *templates* inside this package rather than being hard-coded in the source.

Resolution order for every config file:

1. ``$COMUTPLOTLIB_CONFIG_DIR/<filename>`` — a directory of user overrides,
   if the environment variable is set and the file exists. This lets
   collaborators point at their own adapted copies without touching the
   installed package.
2. The template shipped inside ``comutplotlib/config/``.
3. The ``default`` argument passed to :func:`load_config` (used only if the
   packaged template is also unavailable, e.g. a stripped-down install).
"""

import json
import os
from importlib import resources

ENV_CONFIG_DIR = "COMUTPLOTLIB_CONFIG_DIR"


def config_dir() -> str | None:
    """Return the user override directory, if configured."""
    return os.environ.get(ENV_CONFIG_DIR)


def _packaged_config_text(filename: str) -> str:
    return resources.files(__package__).joinpath(filename).read_text(encoding="utf-8")


def load_config(filename: str, default=None):
    """Load a JSON config file following the documented resolution order."""
    override_dir = config_dir()
    if override_dir:
        path = os.path.join(override_dir, filename)
        if os.path.exists(path):
            with open(path, encoding="utf-8") as f:
                return json.load(f)
    try:
        return json.loads(_packaged_config_text(filename))
    except (FileNotFoundError, ModuleNotFoundError, json.JSONDecodeError):
        return default

