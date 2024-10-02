import sys
from pathlib import Path

from pydustmasker._version import VERSION

if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib


def test_versions_match():
    cargo = Path().absolute() / "Cargo.toml"
    with open(cargo, "rb") as f:
        data = tomllib.load(f)
        cargo_version = data["package"]["version"]

    assert VERSION == cargo_version
