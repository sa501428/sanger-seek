import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sanger_seek.devtools.demogen import generate_demo  # noqa: E402


@pytest.fixture(scope="session")
def demo_dir(tmp_path_factory) -> Path:
    return generate_demo(tmp_path_factory.mktemp("demo"))
