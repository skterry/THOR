"""
In-tree PEP 517 build backend for THOR.

This thin wrapper around setuptools' standard build backend guarantees that the
custom setup steps (compile the Fortran sources, download the data archives) run
during ``pip install -e .`` (which calls ``build_editable``) and ``pip install
.`` (which calls ``build_wheel``).

We use a build backend rather than only a ``cmdclass`` override because modern
setuptools does not invoke ``build_py`` for an editable install of a project
that ships no importable Python packages — so a ``build_py`` hook alone would
never fire. ``build_editable`` always runs, so hooking it is reliable.

Wired up in pyproject.toml via::

    [build-system]
    build-backend = "thor_build"
    backend-path = ["."]
"""

# Re-export every standard PEP 517 hook (get_requires_*, prepare_metadata_*,
# build_sdist, build_wheel, build_editable, ...) so this module is a complete
# backend, then override the two build hooks below.
from setuptools import build_meta as _orig
from setuptools.build_meta import *  # noqa: F401,F403

# run_thor_setup() lives in setup.py; importing it does not start a build
# because the setup() call there is guarded by `if __name__ == "__main__"`.
from setup import run_thor_setup


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    run_thor_setup()
    return _orig.build_wheel(wheel_directory, config_settings, metadata_directory)


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    run_thor_setup()
    return _orig.build_editable(wheel_directory, config_settings, metadata_directory)
