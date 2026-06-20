"""Core infrastructure for Merlino4.

This package is intentionally small at the start of the refactor. It provides
stable helpers that older Merlino3 modules can adopt incrementally.
"""

from .manifest import sha256_file, write_manifest
from .paths import repo_root
from .project import ProjectState, ensure_project_state

__all__ = [
    "ProjectState",
    "ensure_project_state",
    "repo_root",
    "sha256_file",
    "write_manifest",
]
