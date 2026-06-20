"""Core infrastructure for Merlino4.

This package is intentionally small at the start of the refactor. It provides
stable helpers that older Merlino3 modules can adopt incrementally.
"""

from .config import MerlinoConfig, load_config, write_default_config
from .errors import BackendError, InputError, MerlinoError, ParseError, ScientificValidationError
from .manifest import RunManifest, build_run_manifest, file_checksums, sha256_file, write_manifest
from .paths import repo_root
from .project import ProjectState, ensure_project_state
from .workspace import WorkspaceLayout, ensure_workspace, slugify

__all__ = [
    "BackendError",
    "InputError",
    "MerlinoConfig",
    "MerlinoError",
    "ParseError",
    "ProjectState",
    "RunManifest",
    "ScientificValidationError",
    "WorkspaceLayout",
    "build_run_manifest",
    "ensure_project_state",
    "ensure_workspace",
    "file_checksums",
    "load_config",
    "repo_root",
    "sha256_file",
    "slugify",
    "write_default_config",
    "write_manifest",
]
