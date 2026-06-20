"""Merlino4 GUI shell.

The legacy GUI remains under `gui/` during migration. New GUI code starts here
and talks to service layers instead of scientific backends directly.
"""

from .dashboard import DashboardWindow
from .manifest_browser import ManifestBrowserWindow, discover_manifests
from .workflow_registry import WorkflowSpec, default_workflows

__all__ = ["DashboardWindow", "ManifestBrowserWindow", "WorkflowSpec", "default_workflows", "discover_manifests"]
