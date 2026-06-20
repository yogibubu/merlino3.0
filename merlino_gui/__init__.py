"""Merlino4 GUI shell.

The legacy GUI remains under `gui/` during migration. New GUI code starts here
and talks to service layers instead of scientific backends directly.
"""

from .dashboard import DashboardWindow
from .workflow_registry import WorkflowSpec, default_workflows

__all__ = ["DashboardWindow", "WorkflowSpec", "default_workflows"]
