"""Fortran backend discovery and execution contracts for Merlino4."""

from .backends import BackendSpec, backend_executable, resolve_backend

__all__ = ["BackendSpec", "backend_executable", "resolve_backend"]
