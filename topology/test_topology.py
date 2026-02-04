"""Compatibility wrapper: canonical implementation is in merlino_fit.topology.test_topology."""

from merlino_fit.topology import test_topology as _impl
from merlino_fit.topology.test_topology import *  # noqa: F401,F403

# Re-export private helper used by existing tests/importers.
_parse_gaussian_topology_overrides = _impl._parse_gaussian_topology_overrides
