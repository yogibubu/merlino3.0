from __future__ import annotations

import os

import pytest

from merlino_gaussian import (
    GaussianInputError,
    ensure_gjf_input,
    gaussian_completion_message,
    select_latest_log,
)


def test_ensure_gjf_input_prefers_existing_gjf(tmp_path):
    (tmp_path / "gauin").write_text("raw\n", encoding="utf-8")
    gjf = tmp_path / "gauin.gjf"
    gjf.write_text("gjf\n", encoding="utf-8")

    assert ensure_gjf_input(tmp_path) == gjf
    assert gjf.read_text(encoding="utf-8") == "gjf\n"


def test_ensure_gjf_input_copies_raw_gauin(tmp_path):
    (tmp_path / "gauin").write_text("raw\n", encoding="utf-8")

    assert ensure_gjf_input(tmp_path) == tmp_path / "gauin.gjf"
    assert (tmp_path / "gauin.gjf").read_text(encoding="utf-8") == "raw\n"


def test_ensure_gjf_input_requires_input(tmp_path):
    with pytest.raises(GaussianInputError):
        ensure_gjf_input(tmp_path)


def test_select_latest_log_uses_mtime_then_size(tmp_path):
    old_log = tmp_path / "gauin.log"
    new_log = tmp_path / "gauout.log"
    old_log.write_text("old\n", encoding="utf-8")
    new_log.write_text("newer\n", encoding="utf-8")
    os.utime(old_log, (1, 1))
    os.utime(new_log, (2, 2))

    assert select_latest_log(tmp_path) == new_log


def test_gaussian_completion_message_detects_normal_termination(tmp_path):
    (tmp_path / "gauin.log").write_text(
        "Normal termination of Gaussian 16\n",
        encoding="utf-8",
    )

    success, message = gaussian_completion_message(tmp_path, exit_code=1)

    assert success
    assert message == "Gaussian completed successfully"
