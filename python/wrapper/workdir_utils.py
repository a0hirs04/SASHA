from __future__ import annotations

import os
from pathlib import Path


def default_project_work_root(project_root: Path) -> Path:
    override = os.environ.get("PROJECT_NORTHSTAR_WORK_ROOT") or os.environ.get("NORTHSTAR_WORK_ROOT")
    if override:
        return Path(override).expanduser()

    user = os.environ.get("USER")
    if user:
        work_user = Path("/work") / user
        if work_user.exists() and os.access(work_user, os.W_OK):
            return work_user / "PROJECT-NORTHSTAR" / "build"

    return project_root / "build"


def default_reality_check_dir(project_root: Path, check_name: str) -> Path:
    return default_project_work_root(project_root) / check_name
