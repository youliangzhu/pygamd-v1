from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile


def run_gmx_dump(
    tpr_path: str | Path,
    *,
    gmx_bin: str = "gmx",
    keep_path: str | Path | None = None,
) -> Path:
    tpr = Path(tpr_path).resolve()
    cmd = [gmx_bin, "dump", "-s", str(tpr)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        stderr = "\n".join(result.stderr.splitlines()[:20])
        raise RuntimeError(
            "gmx dump failed\n"
            f"  gmx_bin: {gmx_bin}\n"
            f"  tpr_path: {tpr}\n"
            f"  returncode: {result.returncode}\n"
            f"  stderr:\n{stderr}"
        )

    if keep_path is not None:
        out = Path(keep_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(result.stdout, encoding="utf-8")
        return out

    with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".tpr.dmp", delete=False) as fh:
        fh.write(result.stdout)
        return Path(fh.name)
