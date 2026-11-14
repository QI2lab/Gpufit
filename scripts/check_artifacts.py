"""Validate wheel/sdist contents and install artifacts outside the checkout."""

import argparse
import json
import shutil
import subprocess
import tempfile
import tomllib
import zipfile
from pathlib import Path
from tarfile import open as open_tar


def checked(*args: str, cwd: Path) -> None:
    """Run an acceptance command and propagate its failure."""
    subprocess.run(args, cwd=cwd, check=True)


def check_wheel(wheel: Path, test_source: Path) -> dict:
    """Install a locked artifact profile in a fresh temporary directory and test it."""
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        assert "pygpufit/build_manifest.json" in names
        assert "pygpufit/licenses/NVIDIA-CUDA-EULA.txt" in names
        assert any(name.endswith("licenses/LICENSE.txt") for name in names)
        assert not any(".venv/" in name or ".uv-cache/" in name for name in names)
        identity = json.loads(archive.read("pygpufit/build_manifest.json"))
        print(json.dumps(identity, indent=2), flush=True)
    uv = shutil.which("uv")
    if uv is None:
        raise RuntimeError("uv must be installed before artifact acceptance")
    with tempfile.TemporaryDirectory(prefix="gpufit-wheel-acceptance-") as directory:
        target = Path(directory)
        (target / "pyproject.toml").write_text(
            '[project]\nname="gpufit-artifact-acceptance"\nversion="0.0.0"\n'
            'requires-python=">=3.12"\ndependencies=["pyGpufit=='
            + identity["version"]
            + '","pytest==9.1.1"]\n'
            '[tool.uv.sources]\npygpufit={path="' + wheel.as_posix() + '"}\n',
            encoding="utf-8",
        )
        shutil.copy2(test_source, target / "test_cuda_constraints.py")
        checked(uv, "lock", cwd=target)
        checked(uv, "sync", "--locked", cwd=target)
        checked(uv, "sync", "--locked", cwd=target)
        checked(uv, "run", "--locked", "--no-sync", "pytest", "-q", cwd=target)
        checked(
            uv,
            "run",
            "--locked",
            "--no-sync",
            "python",
            "-c",
            "import json,pygpufit.gpufit as gf; "
            "a=gf.get_capabilities(); print(json.dumps(a)); "
            "assert a['source_revision']==" + repr(identity["source_revision"]) + "; "
            "assert a['source_dirty']==" + repr(identity["source_dirty"]),
            cwd=target,
        )
    return identity


def main() -> None:
    """Check the built wheel, then optionally rebuild and install the source archive."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-sdist", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    version = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"][
        "version"
    ]
    wheel = next((root / "dist").glob(f"*-{version}-*.whl"))
    source_archive = next((root / "dist").glob(f"*-{version}.tar.gz"))
    original = check_wheel(wheel, root / "Gpufit/python/tests/test_cuda_constraints.py")
    with open_tar(source_archive) as archive:
        names = archive.getnames()
        assert any(name.endswith("/cmake/source_identity.json") for name in names)
        assert not any(".uv-cache/" in name or ".venv/" in name for name in names)
        assert any(name.endswith("/Gpufit/models/gauss_3d_arb.cuh") for name in names)
        if not args.rebuild_sdist:
            return
        with tempfile.TemporaryDirectory(prefix="gpufit-sdist-acceptance-") as directory:
            target = Path(directory)
            archive.extractall(target, filter="data")
            extracted = next(path for path in target.iterdir() if path.is_dir())
            checked(
                "powershell",
                "-NoProfile",
                "-ExecutionPolicy",
                "Bypass",
                "-File",
                str(extracted / "scripts/build_windows.ps1"),
                "-Architectures",
                "86",
                cwd=extracted,
            )
            rebuilt = check_wheel(
                next((extracted / "dist").glob("*.whl")),
                extracted / "Gpufit/python/tests/test_cuda_constraints.py",
            )
            assert rebuilt["source_revision"] == original["source_revision"]
            assert rebuilt["source_dirty"] == original["source_dirty"]


if __name__ == "__main__":
    main()
