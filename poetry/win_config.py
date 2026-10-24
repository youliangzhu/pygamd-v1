#!/usr/bin/env python3
"""Copy CUDA runtime DLLs next to this script on Windows.

This script is intentionally standalone. It does not read CMakeCache.txt or
any other project file, so it can be copied together with cu_gala.pyd to a new
directory and run there.
"""

import argparse
import os
import shutil
from pathlib import Path


CUDA_DLL_PATTERNS = (
    "cudart64*.dll",
    "cufft64*.dll",
)


def as_dir(value):
    if not value:
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_dir() else None


def unique_dirs(paths):
    result = []
    seen = set()
    for path in paths:
        if not path:
            continue
        path = path.resolve()
        key = os.path.normcase(str(path))
        if path.is_dir() and key not in seen:
            seen.add(key)
            result.append(path)
    return result


def cuda_roots_from_system():
    roots = []

    for name in ("CUDA_PATH", "CUDA_HOME", "CUDA_ROOT"):
        root = as_dir(os.environ.get(name))
        if root:
            roots.append(root)

    nvcc = shutil.which("nvcc")
    if nvcc:
        roots.append(Path(nvcc).resolve().parent.parent)

    for env_name in ("ProgramFiles", "ProgramFiles(x86)"):
        program_files = as_dir(os.environ.get(env_name))
        if not program_files:
            continue
        cuda_parent = program_files / "NVIDIA GPU Computing Toolkit" / "CUDA"
        if cuda_parent.is_dir():
            roots.extend(sorted(cuda_parent.glob("v*"), reverse=True))

    return unique_dirs(roots)


def runtime_dirs(cuda_root):
    return (
        cuda_root / "bin" / "x64",
        cuda_root / "bin",
    )


def dlls_from_root(cuda_root):
    found = {}
    searched = []

    for folder in runtime_dirs(cuda_root):
        searched.append(folder)
        if not folder.is_dir():
            continue
        for pattern in CUDA_DLL_PATTERNS:
            for dll in folder.glob(pattern):
                found.setdefault(dll.name.lower(), dll.resolve())

    return sorted(found.values(), key=lambda item: item.name.lower()), searched


def choose_cuda_runtime(cuda_roots):
    all_searched = []

    for root in unique_dirs(cuda_roots):
        dlls, searched = dlls_from_root(root)
        all_searched.extend(searched)
        names = [dll.name.lower() for dll in dlls]
        has_cudart = any(name.startswith("cudart64") for name in names)
        has_cufft = any(name.startswith("cufft64") for name in names)
        if has_cudart and has_cufft:
            return root, dlls, all_searched

    return None, [], all_searched


def copy_dlls(dlls, destination, force=False, dry_run=False):
    destination.mkdir(parents=True, exist_ok=True)
    copied = []
    skipped = []

    for dll in dlls:
        target = destination / dll.name
        if target.exists() and not force:
            skipped.append(target)
            continue
        copied.append(target)
        if not dry_run:
            shutil.copy2(dll, target)

    return copied, skipped


def parse_args():
    parser = argparse.ArgumentParser(description="Copy CUDA runtime DLLs next to cu_gala.pyd.")
    parser.add_argument("--cuda-root", action="append", default=[],
                        help="CUDA Toolkit root, for example C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v13.3")
    parser.add_argument("--dest", default=None,
                        help="Destination directory. Defaults to the directory containing this script.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing DLLs.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be copied without changing files.")
    return parser.parse_args()


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    destination = Path(args.dest).expanduser().resolve() if args.dest else script_dir

    cuda_roots = [Path(path).expanduser() for path in args.cuda_root]
    cuda_roots.extend(cuda_roots_from_system())

    cuda_root, dlls, searched = choose_cuda_runtime(cuda_roots)
    if not cuda_root:
        print("CUDA runtime DLLs were not found.")
        print("Searched directories:")
        for folder in searched:
            print(f"  {folder}")
        print("Try: python win_config.py --cuda-root \"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v13.3\"")
        return 1

    copied, skipped = copy_dlls(dlls, destination, force=args.force, dry_run=args.dry_run)
    action = "Would copy" if args.dry_run else "Copied"

    print(f"CUDA root: {cuda_root}")
    print(f"Destination: {destination}")
    print(f"Found {len(dlls)} CUDA runtime DLL(s).")
    for target in copied:
        print(f"{action}: {target.name}")
    if skipped:
        print(f"Skipped {len(skipped)} existing DLL(s). Use --force to overwrite.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
