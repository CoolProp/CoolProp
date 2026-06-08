#!/usr/bin/env python
"""
Test script to compare wheel contents between old setup.py and new scikit-build-core builds.

This script:
1. Builds a wheel using the old setup.py approach
2. Builds a wheel using the new scikit-build-core approach
3. Extracts and compares the contents
4. Reports any differences
5. Exits with code 0 if identical, 1 if different
"""

import sys
import os
import tempfile
import subprocess
import zipfile
import shutil
import filecmp
import difflib
from pathlib import Path
from typing import Tuple, List, Set

# Get the repository root
REPO_ROOT = Path(__file__).parent.parent.absolute()

def run_command(cmd: List[str], cwd: Path, description: str) -> subprocess.CompletedProcess:
    """Run a command and handle errors."""
    print(f"  Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {description} failed!")
        print(f"  Command: {' '.join(cmd)}")
        print(f"  Return code: {e.returncode}")
        print(f"  stdout: {e.stdout}")
        print(f"  stderr: {e.stderr}")
        sys.exit(1)

def build_old_wheel(tmpdir: Path) -> Path:
    """Build wheel using old setup.py approach."""
    print("\n" + "="*80)
    print("Building wheel with OLD setup.py approach...")
    print("="*80)

    # Change to wrappers/Python directory
    build_dir = REPO_ROOT / "wrappers" / "Python"

    # Clean any previous builds
    dist_dir = build_dir / "dist"
    if dist_dir.exists():
        shutil.rmtree(dist_dir)

    # Build the wheel
    run_command(
        [sys.executable, "setup.py", "bdist_wheel"],
        cwd=build_dir,
        description="Old wheel build"
    )

    # Find the built wheel
    wheels = list((build_dir / "dist").glob("*.whl"))
    if not wheels:
        print("ERROR: No wheel file found in dist/")
        sys.exit(1)

    wheel_path = wheels[0]
    print(f"  Built wheel: {wheel_path.name}")

    # Copy to temp directory
    dest = tmpdir / wheel_path.name
    shutil.copy2(wheel_path, dest)

    return dest

def build_new_wheel(tmpdir: Path) -> Path:
    """Build wheel using new scikit-build-core approach."""
    print("\n" + "="*80)
    print("Building wheel with NEW scikit-build-core approach...")
    print("="*80)

    # Build from repository root
    run_command(
        [sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
         "-w", str(tmpdir), str(REPO_ROOT)],
        cwd=REPO_ROOT,
        description="New wheel build"
    )

    # Find the built wheel
    wheels = list(tmpdir.glob("*.whl"))
    if not wheels:
        print("ERROR: No wheel file found")
        sys.exit(1)

    wheel_path = wheels[0]
    print(f"  Built wheel: {wheel_path.name}")

    return wheel_path

def extract_wheel(wheel_path: Path, extract_dir: Path) -> Set[str]:
    """Extract wheel and return set of relative file paths."""
    print(f"\n  Extracting {wheel_path.name}...")

    with zipfile.ZipFile(wheel_path, 'r') as zf:
        zf.extractall(extract_dir)

    # Get all files recursively
    files = set()
    for root, dirs, filenames in os.walk(extract_dir):
        for filename in filenames:
            filepath = Path(root) / filename
            relpath = filepath.relative_to(extract_dir)
            files.add(str(relpath))

    print(f"    Extracted {len(files)} files")
    return files

def compare_file_contents(file1: Path, file2: Path) -> Tuple[bool, str]:
    """Compare two files. Returns (identical, diff_text)."""
    # Binary comparison first
    if filecmp.cmp(file1, file2, shallow=False):
        return True, ""

    # If different, try to generate a useful diff
    try:
        with open(file1, 'r', encoding='utf-8', errors='ignore') as f1:
            lines1 = f1.readlines()
        with open(file2, 'r', encoding='utf-8', errors='ignore') as f2:
            lines2 = f2.readlines()

        diff = list(difflib.unified_diff(
            lines1, lines2,
            fromfile=str(file1),
            tofile=str(file2),
            lineterm='',
            n=3
        ))

        if diff:
            return False, '\n'.join(diff[:50])  # Limit diff output
    except:
        pass

    return False, "(binary files differ)"

def compare_wheels(old_dir: Path, new_dir: Path, old_files: Set[str], new_files: Set[str]) -> bool:
    """Compare extracted wheel contents. Returns True if identical."""
    print("\n" + "="*80)
    print("Comparing wheel contents...")
    print("="*80)

    # Find differences in file lists
    only_in_old = old_files - new_files
    only_in_new = new_files - old_files
    common_files = old_files & new_files

    # Filter out metadata files and irrelevant files that are expected to differ
    def is_ignorable(path: str) -> bool:
        parts = Path(path).parts
        filename = Path(path).name
        return (len(parts) > 0 and
                (parts[0].endswith('.dist-info') or
                 parts[0].endswith('.egg-info') or
                 path.endswith('RECORD') or
                 path.endswith('WHEEL') or
                 filename == '.DS_Store' or
                 filename == '.gitignore' or
                 path.startswith('_py_backend/')))

    only_in_old = {f for f in only_in_old if not is_ignorable(f)}
    only_in_new = {f for f in only_in_new if not is_ignorable(f)}

    all_identical = True

    # Report files only in old
    if only_in_old:
        print(f"\n⚠️  Files ONLY in OLD wheel ({len(only_in_old)}):")
        for f in sorted(only_in_old)[:20]:  # Limit output
            print(f"    - {f}")
        if len(only_in_old) > 20:
            print(f"    ... and {len(only_in_old) - 20} more")
        all_identical = False

    # Report files only in new
    if only_in_new:
        print(f"\n⚠️  Files ONLY in NEW wheel ({len(only_in_new)}):")
        for f in sorted(only_in_new)[:20]:  # Limit output
            print(f"    - {f}")
        if len(only_in_new) > 20:
            print(f"    ... and {len(only_in_new) - 20} more")
        all_identical = False

    # Compare common files (excluding binary .so files which are expected to differ)
    non_binary_files = [f for f in common_files
                        if not is_ignorable(f) and not f.endswith('.so')]

    print(f"\n  Comparing {len(non_binary_files)} common non-binary files...")
    print(f"  (Skipping {len([f for f in common_files if f.endswith('.so')])} binary .so files)")
    different_files = []

    for relpath in sorted(non_binary_files):
        file1 = old_dir / relpath
        file2 = new_dir / relpath

        identical, diff_text = compare_file_contents(file1, file2)
        if not identical:
            different_files.append((relpath, diff_text))

    if different_files:
        print(f"\n⚠️  Files with DIFFERENT contents ({len(different_files)}):")
        for relpath, diff_text in different_files[:10]:  # Limit output
            print(f"\n    File: {relpath}")
            if diff_text:
                print("    Diff preview:")
                for line in diff_text.split('\n')[:20]:
                    print(f"      {line}")
        if len(different_files) > 10:
            print(f"    ... and {len(different_files) - 10} more files differ")
        all_identical = False

    # Summary
    print("\n" + "="*80)
    if all_identical:
        print("✅ SUCCESS: Wheel contents are IDENTICAL!")
    else:
        print("❌ FAILURE: Wheel contents DIFFER!")
        print(f"  - Files only in old: {len(only_in_old)}")
        print(f"  - Files only in new: {len(only_in_new)}")
        print(f"  - Files with different content: {len(different_files)}")
    print("="*80)

    return all_identical

def main():
    print("="*80)
    print("CoolProp Wheel Comparison Test")
    print("="*80)
    print(f"Repository: {REPO_ROOT}")

    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir = Path(tmpdir_str)

        # Create subdirectories
        old_build = tmpdir / "old_build"
        new_build = tmpdir / "new_build"
        old_extract = tmpdir / "old_extract"
        new_extract = tmpdir / "new_extract"

        for d in [old_build, new_build, old_extract, new_extract]:
            d.mkdir()

        # Build both wheels
        old_wheel = build_old_wheel(old_build)
        new_wheel = build_new_wheel(new_build)

        # Extract both wheels
        old_files = extract_wheel(old_wheel, old_extract)
        new_files = extract_wheel(new_wheel, new_extract)

        # Compare contents
        identical = compare_wheels(old_extract, new_extract, old_files, new_files)

        # Exit with appropriate code
        sys.exit(0 if identical else 1)

if __name__ == "__main__":
    main()
