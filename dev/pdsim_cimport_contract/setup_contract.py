"""
Build the PDSim cimport-contract extension against the *installed* CoolProp.

This mirrors exactly how ibell/pdsim builds: include dirs come from the
installed CoolProp package (its shipped headers + the .pxd cimport surface),
and crucially there is NO link step against libCoolProp -- the cimported cdef
classes' methods resolve at runtime through Cython's vtable capsule.  If a
future v8 build breaks that, this build (and the test) fails.

    python setup_contract.py build_ext --inplace
"""
import os

import CoolProp
from setuptools import Extension, setup
from Cython.Build import cythonize

pkg = os.path.dirname(CoolProp.__file__)
site_packages = os.path.dirname(pkg)   # cimport root: `CoolProp.State` etc. resolve here
include_dirs = [
    os.path.join(pkg, "include"),   # CoolProp C++ headers (AbstractState.h, ...)
    pkg,
    CoolProp.get_include_directory(),
]


def _find_dep_include(header_relpath, env_var, fetch_glob):
    """STOPGAP for header leakage, not a real PDSim dependency.

    CoolProp's public/cimport surface currently *leaks* fmt as a transitive
    include (fmt via CPstrings.h's inline format() helpers).  It does NOT
    appear in any State/AbstractState API signature.  Until the header is
    de-leaked (see SURFACE.md finding 3), this build needs the extra -I path;
    after that, delete these calls and the contract still compiles with only
    CoolProp.get_include_directory()."""
    import glob
    import subprocess
    # Search build trees in every checkout that shares this repo: this worktree
    # AND the main checkout (a fresh worktree has no build*/_deps of its own).
    roots = [os.path.abspath(os.path.join(os.path.dirname(__file__), *([".."] * 2)))]
    for cmd in (["git", "rev-parse", "--show-toplevel"],
                ["git", "rev-parse", "--git-common-dir"]):
        try:
            out = subprocess.check_output(cmd, cwd=os.path.dirname(__file__),
                                          text=True, stderr=subprocess.DEVNULL).strip()
            if out:
                # --git-common-dir gives <main>/.git; its parent is the checkout
                roots.append(os.path.dirname(out) if out.endswith(".git") else out)
        except Exception:
            # git missing / not a repo / command failed: best-effort only, fall
            # back to the other root candidates (env var, pkg include, system).
            pass
    candidates = []
    if os.environ.get(env_var):
        candidates.append(os.environ[env_var])
    candidates.append(os.path.join(pkg, "include"))          # the v8 goal: bundled
    for root in dict.fromkeys(os.path.abspath(r) for r in roots):
        candidates += glob.glob(os.path.join(root, fetch_glob))
    candidates += ["/opt/homebrew/include", "/usr/local/include"]
    for c in candidates:
        if c and os.path.isfile(os.path.join(c, header_relpath)):
            return c
    raise SystemExit(
        f"PDSim contract: could not find {header_relpath} (a CoolProp header "
        f"dependency).  Set {env_var} to a dir containing it.")


include_dirs.append(_find_dep_include(
    "fmt/format.h", "COOLPROP_FMT_INCLUDE", "build*/_deps/fmt-src/include"))

ext = Extension(
    "pdsim_surface",
    ["pdsim_surface.pyx"],
    language="c++",
    include_dirs=include_dirs,
    extra_compile_args=["-std=c++17"],   # CoolProp public headers use std::filesystem
    # No libraries / library_dirs: the contract is link-free by design.
)

setup(
    name="pdsim_surface",
    # cython's include_path must point at the cimport root (site-packages) so
    # `from CoolProp.State cimport State` and the package's own relative
    # cimports (`from . cimport constants_header`) both resolve.
    ext_modules=cythonize([ext], language_level=3, include_path=[site_packages]),
)
