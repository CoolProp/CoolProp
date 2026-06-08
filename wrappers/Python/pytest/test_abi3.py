"""Guard the abi3 single-wheel invariant (bd CoolProp-r9sq.19).

The v8 nanobind package, when built with ``COOLPROP_NANOBIND=ON`` on Python
>= 3.12, ships as ONE ``cp312-abi3`` wheel that must run unchanged on
3.12/3.13/3.14.  The load-bearing invariant this module enforces:

    if the installed wheel is tagged abi3, every shipped extension module must
    be limited-API (``.abi3.so``).

A version-specific ``.so`` inside an abi3-tagged wheel imports fine on the build
interpreter but fails to load on every *other* minor version -- exactly the
silent cross-version breakage abi3 is meant to prevent, and the kind of mixed
artifact that is otherwise invisible until a user on a different Python hits it.

This is a *consistency* check, valid for any install.  A per-version nanobind
wheel -- a define-only build (``-DCOOLPROP_NANOBIND=ON`` without the env signal
that triggers ``wheel.py-api``), or any build on Python < 3.12 below nanobind's
stable-ABI floor -- is correctly tagged version-specific and legitimately
carries version-specific sidecars, so the invariant does not apply and the test
skips.  That the *shipped* wheel actually IS abi3 is enforced separately by the
``nanobind-abi3`` CI job, which asserts a ``cp312-abi3`` wheel is produced and
then installs + runs it on 3.12, 3.13 AND 3.14.
"""
import sys
import importlib
import importlib.metadata as importlib_metadata

import pytest

# nanobind-only: the capsule State shim is absent in the legacy build, so this
# whole module is skipped there rather than reported as a failure.
pytest.importorskip("CoolProp.State", reason="legacy (non-nanobind) build")

# The nanobind core, the capsule State shim, and the _constants cimport-surface
# backing module -- all three ship inside the wheel.
_MODULES = ["CoolProp.CoolProp", "CoolProp.State", "CoolProp._constants"]


def _installed_wheel_is_abi3():
    """True iff the installed CoolProp distribution's wheel tag is abi3.

    Reads the wheel's recorded ``Tag:`` lines from its dist-info ``WHEEL`` file.
    Returns False when there is no wheel metadata (editable / in-tree install)
    or the dist cannot be found -- in those cases there is no abi3 claim to hold
    the extensions to, so the invariant simply does not apply.
    """
    try:
        wheel_meta = importlib_metadata.distribution("CoolProp").read_text("WHEEL") or ""
    except importlib_metadata.PackageNotFoundError:  # pragma: no cover
        return False
    tags = [
        line.split(":", 1)[1].strip()
        for line in wheel_meta.splitlines()
        if line.startswith("Tag:")
    ]
    return any("abi3" in t for t in tags)


@pytest.mark.skipif(
    not sys.platform.startswith(("linux", "darwin")),
    reason="abi3 SOABI filename check is POSIX-specific (.abi3.so)",
)
@pytest.mark.skipif(
    not _installed_wheel_is_abi3(),
    reason="per-version wheel (or no wheel metadata): version-specific .so are correct here",
)
@pytest.mark.parametrize("modname", _MODULES)
def test_abi3_wheel_ships_only_abi3_extensions(modname):
    """An abi3-tagged wheel must contain only limited-API ``.abi3.so``."""
    mod = importlib.import_module(modname)
    path = getattr(mod, "__file__", "") or ""
    assert path.endswith(".abi3.so"), (
        "%s loaded from %r inside an abi3-tagged wheel; expected a limited-API "
        ".abi3.so.  A version-specific .so in an abi3 wheel imports on the build "
        "interpreter but breaks on other minor versions (r9sq.19)." % (modname, path)
    )
