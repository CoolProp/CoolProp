"""
Tests: REFPROP backend resolves every CoolProp fluid alias correctly.

Every test checks only that AbstractState construction succeeds (or raises),
plus CAS-number identity where fluid identity needs to be confirmed.
No fluid-property calls are made.

Test groups
-----------
1. Parametrized alias coverage (TestREFPROPAliasResolution)
   For each CoolProp fluid whose REFPROP_NAME matches a .FLD/.PPF file on
   disk, every alias (including the canonical CoolProp name) must construct
   an AbstractState and report the same CAS number as the canonical REFPROP
   name.

2. REFPROP-only fluids (TestREFPROPOnlyFluids)
   Fluids present in REFPROP's FLUIDS/ directory that have no entry in
   CoolProp's fluid library.  These pass through unchanged and must
   construct successfully.

3. .FLD-suffixed names (TestDirectFLDPath)
   Names ending with a REFPROP file extension (.FLD, .fld, .PPF, .ppf)
   — with either the canonical REFPROP stem or a CoolProp alias — must
   construct successfully.  Absolute paths are not supported and must
   raise ValueError.

4. Mixture access patterns (TestMixtures)
   Combinations of alias / canonical REFPROP name / .FLD-suffixed name
   across multi-component mixtures must construct successfully.

5. Pseudopure (.PPF) fluids (TestPPFFluids)
   AIR, R410A, R407C, R507A, R404A — stored as .PPF in REFPROP.
   Aliases defined in CoolProp must construct successfully and report the
   correct CAS number.

6. Edge cases and negative tests (TestEdgeCases)
   Unknown fluids must raise ValueError; case-insensitive aliases must
   resolve to the correct CAS number.

7. Discovery smoke test (TestDiscovery)

All tests are skipped if REFPROP is not installed or its FLUIDS directory
cannot be found.
"""

from __future__ import annotations

import json
import os

import pytest

import CoolProp
import CoolProp.CoolProp as CP


# ---------------------------------------------------------------------------
# Shared infrastructure
# ---------------------------------------------------------------------------

def _find_refprop_fluids_dir() -> str | None:
    """Return the REFPROP FLUIDS directory, or None if not found."""
    candidates = [
        os.path.join(os.getenv("RPPREFIX", ""), "FLUIDS"),
        os.path.join(os.getenv("COOLPROP_REFPROP_ROOT", ""), "FLUIDS"),
        "/opt/refprop/FLUIDS",
        "/home/francesco/refprop/REFPROP/FLUIDS",
    ]
    for path in candidates:
        if path and os.path.isdir(path):
            return path
    return None


_FLUIDS_DIR: str | None = _find_refprop_fluids_dir()
_REFPROP_VERSION: str = CP.get_global_param_string("REFPROP_version")
_SKIP = not _REFPROP_VERSION or _FLUIDS_DIR is None
_SKIP_REASON = (
    "REFPROP not available"
    if not _REFPROP_VERSION
    else "REFPROP FLUIDS directory not found"
)


def _fld_stems() -> set[str]:
    """Upper-case stems of all .FLD files in the REFPROP FLUIDS directory."""
    if _FLUIDS_DIR is None:
        return set()
    return {
        os.path.splitext(f)[0].upper()
        for f in os.listdir(_FLUIDS_DIR)
        if f.upper().endswith(".FLD")
    }


def _ppf_stems() -> set[str]:
    """Upper-case stems of all .PPF files in the REFPROP FLUIDS directory."""
    if _FLUIDS_DIR is None:
        return set()
    return {
        os.path.splitext(f)[0].upper()
        for f in os.listdir(_FLUIDS_DIR)
        if f.upper().endswith(".PPF")
    }


def _cp_rp_index() -> dict[str, tuple[str, str, list[str]]]:
    """
    Return {cp_canonical: (rp_name, cp_canonical, aliases)} for every CoolProp
    fluid whose REFPROP_NAME appears in the FLUIDS directory.
    """
    available = _fld_stems() | _ppf_stems()
    index: dict[str, tuple[str, str, list[str]]] = {}
    for cp_name in CP.get_global_param_string("FluidsList").split(","):
        try:
            info = json.loads(CP.get_fluid_param_string(cp_name, "JSON"))[0]["INFO"]
        except Exception:
            continue
        rp_name: str = info["REFPROP_NAME"]
        if rp_name.upper() in available:
            index[cp_name] = (rp_name, info["NAME"], info["ALIASES"])
    return index


def _build_alias_cases() -> list[tuple[str, str, str]]:
    """
    Return (cp_canonical, rp_name, alias) for every fluid×alias combination.
    Each alias is tested independently.
    """
    cases: list[tuple[str, str, str]] = []
    for cp_name, (rp_name, canonical, aliases) in _cp_rp_index().items():
        for alias in [canonical] + aliases:
            cases.append((cp_name, rp_name, alias))
    return cases


_ALIAS_CASES = _build_alias_cases()
_CP_RP_INDEX = _cp_rp_index()

skip = pytest.mark.skipif(_SKIP, reason=_SKIP_REASON)


# ---------------------------------------------------------------------------
# 1. Parametrized alias coverage
# ---------------------------------------------------------------------------

@skip
@pytest.mark.parametrize(
    "cp_name, rp_name, alias",
    _ALIAS_CASES,
    ids=[f"{c[0]}::{c[2]}" for c in _ALIAS_CASES],
)
class TestREFPROPAliasResolution:
    """One instance per (fluid, alias) — covers every alias for every fluid."""

    def test_abstractstate_construction_and_identity(self, cp_name: str, rp_name: str, alias: str):
        """
        AbstractState('REFPROP', alias) must construct without raising (proving
        the alias was resolved to a valid .FLD name) and must report the same
        CAS number as the canonical REFPROP name (proving the correct fluid was
        loaded, not just any fluid).

        We use CAS rather than fluid_param_string('name') because REFPROP
        returns a 12-character truncated chemical name from NAMEdll that does
        not necessarily match the .FLD filename stem.
        """
        ref_cas  = CoolProp.AbstractState("REFPROP", rp_name).fluid_param_string("CAS")
        test_cas = CoolProp.AbstractState("REFPROP", alias).fluid_param_string("CAS")
        assert test_cas == ref_cas, (
            f"AbstractState('REFPROP', {alias!r}) CAS={test_cas!r}, "
            f"expected {ref_cas!r} (canonical: {rp_name!r}, CoolProp: {cp_name!r})"
        )


# ---------------------------------------------------------------------------
# 2. REFPROP-only fluids (no CoolProp alias)
# ---------------------------------------------------------------------------

def _refprop_only_stems() -> list[str]:
    """FLD stems present in REFPROP but absent from CoolProp's fluid library."""
    cp_rp_names = {rp for rp, _, _ in _CP_RP_INDEX.values()}
    return sorted(_fld_stems() - {n.upper() for n in cp_rp_names})


_RP_ONLY = _refprop_only_stems()


@skip
@pytest.mark.parametrize("rp_name", _RP_ONLY, ids=_RP_ONLY)
class TestREFPROPOnlyFluids:
    """
    Fluids that exist in REFPROP but have no CoolProp entry.
    The alias-resolution code must pass these through unchanged so REFPROP
    can load them directly by their native filename.
    """

    def test_abstractstate_constructs_and_has_cas(self, rp_name: str):
        """
        AbstractState must construct without raising (proving the name passes
        through the alias-resolution step unchanged and REFPROP finds the .FLD)
        and must report a non-empty CAS number.  We do not compare the CAS to
        anything from CoolProp because these fluids have no CoolProp entry, and
        we do not compare REFPROP's NAMEdll output to the .FLD stem because
        REFPROP returns a 12-character truncated chemical name that does not
        necessarily match the filename (e.g. '1,3-Butadien' for 13BUTADIENE).
        """
        state = CoolProp.AbstractState("REFPROP", rp_name)
        cas = state.fluid_param_string("CAS")
        assert cas, f"fluid_param_string('CAS') empty for REFPROP-only fluid {rp_name!r}"


# ---------------------------------------------------------------------------
# 3. Path-based fluid specification
# ---------------------------------------------------------------------------

@skip
class TestDirectFLDPath:
    """
    Tests for fluid names that include a REFPROP file extension (.FLD / .fld /
    .PPF / .ppf).

    When a name ends with a known REFPROP extension the alias-resolution step
    strips the extension, resolves the remaining stem through the CoolProp
    alias table, and then lets the normal file-search loop reconstruct the
    correct path.  This means both canonical REFPROP names and CoolProp aliases
    can be suffixed with .FLD and still load correctly.

    Absolute paths (e.g. /opt/refprop/FLUIDS/BUTANE.FLD) are NOT supported:
    the implementation prepends the REFPROP fluids-directory prefix to every
    component, which produces a double-prefixed path for absolute inputs.
    """

    def test_canonical_name_with_fld_extension(self):
        """BUTANE.FLD constructs AbstractState successfully."""
        state = CoolProp.AbstractState("REFPROP", "BUTANE.FLD")
        assert state.fluid_param_string("CAS") == "106-97-8"

    def test_canonical_name_lowercase_fld_extension(self):
        """butane.fld (all lowercase) also works."""
        state = CoolProp.AbstractState("REFPROP", "butane.fld")
        assert state.fluid_param_string("CAS") == "106-97-8"

    def test_alias_with_fld_extension(self):
        """R600.FLD resolves the alias before looking up the file."""
        state = CoolProp.AbstractState("REFPROP", "R600.FLD")
        assert state.fluid_param_string("CAS") == "106-97-8"

    def test_alias_lowercase_mixed_extension(self):
        """r600.fld (all lowercase) resolves correctly via uppercased stem lookup."""
        state = CoolProp.AbstractState("REFPROP", "r600.fld")
        assert state.fluid_param_string("CAS") == "106-97-8"

    def test_canonical_name_alias_identity(self):
        """BUTANE.FLD and R600.FLD load the same fluid."""
        canonical = CoolProp.AbstractState("REFPROP", "BUTANE.FLD").fluid_param_string("CAS")
        alias     = CoolProp.AbstractState("REFPROP", "R600.FLD").fluid_param_string("CAS")
        assert alias == canonical

    def test_unknown_name_with_fld_extension_raises(self):
        """An unresolvable .FLD name raises ValueError."""
        with pytest.raises(ValueError, match="Could not load"):
            CoolProp.AbstractState("REFPROP", "DOES_NOT_EXIST.FLD")

    def test_absolute_path_raises(self):
        """
        Absolute paths are not supported: the implementation prepends the
        fluids-directory prefix, producing an invalid double-prefixed path.
        """
        path = os.path.join(_FLUIDS_DIR, "BUTANE.FLD")
        with pytest.raises(ValueError, match="Could not load"):
            CoolProp.AbstractState("REFPROP", path)

    def test_mixture_alias_and_fld_name(self):
        """Mixing a CoolProp alias with a .FLD-suffixed name in a mixture constructs successfully."""
        CoolProp.AbstractState("REFPROP", "R600&R290.FLD")


# ---------------------------------------------------------------------------
# 4. Mixture access patterns
# ---------------------------------------------------------------------------

@skip
class TestMixtures:
    """Mixtures combining different name forms — verify construction succeeds."""

    def test_both_aliases(self):
        """Both components specified as CoolProp aliases."""
        CoolProp.AbstractState("REFPROP", "R600&R290")

    def test_alias_and_refprop_name(self):
        """First component as alias, second as canonical REFPROP name."""
        CoolProp.AbstractState("REFPROP", "R600&PROPANE")

    def test_refprop_name_and_alias(self):
        """First component as canonical REFPROP name, second as alias."""
        CoolProp.AbstractState("REFPROP", "BUTANE&R290")

    def test_four_component_mixed_names(self):
        """Four-component mixture: alias / canonical / alias / canonical."""
        CoolProp.AbstractState("REFPROP", "R600&R290&R134a&Ethane")

    def test_mixture_with_refprop_only_fluid(self):
        """One REFPROP-only fluid combined with a CoolProp alias."""
        CoolProp.AbstractState("REFPROP", "13BUTADIENE&R290")


# ---------------------------------------------------------------------------
# 5. Pseudopure (.PPF) fluids
# ---------------------------------------------------------------------------

_PPF_CASES = [
    ("Air",   "AIR",   ["air", "AIR", "R729"]),
    ("R410A", "R410A", ["R410a"]),
    ("R407C", "R407C", []),
    ("R507A", "R507A", ["R507a"]),
    ("R404A", "R404A", ["R404a"]),
]
# Filter to those whose CoolProp entry + PPF file both exist
_PPF_CASES_FILTERED = [
    (cp, rp, aliases)
    for cp, rp, aliases in _PPF_CASES
    if rp.upper() in _ppf_stems()
]


@skip
class TestPPFFluids:
    """Pseudopure fluids stored as .PPF; aliases must resolve the same way."""

    @pytest.mark.parametrize(
        "cp_name, rp_name, aliases",
        _PPF_CASES_FILTERED,
        ids=[c[0] for c in _PPF_CASES_FILTERED],
    )
    def test_ppf_aliases_resolve(self, cp_name: str, rp_name: str, aliases: list[str]):
        """Every alias for a PPF fluid must load the same fluid as the canonical REFPROP name."""
        ref_cas = CoolProp.AbstractState("REFPROP", rp_name).fluid_param_string("CAS")
        for alias in [cp_name] + aliases:
            cas = CoolProp.AbstractState("REFPROP", alias).fluid_param_string("CAS")
            assert cas == ref_cas, (
                f"AbstractState('REFPROP', {alias!r}) CAS={cas!r}, expected {ref_cas!r}"
            )


# ---------------------------------------------------------------------------
# 6. Edge cases and negative tests
# ---------------------------------------------------------------------------

@skip
class TestEdgeCases:
    def test_unknown_name_raises_valueerror(self):
        """A completely unknown name must raise ValueError, not crash silently."""
        with pytest.raises(ValueError):
            CoolProp.AbstractState("REFPROP", "DOES_NOT_EXIST_FLUID")

    def test_case_insensitive_alias_uppercase(self):
        """Uppercase alias 'WATER' must load the same fluid as canonical 'Water'."""
        ref_cas = CoolProp.AbstractState("REFPROP", "WATER").fluid_param_string("CAS")
        cas     = CoolProp.AbstractState("REFPROP", "Water").fluid_param_string("CAS")
        assert cas == ref_cas

    def test_case_insensitive_alias_lowercase(self):
        """Lowercase alias 'water' must load the same fluid as canonical 'WATER'."""
        ref_cas = CoolProp.AbstractState("REFPROP", "WATER").fluid_param_string("CAS")
        cas     = CoolProp.AbstractState("REFPROP", "water").fluid_param_string("CAS")
        assert cas == ref_cas


# ---------------------------------------------------------------------------
# 7. Discovery smoke test
# ---------------------------------------------------------------------------

@skip
class TestDiscovery:
    def test_at_least_50_alias_cases_found(self):
        """Sanity check: discovery must produce at least 50 (fluid, alias) pairs."""
        assert len(_ALIAS_CASES) >= 50, (
            f"Only {len(_ALIAS_CASES)} alias cases found; "
            "check that the REFPROP FLUIDS directory is populated."
        )

    def test_known_aliases_present(self):
        """R600, n-Butane, and Water must appear among the discovered aliases."""
        found = {alias for _, _, alias in _ALIAS_CASES}
        for expected in ("R600", "n-Butane", "Water"):
            assert expected in found, f"Alias {expected!r} not discovered"

    def test_refprop_only_fluids_found(self):
        """At least one REFPROP-only fluid (no CoolProp entry) must exist."""
        assert len(_RP_ONLY) >= 1, (
            "Expected at least one REFPROP-only fluid, found none. "
            "Check that the FLUIDS directory contains fluids beyond CoolProp's set."
        )
