"""Sanity checks on the raw reference-data loaders.

Companions to DATA_AUDIT.md: the audit found exactly one ordering hazard in
the corpus (the HFE-7100 tables store temperature descending), and the
loaders now sort every grid by its own axis values instead of trusting file
order. These tests pin that behavior for every data-backed fluid.

Needs numpy and scipy (the loaders import both); skips cleanly without either.
"""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")  # SecCoolFluids imports it transitively at module load

from CPIncomp.SecCoolFluids import SecCoolSolutionData

PROPERTY_IDS = ["Rho", "Cp", "Mu", "Cond"]


def _loaded_axes(obj, dataID):
    import os

    if not os.path.isfile(obj.getFile(dataID)):
        return None, None  # not every fluid carries every property table
    arr = obj.getFromFile(dataID)
    if arr.shape == (1, 1):
        return None, None
    return arr[1:, 0], arr[0, 1:]


@pytest.fixture(scope="module")
def seccool_fluids():
    # factory() instantiates and file-loads every SecCool fluid (and prints
    # progress); build it once for the module rather than per test.
    return SecCoolSolutionData.factory()


def test_hfe7100_descending_source_loads_ascending(seccool_fluids):
    # The one real ordering hazard found by the audit: raw HFE-7100 files
    # store T from +64 down to -80 degC.
    fluid = next((o for o in seccool_fluids if getattr(o, "sFile", "") == "HFE-7100"), None)
    assert fluid is not None, "HFE-7100 not found among the SecCool fluids"
    for dataID in PROPERTY_IDS:
        T, _x = _loaded_axes(fluid, dataID)
        assert T is not None, dataID
        assert np.all(np.diff(T) > 0), (dataID, T)


def test_all_seccool_grids_load_with_increasing_axes(seccool_fluids):
    for fluid in seccool_fluids:
        if not hasattr(fluid, "sFile"):
            continue
        if type(fluid).getFromFile is not SecCoolSolutionData.getFromFile:
            # Subclasses with their own loaders (SecCoolIceData) read csv
            # tables of which only Hfusion is in production use; the unused
            # Cond/Mu csvs are latin-1 encoded and not loadable as-is (see
            # DATA_AUDIT.md).
            continue
        for dataID in PROPERTY_IDS:
            T, x = _loaded_axes(fluid, dataID)
            if T is None:
                continue  # no file for this property
            Tf, xf = T[np.isfinite(T)], x[np.isfinite(x)]
            assert np.all(np.diff(Tf) > 0), (fluid.name, dataID, Tf)
            assert np.all(np.diff(xf) >= 0), (fluid.name, dataID, xf)


def test_all_toplevel_grids_load_with_increasing_temperature():
    import glob
    import os

    from CPIncomp.DataObjects import DigitalData

    data_dir = os.path.join(os.path.dirname(__file__), "CPIncomp", "data")
    loader = DigitalData()
    for path in sorted(glob.glob(os.path.join(data_dir, "*.txt"))):
        base = os.path.basename(path)
        name, _, dataID = base.rpartition("_")
        loader.name = name
        arr = loader.getFromFile(dataID[:-4])
        T = arr[1:, 0]
        T = T[np.isfinite(T)]
        assert np.all(np.diff(T) > 0), (base, T)


def test_seccool_ice_hfusion_grid_loads_with_increasing_axes(seccool_fluids):
    # SecCoolIceData overrides getFromFile, so the test above skips it and its
    # sortGridAxes call went unexercised. Only the Hfusion table is in
    # production use (the Cond/Mu csvs for these fluids are latin-1 encoded and
    # not loadable as-is -- see DATA_AUDIT.md), so pin that one directly.
    #
    # These tables happen to ship ascending already, so this asserts the
    # contract for the real production path rather than exercising the sort
    # itself; test_sort_grid_axes_orders_a_shuffled_grid below does that.
    from CPIncomp.SecCoolFluids import SecCoolIceData

    ice = [o for o in seccool_fluids if isinstance(o, SecCoolIceData)]
    assert ice, "no SecCoolIceData fluids found"
    for fluid in ice:
        T, x = _loaded_axes(fluid, "Hfusion")
        assert T is not None, fluid.name
        assert np.all(np.diff(T) > 0), (fluid.name, "temperature axis not ascending", T)
        if x is not None and np.size(x) > 1:
            assert np.all(np.diff(x) > 0), (fluid.name, "composition axis not ascending", x)


def test_sort_grid_axes_orders_a_shuffled_grid():
    # Direct check on the shared helper, with a grid that is genuinely out of
    # order in both directions -- the shipped tables are mostly ascending
    # already, so without this the sorting logic could regress unnoticed.
    from CPIncomp.BaseObjects import IncompressibleFitter

    grid = np.array([
        [0.0, 0.30, 0.10, 0.20],   # composition axis, shuffled
        [300.0, 13.0, 11.0, 12.0],
        [100.0, 3.0, 1.0, 2.0],    # temperature axis, descending
        [200.0, 8.0, 6.0, 7.0],
    ])
    out = IncompressibleFitter.sortGridAxes(grid.copy())

    assert np.all(np.diff(out[1:, 0]) > 0), out[1:, 0]
    assert np.all(np.diff(out[0, 1:]) > 0), out[0, 1:]
    # Values must travel with their axes, not just be re-sorted independently.
    expected = np.array([
        [0.0, 0.10, 0.20, 0.30],
        [100.0, 1.0, 2.0, 3.0],
        [200.0, 6.0, 7.0, 8.0],
        [300.0, 11.0, 12.0, 13.0],
    ])
    assert np.array_equal(out, expected), out
