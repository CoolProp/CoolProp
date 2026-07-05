"""Sanity checks on the raw reference-data loaders.

Companions to DATA_AUDIT.md: the audit found exactly one ordering hazard in
the corpus (the HFE-7100 tables store temperature descending), and the
loaders now sort every grid by its own axis values instead of trusting file
order. These tests pin that behavior for every data-backed fluid.

Needs numpy (the loaders do); skips cleanly without it.
"""

import numpy as np
import pytest

pytest.importorskip("numpy")

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


def test_hfe7100_descending_source_loads_ascending():
    # The one real ordering hazard found by the audit: raw HFE-7100 files
    # store T from +64 down to -80 degC.
    fluid = next(o for o in SecCoolSolutionData.factory() if getattr(o, "sFile", "") == "HFE-7100")
    for dataID in PROPERTY_IDS:
        T, _x = _loaded_axes(fluid, dataID)
        assert T is not None, dataID
        assert np.all(np.diff(T) > 0), (dataID, T)


def test_all_seccool_grids_load_with_increasing_axes():
    for fluid in SecCoolSolutionData.factory():
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
