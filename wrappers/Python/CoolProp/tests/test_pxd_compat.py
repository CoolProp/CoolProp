"""
Backward-compatibility tests: verify Cython .pxd files remain usable for
PDSim-style cimport, and that the Cython module is still accessible alongside
the new pybind11 module.
"""
import os
import sys
import pytest


def test_pxd_files_installed():
    import CoolProp
    pkg_dir = os.path.dirname(CoolProp.__file__)
    pxd_files = [f for f in os.listdir(pkg_dir) if f.endswith('.pxd')]
    assert pxd_files, f"No .pxd files found in {pkg_dir}"


def test_cython_module_still_accessible():
    from CoolProp.CoolProp import AbstractState
    assert hasattr(AbstractState, '__pyx_vtable__'), \
        "CoolProp.CoolProp.AbstractState is not a Cython type (missing __pyx_vtable__)"


def test_pybind11_module_accessible():
    pytest.importorskip('CoolProp._CoolProp_pybind11',
                        reason="_CoolProp_pybind11 not built")
    from CoolProp._CoolProp_pybind11 import AbstractState
    AS = AbstractState('HEOS', 'Water')
    assert AS is not None


def test_pdsim_cython_step():
    Cython = pytest.importorskip('Cython', reason="Cython not installed")
    pdsim_src = os.environ.get('PDSIM_SOURCE_DIR', '')
    if not pdsim_src:
        pytest.skip("PDSIM_SOURCE_DIR env var not set")
    from Cython.Build import cythonize
    pyx_files = [
        os.path.join(pdsim_src, f)
        for f in os.listdir(pdsim_src)
        if f.endswith('.pyx')
    ]
    if not pyx_files:
        pytest.skip(f"No .pyx files found in {pdsim_src}")
    cythonize(pyx_files, include_path=[os.path.dirname(__import__('CoolProp').__file__)])


def test_pdsim_full_build():
    pdsim_src = os.environ.get('PDSIM_SOURCE_DIR', '')
    if not pdsim_src:
        pytest.skip("PDSIM_SOURCE_DIR env var not set")
    import subprocess
    result = subprocess.run(
        [sys.executable, 'setup.py', 'build_ext', '--inplace'],
        cwd=pdsim_src,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, \
        f"PDSim build failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
