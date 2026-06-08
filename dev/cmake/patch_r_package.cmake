# Rewrite the SWIG-generated R proxy's PACKAGE= references so they match the
# renamed shared library (CoolProp -> CoolPropR, see #1674). R resolves the
# PACKAGE argument of .Call() against the dyn.load basename, not the SWIG
# %module name, so without this the renamed R wrapper cannot find its routines.
file(READ "${RFILE}" _contents)
string(REPLACE "PACKAGE='CoolProp'" "PACKAGE='CoolPropR'" _contents "${_contents}")
string(REPLACE "PACKAGE=\"CoolProp\"" "PACKAGE=\"CoolPropR\"" _contents "${_contents}")
file(WRITE "${RFILE}" "${_contents}")
