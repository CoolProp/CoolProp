import os
import ctREFPROP.ctREFPROP as ct
root  = os.getenv("COOLPROP_REFPROP_ROOT")
RP = ct.REFPROPFunctionLibrary(root)
RP.SETPATHdll(root)
r = RP.REFPROPdll('WATER','PQ','T',RP.MOLAR_BASE_SI,0,0,101325,0,[1.0])
if r.ierr != 0:
    raise ValueError(r.herr)