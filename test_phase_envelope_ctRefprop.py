

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os


def test_phase_envelope_R410A():
    print("Testing mass quality calculations with R410A.MIX at 278.15 K and 0.5 quality (ctRefprop)")
    print("=" * 49)
    
    try:

        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        RP.RPVersion()  

        r = RP.SETMIXTUREdll("R410A.MIX")
        mix = r[0]
        nc = len([m for m in mix if m > 0.0])
        mw =RP.WMOLdll(mix) # kg/kmol

        r = RP.TQFLSHdll(278.15, 0.5, mix, 1)
        print(f"  Saturation pressure (molar basis): {r.P*1000:.12f} Pa")
        print(f"  Mole density: {r.D:.12f} kmol/m3")
        print(f"  Mass density: {r.D*mw:.12f} kg/m3")
        print(f"  Molar quality: {0.5:.12f} -")
        _r = RP.QMASSdll(0.5, r.x, r.y)
        print(f"  Mass quality: {_r.qkg:.12f} -")
        print(f"  Mole fractions: {' '.join(f'{m:.12f}' for m in mix[:nc])}")
        _r = RP.XMASSdll(mix)
        print(f"  Mass fractions: {' '.join(f'{m:.12f}' for m in _r.xkg[:nc])}")

        r = RP.TQFLSHdll(278.15, 0.5, mix, 2)
        print(f"  Saturation pressure (mass basis): {r.P*1000:.12f} Pa")
        print(f"  Mole density: {r.D:.12f} kmol/m3")
        print(f"  Mass density: {r.D*mw:.12f} kg/m3")
        ykg = RP.XMASSdll(r.y).xkg
        xkg = RP.XMASSdll(r.x).xkg
        _r = RP.QMOLEdll(0.5, xkg, ykg)
        print(f"  Molar quality: {_r.qmol:.12f} -")
        print(f"  Mass quality: {0.5:.12f} -")

        r = RP.PQFLSHdll(1e2, 0.5, mix, 1)
        print(f"  Saturation temperature (molar basis): {r.T:.12f} K")
        print(f"  Mole density: {r.D:.12f} kmol/m3")
        print(f"  Mass density: {r.D*mw:.12f} kg/m3")
        print(f"  Molar quality: {0.5:.12f} -")
        _r = RP.QMASSdll(0.5, r.x, r.y)
        print(f"  Mass quality: {_r.qkg:.12f} -")

        r = RP.PQFLSHdll(1e2, 0.5, mix, 2)
        print(f"  Saturation temperature (mass basis): {r.T:.12f} K")
        print(f"  Mole density: {r.D:.12f} kmol/m3")
        print(f"  Mass density: {r.D*mw:.12f} kg/m3")
        ykg = RP.XMASSdll(r.y).xkg
        xkg = RP.XMASSdll(r.x).xkg
        _r = RP.QMOLEdll(0.5, xkg, ykg)
        print(f"  Molar quality: {_r.qmol:.12f} -")
        print(f"  Mass quality: {0.5:.12f} -")



    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_phase_envelope_R410A()



    # RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    # RP.SETPATHdll(os.environ['RPPREFIX'])
    # RP.RPVersion()  


    # # Legacy
    # r = RP.SETMIXTUREdll("R410A.MIX")
    # mix = r[0]
    # mw =RP.WMOLdll(mix)
    # r = RP.TPRHOdll(250, 100, mix, 2, 10, 10)
    # r = RP.THERMdll(250, r.D, mix)
    # r = RP.TPFLSHdll(250, 1e2, mix)
    # # THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);

    # r = RP.SATSPLNdll(mix)
    # r = RP.SATTdll(278.15, mix, 2)
    # print(f"Pv={r[0]:.6f} kPa, Dv={r[1]:.12f} mol/L, Dl={r[2]:.6f} mol/L")
    # r = RP.THERMdll(278.15, r.Dv, mix)
    # print(f"P={r[0]:.6f} kPa")
    # r = RP.PQFLSHdll(1e2, 0.5, mix, 2)
    # #print(f"P={r[0]:.6f} kPa")

    # #print(f"  Saturation pressure (Q=1): {r.P:.10f} kPa")
    # print(f"  Vapor density: {r.D:.10f} mol/m³")
    # print(f"  Vapor density: {r.D*mw:.10f} kg/m³")        
