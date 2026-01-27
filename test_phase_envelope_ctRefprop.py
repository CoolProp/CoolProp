

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os


def test_phase_envelope_R410A():
    print("Testing phase envelope for R410A.MIX at 278.15 K")
    print("=" * 60)
    
    try:

        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        RP.RPVersion()  


        # Legacy
        r = RP.SETMIXTUREdll("R410A.MIX")
        mix = r[0]
        mw =RP.WMOLdll(mix)
        r = RP.TPRHOdll(250, 100, mix, 2, 10, 10)
        r = RP.THERMdll(250, r.D, mix)
        r = RP.TPFLSHdll(250, 1e2, mix)
        # THERMdll(&_T, &rho_mol_L, &(mole_fractions[0]), &p_kPa, &emol, &hmol, &smol, &cvmol, &cpmol, &w, &hjt);

        r = RP.SATSPLNdll(mix)
        r = RP.SATTdll(278.15, mix, 2)
        print(f"Pv={r[0]:.6f} kPa, Dv={r[1]:.12f} mol/L, Dl={r[2]:.6f} mol/L")
        r = RP.THERMdll(278.15, r.Dv, mix)
        print(f"P={r[0]:.6f} kPa")
        r = RP.PQFLSHdll(1e2, 0.5, mix, 2)
        #print(f"P={r[0]:.6f} kPa")

        #print(f"  Saturation pressure (Q=1): {r.P:.10f} kPa")
        print(f"  Vapor density: {r.D:.10f} mol/m³")
        print(f"  Vapor density: {r.D*mw:.10f} kg/m³")        


        print("done")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_phase_envelope_R410A()


