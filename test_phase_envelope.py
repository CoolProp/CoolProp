"""
Test phase envelope building for R410A.MIX at 278.15 K
"""
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

def test_phase_envelope_R410A():
    print("Testing phase envelope for R410A.MIX at 278.15 K")
    print("=" * 60)
    
    try:
        # Create abstract state for R410A.MIX
        AS: AbstractState = AbstractState("REFPROP", "R410A.MIX")
        AS.build_phase_envelope("")

      
        target_T = 278.15
        print(f"\n\nTesting property calculation at T = {target_T} K:")
        try:
            #AS.update(CP.QT_INPUTS, 0.5, target_T)  
            AS.update(CP.PQ_INPUTS, 1e5, 0.5)  
            
            print(f"  Saturation pressure (Q=1): {AS.p():.10f} Pa")
            print(f"  Vapor density: {AS.rhomolar():.10f} mol/m³")
            print(f"  Vapor density: {AS.rhomass():.10f} kg/m³")

        except Exception as e:
            print(f"  Error during property calculation: {e}")
        
        return True
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_phase_envelope_R410A()
    exit(0 if success else 1)
