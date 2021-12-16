import CoolProp.CoolProp as cp

def main():
    
    fluid = "helium"
    
    heos = cp.AbstractState("HEOS", fluid)
    bicubic = cp.AbstractState("BICUBIC&HEOS", fluid)
    
    heos.update(cp.PT_INPUTS, 1e5, 300.0)
    bicubic.update(cp.DmolarUmolar_INPUTS, heos.rhomolar(), heos.umolar())
    
    print(bicubic.p())
    print(bicubic.T())
    
if __name__== "__main__":
    
    main()