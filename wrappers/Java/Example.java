// Example for Java
// Ian Bell, 2013

public class Example {
    static {
        System.loadLibrary("CoolProp");
    }
    
    public static void main(String argv[]){
        
        double T, h, p, D;
        int SI,KSI;
        System.out.println("CoolProp version: " + CoolProp.get_global_param_string("version"));
        System.out.println("CoolProp gitrevision: " + CoolProp.get_global_param_string("gitrevision"));
        System.out.println("CoolProp fluids: " + CoolProp.get_global_param_string("FluidsList"));

        System.out.println(" ");
        System.out.println("************ USING EOS *************");
        System.out.println(" ");
        System.out.println("FLUID STATE INDEPENDENT INPUTS");
        //System.out.println("Critical Density Propane: " + CoolProp.Props1SI("Propane", "rhocrit") + "kg/m^3");
        System.out.println("TWO PHASE INPUTS (Pressure)");
        System.out.println("Density of saturated liquid Propane at 101325 kPa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "Propane") + " kg/m^3");
        System.out.println("Density of saturated vapor R290 at 101325 kPa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 1, "R290") + " kg/m^3");
        System.out.println("TWO PHASE INPUTS (Temperature)");
        System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "Propane") + " kg/m^3");
        System.out.println("Density of saturated vapor R290 at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 1, "R290") + " kg/m^3");
        System.out.println("SINGLE PHASE CYCLE (propane)");
        p = CoolProp.PropsSI("P", "T", 300, "D", 1, "Propane");
        h = CoolProp.PropsSI("H", "T", 300, "D", 1, "Propane");
        System.out.println("T,D -> P,H " + 300 + "," + 1 + " --> " + p + "," + h);
        T = CoolProp.PropsSI("T", "P", p, "H", h, "Propane");
        D = CoolProp.PropsSI("D", "P", p, "H", h, "Propane");
        System.out.println("P,H -> T,D " + p + "," + h + " --> " + T + "," + D);

        //~ System.out.println(" ");
        //~ System.out.println("************ USING TTSE ***************");
        //~ System.out.println(" ");
        //~ CoolProp.enable_TTSE_LUT("Propane");
        //~ System.out.println("TWO PHASE INPUTS (Pressure)");
        //~ System.out.println("Density of saturated liquid Propane at 101325 kPa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "Propane") + " kg/m^3");
        //~ System.out.println("Density of saturated vapor R290 at 101325 kPa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 1, "R290") + " kg/m^3");
        //~ System.out.println("TWO PHASE INPUTS (Temperature)");
        //~ System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "Propane") + " kg/m^3");
        //~ System.out.println("Density of saturated vapor R290 at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 1, "R290") + " kg/m^3");
        //~ System.out.println("SINGLE PHASE CYCLE (propane)");
        //~ p = CoolProp.PropsSI("P", "T", 300, "D", 1, "Propane");
        //~ h = CoolProp.PropsSI("H", "T", 300, "D", 1, "Propane");
        //~ System.out.println("T,D -> P,H " + 300 + ","+ 1+ " --> " + p + "," + h);
        //~ T = CoolProp.PropsSI("T", "P", p, "H", h, "Propane");
        //~ D = CoolProp.PropsSI("D", "P", p, "H", h, "Propane");
        //~ System.out.println("P,H -> T,D " + p + "," + h + " --> " + T + "," + D);
        //~ CoolProp.disable_TTSE_LUT("Propane");

        try
        {
            System.out.println(" ");
            System.out.println("************ USING REFPROP ***************");
            System.out.println(" ");
            System.out.println("TWO PHASE INPUTS (Pressure)");
            System.out.println("Density of saturated liquid Propane at 101325 kPa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "REFPROP::Propane") + " kg/m^3");
            System.out.println("TWO PHASE INPUTS (Temperature)");
            System.out.println("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "REFPROP::Propane") + " kg/m^3");
            System.out.println("SINGLE PHASE CYCLE (propane)");
            p = CoolProp.PropsSI("P","T",300,"D",1,"REFPROP::Propane"); 
            h = CoolProp.PropsSI("H","T",300,"D",1,"REFPROP::Propane");
            System.out.println("T,D -> P,H " + 300 + "," + 1 + " --> " + p + "," + h);
            T = CoolProp.PropsSI("T","P",p,"H",h,"REFPROP::Propane"); 
            D = CoolProp.PropsSI("D","P",p,"H",h,"REFPROP::Propane");
            System.out.println("P,H -> T,D " + p + "," + h + " --> " + T + "," + D);
        }
        catch(Exception e)
        {
            System.out.println(" ");
            System.out.println("************ CANT USE REFPROP ************");
            System.out.println(" ");
        }

        System.out.println(" ");
        System.out.println("************ BRINES AND SECONDARY WORKING FLUIDS *************");
        System.out.println(" ");
        System.out.println("Density of 50% (mass) ethylene glycol/water at 300 K, 101325 kPa: " + CoolProp.PropsSI("D", "T", 300, "P", 101325, "EG-50%") + "kg/m^3");
        System.out.println("Viscosity of Therminol D12 at 350 K, 101325 kPa: " + CoolProp.PropsSI("V", "T", 350, "P", 101325, "TD12") + "Pa-s");

        //~ System.out.println(" ");
        //~ System.out.println("************ HUMID AIR PROPERTIES *************");
        //~ System.out.println(" ");
        //~ System.out.println("Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: " + CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5) + " kg_w/kg_da");
        //~ System.out.println("Relative humidity from last calculation: " + CoolProp.HAProps("R", "T", 300, "P", 101.325, "W", CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5)) + "(fractional)");
    }
}