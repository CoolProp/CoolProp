using System;
using System.Collections.Generic;
using System.Linq;
using System.Text; 

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {
            AbstractState State = AbstractState.factory("HEOS","Water");
            State.update((int)input_pairs.PT_INPUTS, 1e5, 300);
            double hmol = State.hmolar();
            Console.Write("Hmol: " + hmol + " J/kg" + "\n");
            
            double T, h, p, D;
            //Console.Write("CoolProp version: " + CoolProp.get_global_param_string("version") + "\n");
            //Console.Write("CoolProp gitrevision: " + CoolProp.get_global_param_string("gitrevision") + "\n");
            //Console.Write("CoolProp fluids: " + CoolProp.get_global_param_string("FluidsList") + "\n");

            Console.Write(" " + "\n");
            Console.Write("************ USING EOS *************" + "\n");
            Console.Write(" " + "\n");
            Console.Write("FLUID STATE INDEPENDENT INPUTS" + "\n");
            //Console.Write("Critical Density Propane: " + CoolProp.Props1SI("Propane", "rhocrit") + " kg/m^3" + "\n");
            Console.Write("TWO PHASE INPUTS (Pressure)" + "\n");
            Console.Write("Density of saturated liquid Propane at 101325 Pa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "Propane") + " kg/m^3" + "\n");
            Console.Write("Density of saturated vapor R290 at 101325 Pa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 1, "R290") + " kg/m^3" + "\n");
            Console.Write("TWO PHASE INPUTS (Temperature)" + "\n");
            Console.Write("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "Propane") + " kg/m^3" + "\n");
            Console.Write("Density of saturated vapor R290 at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 1, "R290") + " kg/m^3" + "\n");
            Console.Write("SINGLE PHASE CYCLE (propane)" + "\n");
            p = CoolProp.PropsSI("P", "T", 300, "D", 1, "Propane");
            h = CoolProp.PropsSI("H", "T", 300, "D", 1, "Propane");
            Console.Write("T,D -> P,H " + 300 + "," + 1 + " --> " + p + "," + h + "\n");
            T = CoolProp.PropsSI("T", "P", p, "H", h, "Propane");
            D = CoolProp.PropsSI("D", "P", p, "H", h, "Propane");
            Console.Write("P,H -> T,D " + p + "," + h + " --> " + T + "," + D + "\n");

            //~ Console.Write(" " + "\n");
            //~ Console.Write("************ USING TTSE ***************" + "\n");
            //~ Console.Write(" " + "\n");
            //~ //CoolProp.enable_TTSE_LUT("Propane");
            //~ Console.Write("TWO PHASE INPUTS (Pressure)" + "\n");
            //~ Console.Write("Density of saturated liquid Propane at 101325 Pa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "Propane") + " kg/m^3" + "\n");
            //~ Console.Write("Density of saturated vapor R290 at 101325 Pa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 1, "R290") + " kg/m^3" + "\n");
            //~ Console.Write("TWO PHASE INPUTS (Temperature)" + "\n");
            //~ Console.Write("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "Propane") + " kg/m^3" + "\n");
            //~ Console.Write("Density of saturated vapor R290 at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 1, "R290") + " kg/m^3" + "\n");
            //~ Console.Write("SINGLE PHASE CYCLE (propane)" + "\n");
            //~ p = CoolProp.PropsSI("P", "T", 300, "D", 1, "Propane");
            //~ h = CoolProp.PropsSI("H", "T", 300, "D", 1, "Propane");
            //~ Console.Write("T,D -> P,H " + 300 + ","+ 1+ " --> " + p + "," + h + "\n");
            //~ T = CoolProp.PropsSI("T", "P", p, "H", h, "Propane");
            //~ D = CoolProp.PropsSI("D", "P", p, "H", h, "Propane");
            //~ Console.Write("P,H -> T,D " + p + "," + h + " --> " + T + "," + D + "\n");
            //CoolProp.disable_TTSE_LUT("Propane");

            try
            {
                Console.Write(" " + "\n");
                Console.Write("************ USING REFPROP ***************" + "\n");
                Console.Write(" " + "\n");
                Console.Write("TWO PHASE INPUTS (Pressure)" + "\n");
                Console.Write("Density of saturated liquid Propane at 101325 Pa: " + CoolProp.PropsSI("D", "P", 101325, "Q", 0, "REFPROP::Propane") + " kg/m^3" + "\n");
                Console.Write("TWO PHASE INPUTS (Temperature)" + "\n");
                Console.Write("Density of saturated liquid Propane at 300 K: " + CoolProp.PropsSI("D", "T", 300, "Q", 0, "REFPROP::Propane") + " kg/m^3" + "\n");
                Console.Write("SINGLE PHASE CYCLE (propane)" + "\n");
                p = CoolProp.PropsSI("P","T",300,"D",1,"REFPROP::Propane"); 
                h = CoolProp.PropsSI("H","T",300,"D",1,"REFPROP::Propane");
                Console.Write("T,D -> P,H " + 300 + "," + 1 + " --> " + p + "," + h + "\n");
                T = CoolProp.PropsSI("T","P",p,"H",h,"REFPROP::Propane"); 
                D = CoolProp.PropsSI("D","P",p,"H",h,"REFPROP::Propane");
                Console.Write("P,H -> T,D " + p + "," + h + " --> " + T + "," + D + "\n");
            }
            catch
            {
                Console.Write(" " + "\n");
                Console.Write("************ CANT USE REFPROP ************" + "\n");
                Console.Write(" " + "\n");
            }

            Console.Write(" " + "\n");
            Console.Write("************ BRINES AND SECONDARY WORKING FLUIDS *************" + "\n");
            Console.Write(" " + "\n");
            Console.Write("Density of 50% (mass) ethylene glycol/water at 300 K, 101325 Pa: " + CoolProp.PropsSI("D", "T", 300, "P", 101325, "EG-50%") + " kg/m^3" + "\n");
            Console.Write("Viscosity of Therminol D12 at 350 K, 101325 kPa: " + CoolProp.PropsSI("V", "T", 350, "P", 101325, "INCOMP::TD12") + " Pa-s" + "\n");

            Console.Write(" " + "\n");
            Console.Write("************ HUMID AIR PROPERTIES *************" + "\n");
            Console.Write(" " + "\n");
            //Console.Write("Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: " + CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5) + " kg_w/kg_da" + "\n");
            //Console.Write("Relative humidity from last calculation: " + CoolProp.HAProps("R", "T", 300, "P", 101.325, "W", CoolProp.HAProps("W", "T", 300, "P", 101.325, "R", 0.5)) + "(fractional)" + "\n");

            //Console.Write("Enter to quit");
            //Console.ReadLine();
        }
    }
}
