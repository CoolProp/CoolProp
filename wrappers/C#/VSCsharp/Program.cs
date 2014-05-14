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
            double T, rho, h;
            T = CoolProp.Props1("R410A","Tcrit");

            Console.Write("The critical temperature of R410A is: " + T + " K\n");
            
            rho = CoolProp.Props("D",'T',298.15,'P',101.325,"Air");
            Console.Write("Density of air at STP is "+rho+" kg/m3\n");
            
            h = CoolProp.Props("D",'T',275,'Q',1.0,"R290");
            Console.Write("The saturated vapor enthalpy of Propane at 275 K is: "+h+" kJ/kg\n");

            Console.Write("Enter to quit");
            Console.ReadLine();
        }
    }
}
