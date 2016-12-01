using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_HAProps : IFunction
    {
        // double HAPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Name3, double Prop3);
        [DllImport(
            "CoolProp_x86", EntryPoint = "HAPropsSI",
            CharSet = CharSet.Ansi)]
        internal static extern double CoolPropDLLfunc_x86(
            string Output,
            string Name1,
            double Prop1,
            string Name2,
            double Prop2,
            string Name3,
            double Prop3);
        [DllImport(
            "CoolProp_x64", EntryPoint = "HAPropsSI",
            CharSet = CharSet.Ansi)]
        internal static extern double CoolPropDLLfunc_x64(
            string Output,
            string Name1,
            double Prop1,
            string Name2,
            double Prop2,
            string Name3,
            double Prop3);
        internal static double CoolPropDLLfunc(
            string Output,
            string Name1,
            double Prop1,
            string Name2,
            double Prop2,
            string Name3,
            double Prop3)
        {
            switch (System.IntPtr.Size) {
                case 4:
                    return CoolPropDLLfunc_x86(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3);
                case 8:
                    return CoolPropDLLfunc_x64(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new[] { 7 };

        public CoolProp_HAProps(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        TermInfo IFunction.GetTermInfo(string lang)
        {
            string funcInfo = "(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3) Return a humid air property\r\n" +
                "Output: The output parameter, one of \"T\", \"D\", \"H\", etc...\r\n" +
                "Name1: The first state variable name, one of \"T\", \"D\", \"H\", etc...\r\n" +
                "Prop1: The first state variable value\r\n" +
                "Name2: The second state variable name, one of \"T\", \"D\", \"H\", etc...\r\n" +
                "Prop2: The second state variable value\r\n" +
                "Name3: The third state variable name, one of \"T\", \"D\", \"H\", etc...\r\n" +
                "Prop3: The third state variable value";

            var argsInfos = new [] {
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.RealNumber),
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.RealNumber),
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.RealNumber),
            };

            return new TermInfo(inf.Text,
                                inf.Type,
                                funcInfo,
                                FunctionSections.Unknown,
                                true,
                                argsInfos);
        }

        public bool TryEvaluateExpression(Entry value, SMath.Math.Store context, out Entry result)
        {
            // Possible inputs:
            // "Omega"="HumRat"="W"        = Humidity ratio
            // "Tdp"="T_dp"="DewPoint"="D" = Dewpoint temperature
            // "Twb"="T_wb"="WetBulb"="B"  = Wet bulb temperature
            // "Enthalpy"="H"              = Enthalpy
            // "Entropy"="S"               = Entropy
            // "RH"="RelHum"="R"           = Relative humidity
            // "Tdb"="T_db"="T"            = Dry-bulb temperature
            // "P"                         = Pressure
            // "V"="Vda"                   = Volume of dry air
            // "mu"="Visc"="M"             = Viscosity
            // "k"="Conductivity"="K"      = Conductivity

            // Output:
            // "Vda"="V"[m^3/kg_da]
            // "Vha"[m^3/kg_ha]
            // "Y"[mol_w/mol]
            // "Hda"="H"
            // "Hha"[kJ/kg_ha]
            // "S"="Entropy"[kJ/kg_da]
            // "C"="cp"[kJ/kg_da]
            // "Cha"="cp_ha"[kJ/kg_da]
            // "Tdp"="D"[K]
            // "Twb"="T_wb"="WetBulb"="B"[K]
            // "Omega"="HumRat"="W"
            // "RH"="RelHum"="R"
            // "mu"="Visc"="M"
            // "k"="Conductivity"="K"

            var Output = coolpropPlugin.GetStringParam(value.Items[0], context);
            var Name1 = coolpropPlugin.GetStringParam(value.Items[1], context);
            var Prop1 = coolpropPlugin.GetNumberParam(value.Items[2], context);
            Unit.matchHA(Name1, Prop1, context);
            var Name2 = coolpropPlugin.GetStringParam(value.Items[3], context);
            var Prop2 = coolpropPlugin.GetNumberParam(value.Items[4], context);
            Unit.matchHA(Name2, Prop2, context);
            var Name3 = coolpropPlugin.GetStringParam(value.Items[5], context);
            var Prop3 = coolpropPlugin.GetNumberParam(value.Items[6], context);
            Unit.matchHA(Name3, Prop3, context);
            var Result = CoolPropDLLfunc(Output, Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), Name3, Prop3.obj.ToDouble());
            coolpropPlugin.LogInfo("[INFO]",
                                   "Output = {0}, Name1 = {1}, Prop1 = {2}, Name2 = {3}, Prop2 = {4}, Name3 = {5}, Prop3 = {6}, Result = {7}",
                                   Output, Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), Name3, Prop3.obj.ToDouble(), Result);
            result = coolpropPlugin.MakeDoubleResult(Result, Unit.FindHA(Output));

            return true;
        }
    }
}
