using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_set_reference_stateD : IFunction
    {
        // int set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
        [DllImport(
            "CoolProp_x86", EntryPoint = "set_reference_stateD",
            CharSet = CharSet.Ansi)]
        internal static extern int CoolPropDLLfunc_x86(
            string Ref,
            double T,
            double rho,
            double h0,
            double s0);
        [DllImport(
            "CoolProp_x64", EntryPoint = "set_reference_stateD",
            CharSet = CharSet.Ansi)]
        internal static extern int CoolPropDLLfunc_x64(
            string Ref,
            double T,
            double rho,
            double h0,
            double s0);
        internal static int CoolPropDLLfunc(
            string Ref,
            double T,
            double rho,
            double h0,
            double s0)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(Ref, T, rho, h0, s0);
                case 8:
                    return CoolPropDLLfunc_x64(Ref, T, rho, h0, s0);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new [] {5};

        public CoolProp_set_reference_stateD(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        public TermInfo GetTermInfo(string lang)
        {
            string funcInfo = "(FluidName, T, RhoMolar, h0, s0) Set the reference state based on a thermodynamic state point specified by temperature and molar density\r\n" +
                "FluidName: The name of the fluid\r\n" +
                "T: Temperature at reference state [K]\r\n" +
                "RhoMolar: Density at reference state [mol/m^3]\r\n" +
                "h0: Enthalpy at reference state [J/mol]\r\n" +
                "s0: Entropy at references state [J/mol/K]";

            var argsInfos = new [] {
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.RealNumber),
                new ArgumentInfo(ArgumentSections.RealNumber),
                new ArgumentInfo(ArgumentSections.RealNumber),
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
            var Ref = coolpropPlugin.GetStringParam(value.Items[0], context);
            var T   = coolpropPlugin.GetNumberParam(value.Items[1], context);
            var rho = coolpropPlugin.GetNumberParam(value.Items[2], context);
            var h0  = coolpropPlugin.GetNumberParam(value.Items[3], context);
            var s0  = coolpropPlugin.GetNumberParam(value.Items[4], context);
            var Result = CoolPropDLLfunc(Ref, T.obj.ToDouble(), rho.obj.ToDouble(), h0.obj.ToDouble(), s0.obj.ToDouble());
            coolpropPlugin.LogInfo("[INFO]",
                                   "Ref = {0} T = {1} rho = {2} h0 = {3} s0 = {4} Result = {5}",
                                   Ref, T.obj.ToDouble(), rho.obj.ToDouble(), h0.obj.ToDouble(), s0.obj.ToDouble(), Result);
            result = coolpropPlugin.MakeDoubleResult(Result, Unit.unitless);

            return true;
        }
    }
}
