using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_Props1 : IFunction
    {
        // double Props1SI(const char *FluidName, const char* Output);
        [DllImport(
            "CoolProp_x86", EntryPoint = "Props1SI",
            CharSet = CharSet.Ansi)]
        internal static extern double CoolPropDLLfunc_x86(
            string FluidName,
            string Output);
        [DllImport(
            "CoolProp_x64", EntryPoint = "Props1SI",
            CharSet = CharSet.Ansi)]
        internal static extern double CoolPropDLLfunc_x64(
            string FluidName,
            string Output);
        internal static double CoolPropDLLfunc(
            string FluidName,
            string Output)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(FluidName, Output);
                case 8:
                    return CoolPropDLLfunc_x64(FluidName, Output);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new[] { 2 };

        public CoolProp_Props1(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        TermInfo IFunction.GetTermInfo(string lang)
        {
            string funcInfo = "(FluidName, Output) Return a value that does not depends on the thermodynamic state\r\n" +
                "FluidName: The fluid name\r\n" +
                "Output: The output parameter, one of \"Tcrit\", \"D\", \"H\", etc...";

            var argsInfos = new [] {
                new ArgumentInfo(ArgumentSections.String),
                new ArgumentInfo(ArgumentSections.String)
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
            string FluidName = coolpropPlugin.GetStringParam(value.Items[0], context),
            Output    = coolpropPlugin.GetStringParam(value.Items[1], context);
            var Result = CoolPropDLLfunc(FluidName, Output);
            coolpropPlugin.LogInfo("[INFO]", "FluidName = {0}, Output = {1}, Result = {2}", FluidName, Output, Result);
            var ResUnit = Unit.Find(Output);
            // Props1SI may take parameters in either order; so, Output may actually be in FluidName
            if (ResUnit == Unit.unitless)
                ResUnit = Unit.Find(FluidName);
            result = coolpropPlugin.MakeDoubleResult(Result, ResUnit);

            return true;
        }
    }
}
