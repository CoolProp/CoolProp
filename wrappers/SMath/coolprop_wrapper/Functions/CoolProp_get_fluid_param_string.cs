using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_get_fluid_param_string : IFunction
    {
        // long get_fluid_param_string(const char *fluid, const char *param, char *Output, int n);
        [DllImport(
            "CoolProp_x86", EntryPoint = "get_fluid_param_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x86(
            string fluid,
            string param,
            System.Text.StringBuilder Output,
            int n);
        [DllImport(
            "CoolProp_x64", EntryPoint = "get_fluid_param_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x64(
            string fluid,
            string param,
            System.Text.StringBuilder Output,
            int n);
        internal static long CoolPropDLLfunc(
            string fluid,
            string param,
            System.Text.StringBuilder Output,
            int n)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(fluid, param, Output, n);
                case 8:
                    return CoolPropDLLfunc_x64(fluid, param, Output, n);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new [] {2};

        public CoolProp_get_fluid_param_string(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        public TermInfo GetTermInfo(string lang)
        {
            string funcInfo = "(FluidName, ParamName) Get a string for a value from a fluid\r\n" +
                "FluidName: The name of the fluid that is part of CoolProp, for instance \"n-Propane\"\r\n" +
                "ParamName: A string, can be in one of \"aliases\", \"CAS\", \"CAS_number\", \"ASHRAE34\", \"REFPROPName\", \"REFPROP_name\"";

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
            var fluid = coolpropPlugin.GetStringParam(value.Items[0], context);
            var param = coolpropPlugin.GetStringParam(value.Items[1], context);
            var output = new System.Text.StringBuilder(10000);
            var Result = CoolPropDLLfunc(fluid, param, output, output.Capacity);
            coolpropPlugin.LogInfo("[INFO]", "fluid = {0} param = {1} output = {2} Result = {3}", fluid, param, output.ToString(), Result);
            if (Result != 1)
                coolpropPlugin.CoolPropError();
            result = coolpropPlugin.MakeStringResult(output.ToString());

            return true;
        }
    }
}
