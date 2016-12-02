using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_get_global_param_string : IFunction
    {
        // long get_global_param_string(const char *param, char *Output, int n);
        [DllImport(
            "CoolProp_x86", EntryPoint = "get_global_param_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x86(
            string param,
            System.Text.StringBuilder Output,
            int n);
        [DllImport(
            "CoolProp_x64", EntryPoint = "get_global_param_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x64(
            string param,
            System.Text.StringBuilder Output,
            int n);
        internal static bool CoolPropDLLfunc(string param, out string resultStr)
        {
            var output = new System.Text.StringBuilder(10000);
            long Result = 0;
            switch (System.IntPtr.Size)
            {
                case 4:
                    Result = CoolPropDLLfunc_x86(param, output, output.Capacity);
                    break;
                case 8:
                    Result = CoolPropDLLfunc_x64(param, output, output.Capacity);
                    break;
                default:
                    throw new System.Exception("Unknown platform!");
            }
            coolpropPlugin.LogInfo("[INFO]", "param = {0} output = {1} Result = {2}", param, resultStr = output.ToString(), Result);
            return (Result == 1);
        }

        Term inf;
        public static int[] Arguments = new [] {1};

        public CoolProp_get_global_param_string(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        TermInfo IFunction.GetTermInfo(string lang)
        {
            string funcInfo = "(ParamName) Get a globally-defined string\r\n" +
                "ParamName: A string, one of \"version\", \"errstring\", \"warnstring\", \"gitrevision\", \"FluidsList\", \"fluids_list\", \"parameter_list\",\"predefined_mixtures\"";

            var argsInfos = new [] {
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
            var param = coolpropPlugin.GetStringParam(value.Items[0], context);
            string resultStr;
            if (!CoolPropDLLfunc(param, out resultStr))
                coolpropPlugin.CoolPropError();
            result = coolpropPlugin.MakeStringResult(resultStr);
            return true;
        }
    }
}
