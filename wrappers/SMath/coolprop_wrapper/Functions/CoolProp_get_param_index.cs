using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_get_param_index : IFunction
    {
        // long get_param_index(const char *param);
        [DllImport(
            "CoolProp_x86", EntryPoint = "get_param_index",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x86(
            string param);
        [DllImport(
            "CoolProp_x64", EntryPoint = "get_param_index",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x64(
            string param);
        internal static long CoolPropDLLfunc(
            string param)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(param);
                case 8:
                    return CoolPropDLLfunc_x64(param);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new[] { 1 };

        public CoolProp_get_param_index(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        TermInfo IFunction.GetTermInfo(string lang)
        {
            string funcInfo = "(Name) Returns the index of a parameter\r\n" +
                "Name: The parameter name, one of \"Tcrit\", \"D\", \"H\", etc...";

            var argsInfos = new [] {
                new ArgumentInfo(ArgumentSections.String),
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
            var Result = CoolPropDLLfunc(param);
            coolpropPlugin.LogInfo("[INFO]", "param = {0}, Result = {1}", param, Result);
            result = coolpropPlugin.MakeDoubleResult(Result, Unit.unitless);

            return true;
        }
    }
}
