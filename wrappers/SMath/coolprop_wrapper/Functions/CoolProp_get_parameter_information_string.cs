using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_get_parameter_information_string : IFunction
    {
        // long get_parameter_information_string(const char *key, char *Output, int n);
        [DllImport(
            "CoolProp_x86", EntryPoint = "get_parameter_information_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x86(
            string key,
            System.Text.StringBuilder Output,
            int n);
        [DllImport(
            "CoolProp_x64", EntryPoint = "get_parameter_information_string",
            CharSet = CharSet.Ansi)]
        internal static extern long CoolPropDLLfunc_x64(
            string key,
            System.Text.StringBuilder Output,
            int n);
        internal static long CoolPropDLLfunc(
            string key,
            System.Text.StringBuilder Output,
            int n)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(key, Output, n);
                case 8:
                    return CoolPropDLLfunc_x64(key, Output, n);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new [] {2};

        public CoolProp_get_parameter_information_string(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        public TermInfo GetTermInfo(string lang)
        {
            string funcInfo = "(Key, Output) Get a parameter information string\r\n" +
                "Key: A string\r\n" +
                "Output: A string, one of \"IO\", \"short\", \"long\", \"units\"";

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
            var key    = coolpropPlugin.GetStringParam(value.Items[0], context);
            var Output = coolpropPlugin.GetStringParam(value.Items[1], context);
            var output = new System.Text.StringBuilder(Output, 10000);
            var Result = CoolPropDLLfunc(key, output, output.Capacity);
            coolpropPlugin.LogInfo("[INFO]",
                                   "key = {0} output(in) = {1} output(out) = {2} Result = {3}",
                                   key, Output, output.ToString(), Result);
            if (Result != 1)
                coolpropPlugin.CoolPropError();
            result = coolpropPlugin.MakeStringResult(output.ToString());

            return true;
        }
    }
}
