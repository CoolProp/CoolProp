using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
    class CoolProp_set_reference_stateS : IFunction
    {
        // int set_reference_stateS(const char *Ref, const char *reference_state);
        [DllImport(
            "CoolProp_x86", EntryPoint = "set_reference_stateS",
            CharSet = CharSet.Ansi)]
        internal static extern int CoolPropDLLfunc_x86(
            string Ref,
            string reference_state);
        [DllImport(
            "CoolProp_x64", EntryPoint = "set_reference_stateS",
            CharSet = CharSet.Ansi)]
        internal static extern int CoolPropDLLfunc_x64(
            string Ref,
            string reference_state);
        internal static int CoolPropDLLfunc(
            string Ref,
            string reference_state)
        {
            switch (System.IntPtr.Size)
            {
                case 4:
                    return CoolPropDLLfunc_x86(Ref, reference_state);
                case 8:
                    return CoolPropDLLfunc_x64(Ref, reference_state);
            }
            throw new System.Exception("Unknown platform!");
        }

        Term inf;
        public static int[] Arguments = new [] {2};

        public CoolProp_set_reference_stateS(int argsCount)
        {
            inf = new Term(this.GetType().Name, TermType.Function, argsCount);
        }

        Term IFunction.Info { get { return inf; } }

        public TermInfo GetTermInfo(string lang)
        {
            string funcInfo = "(FluidName, ReferenceState) Set the reference state based on a string representation\r\n" +
                "FluidName: The name of the fluid\r\n" +
                "ReferenceState: The reference state to use, one of \"IIR\", \"ASHRAE\", \"NBP\", \"DEF\", \"RESET\"";

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
            var Ref = coolpropPlugin.GetStringParam(value.Items[0], context);
            var reference_state = coolpropPlugin.GetStringParam(value.Items[1], context);
            var Result = CoolPropDLLfunc(Ref, reference_state);
            coolpropPlugin.LogInfo("[INFO]", "Ref = {0}, reference_state = {1} Result = {2}", Ref, reference_state, Result);
            result = coolpropPlugin.MakeDoubleResult(Result, Unit.unitless);

            return true;
        }
    }
}
