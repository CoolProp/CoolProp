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

    public CoolProp_set_reference_stateS(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    public TermInfo GetTermInfo(string lang)
    {
      string funcInfo = "(FluidName, reference_state) Set the reference state based on a string representation\r\n" +
        "FluidName The name of the fluid\r\n" +
        "reference_state The reference state to use, one of \"IIR\", \"ASHRAE\", \"NBP\", \"DEF\", \"RESET\"";

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

    bool IFunction.ExpressionEvaluation(Term root, Term[][] args, ref SMath.Math.Store context, ref Term[] result)
    {
      var Ref = coolpropPlugin.GetStringParam(args[0], ref context);
      var reference_state = coolpropPlugin.GetStringParam(args[1], ref context);
      var Result = CoolPropDLLfunc(Ref, reference_state);
      coolpropPlugin.LogInfo("[INFO ]", "Ref = {0}, reference_state = {1} Result = {2}", Ref, reference_state, Result);
      result = coolpropPlugin.MakeDoubleResult(Result, Unit.unitless);

      return true;
    }
  }
}
