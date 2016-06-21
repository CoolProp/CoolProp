using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
  class CoolProp_Phase : IFunction
  {
    // long PhaseSI(const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref, char *phase, int n);
    [DllImport(
      "CoolProp_x86", EntryPoint="PhaseSI",
      CharSet = CharSet.Ansi)]
    internal static extern long CoolPropDLLfunc_x86(
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName,
      System.Text.StringBuilder phase,
      int n);
    [DllImport(
      "CoolProp_x64", EntryPoint = "PhaseSI",
      CharSet = CharSet.Ansi)]
    internal static extern long CoolPropDLLfunc_x64(
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName,
      System.Text.StringBuilder phase,
      int n);
    internal static long CoolPropDLLfunc(
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName,
      System.Text.StringBuilder phase,
      int n)
    {
      switch (System.IntPtr.Size)
      {
        case 4:
          return CoolPropDLLfunc_x86(Name1, Prop1, Name2, Prop2, FluidName, phase, n);
        case 8:
          return CoolPropDLLfunc_x64(Name1, Prop1, Name2, Prop2, FluidName, phase, n);
      }
      throw new System.Exception("Unknown platform!");
    }

    Term inf;
    public static int[] Arguments = new [] {5};

    public CoolProp_Phase(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    TermInfo IFunction.GetTermInfo(string lang)
    {
      string funcInfo = "(Name1, Prop1, Name2, Prop2, FluidName) Return a string representation of the phase\r\n" +
                        "Name1 The first state variable name, one of \"T\",\"D\",\"H\",etc.\r\n" +
                        "Prop1 The first state variable value\r\n" +
                        "Name2 The second state variable name, one of \"T\",\"D\",\"H\",etc.\r\n" +
                        "Prop2 The second state variable value\r\n" +
                        "FluidName The fluid name";

      var argsInfos = new [] {
        new ArgumentInfo(ArgumentSections.String),
        new ArgumentInfo(ArgumentSections.RealNumber),
        new ArgumentInfo(ArgumentSections.String),
        new ArgumentInfo(ArgumentSections.RealNumber),
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
      var Name1 = coolpropPlugin.GetStringParam(args[0], ref context);
      var Prop1 = coolpropPlugin.GetNumberParam(args[1], ref context);
      Unit.match(Name1, Prop1, context);
      var Name2 = coolpropPlugin.GetStringParam(args[2], ref context);
      var Prop2 = coolpropPlugin.GetNumberParam(args[3], ref context);
      Unit.match(Name2, Prop2, context);
      var FluidName = coolpropPlugin.GetStringParam(args[4], ref context);
      var phase = new System.Text.StringBuilder(10000);
      var Result = CoolPropDLLfunc(Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), FluidName, phase, phase.Capacity);
      coolpropPlugin.LogInfo("[INFO ]",
        "Name1 = {0}, Prop1 = {1}, Name2 = {2}, Prop2 = {3}, FluidName = {4}, phase = {5}, Result = {6}",
        Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), FluidName, phase.ToString(), Result);
      if (Result != 1)
          coolpropPlugin.CoolPropError();
      result = coolpropPlugin.MakeStringResult(phase.ToString());

      return true;
    }
  }
}
