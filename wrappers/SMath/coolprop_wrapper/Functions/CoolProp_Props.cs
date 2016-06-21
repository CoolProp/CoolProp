using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
  class CoolProp_Props : IFunction
  {
    //  double PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
    [DllImport(
      "CoolProp_x86", EntryPoint = "PropsSI",
      CharSet = CharSet.Ansi)]
    internal static extern double CoolPropDLLfunc_x86(
      string Output,
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName);
    [DllImport(
      "CoolProp_x64", EntryPoint = "PropsSI",
      CharSet = CharSet.Ansi)]
    internal static extern double CoolPropDLLfunc_x64(
      string Output,
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName);
    internal static double CoolPropDLLfunc(
      string Output,
      string Name1,
      double Prop1,
      string Name2,
      double Prop2,
      string FluidName)
    {
      switch (System.IntPtr.Size) {
        case 4:
          return CoolPropDLLfunc_x86(Output, Name1, Prop1, Name2, Prop2, FluidName);
        case 8:
          return CoolPropDLLfunc_x64(Output, Name1, Prop1, Name2, Prop2, FluidName);
      }
      throw new System.Exception("Unknown platform!");
    }

    Term inf;
    public static int[] Arguments = new[] { 6 };

    public CoolProp_Props(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    TermInfo IFunction.GetTermInfo(string lang)
    {
      string funcInfo = "(Output, Name1, Prop1, Name2, Prop2, FluidName) Return a value that depends on the thermodynamic state\r\n" +
                        "Output The output parameter, one of \"T\",\"D\",\"H\",etc.\r\n" +
                        "Name1 The first state variable name, one of \"T\",\"D\",\"H\",etc.\r\n" +
                        "Prop1 The first state variable value\r\n" +
                        "Name2 The second state variable name, one of \"T\",\"D\",\"H\",etc.\r\n" +
                        "Prop2 The second state variable value\r\n" +
                        "FluidName The fluid name";

      var argsInfos = new [] {
        new ArgumentInfo(ArgumentSections.String),
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
      var Output = coolpropPlugin.GetStringParam(args[0], ref context);
      var Name1 = coolpropPlugin.GetStringParam(args[1], ref context);
      var Prop1 = coolpropPlugin.GetNumberParam(args[2], ref context);
      Unit.match(Name1, Prop1, context);
      var Name2 = coolpropPlugin.GetStringParam(args[3], ref context);
      var Prop2 = coolpropPlugin.GetNumberParam(args[4], ref context);
      Unit.match(Name2, Prop2, context);
      var FluidName = coolpropPlugin.GetStringParam(args[5], ref context);
      var Result = CoolPropDLLfunc(Output, Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), FluidName);
      coolpropPlugin.LogInfo("[INFO ]",
        "Output = {0}, Name1 = {1}, Prop1 = {2}, Name2 = {3}, Prop2 = {4}, FluidName = {5}, Result = {6}",
        Output, Name1, Prop1.obj.ToDouble(), Name2, Prop2.obj.ToDouble(), FluidName, Result);
      result = coolpropPlugin.MakeDoubleResult(Result, Unit.Find(Output));

      return true;
    }
  }
}
