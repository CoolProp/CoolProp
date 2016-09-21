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

    public CoolProp_Props1(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    TermInfo IFunction.GetTermInfo(string lang)
    {
      string funcInfo = "(FluidName, Output) Return a value that does not depend on the thermodynamic state\r\n" +
                        "FluidName The fluid name\r\n" +
                        "Output The output parameter, one of \"Tcrit\",\"D\",\"H\",etc.";

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
      string FluidName = coolpropPlugin.GetStringParam(args[0], ref context),
             Output    = coolpropPlugin.GetStringParam(args[1], ref context);
      var Result = CoolPropDLLfunc(FluidName, Output);
      coolpropPlugin.LogInfo("[INFO ]", "FluidName = {0}, Output = {1}, Result = {2}", FluidName, Output, Result);
      var ResUnit = Unit.Find(Output);
      // Props1SI may take parameters in either order; so, Output may actually be in FluidName
      if (ResUnit == Unit.unitless)
        ResUnit = Unit.Find(FluidName);
      result = coolpropPlugin.MakeDoubleResult(Result, ResUnit);

      return true;
    }
  }
}
