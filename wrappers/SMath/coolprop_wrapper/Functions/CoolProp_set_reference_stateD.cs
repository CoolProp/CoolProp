using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
  class CoolProp_set_reference_stateD : IFunction
  {
    // int set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
    [DllImport(
      "CoolProp_x86", EntryPoint = "set_reference_stateD",
      CharSet = CharSet.Ansi)]
    internal static extern int CoolPropDLLfunc_x86(
      string Ref,
      double T,
      double rho,
      double h0,
      double s0);
    [DllImport(
      "CoolProp_x64", EntryPoint = "set_reference_stateD",
      CharSet = CharSet.Ansi)]
    internal static extern int CoolPropDLLfunc_x64(
      string Ref,
      double T,
      double rho,
      double h0,
      double s0);
    internal static int CoolPropDLLfunc(
      string Ref,
      double T,
      double rho,
      double h0,
      double s0)
    {
      switch (System.IntPtr.Size)
      {
        case 4:
          return CoolPropDLLfunc_x86(Ref, T, rho, h0, s0);
        case 8:
          return CoolPropDLLfunc_x64(Ref, T, rho, h0, s0);
      }
      throw new System.Exception("Unknown platform!");
    }

    Term inf;
    public static int[] Arguments = new [] {5};

    public CoolProp_set_reference_stateD(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    public TermInfo GetTermInfo(string lang)
    {
      string funcInfo = "(FluidName, T, rhomolar, h0, s0) Set the reference state based on a thermodynamic state point specified by temperature and molar density\r\n" +
        "FluidName The name of the fluid\r\n" +
        "T Temperature at reference state [K]\r\n" +
        "rhomolar Density at reference state [mol/m^3]\r\n" +
        "h0 Enthalpy at reference state [J/mol]\r\n" +
        "s0 Entropy at references state [J/mol/K]";

      var argsInfos = new [] {
        new ArgumentInfo(ArgumentSections.String),
        new ArgumentInfo(ArgumentSections.RealNumber),
        new ArgumentInfo(ArgumentSections.RealNumber),
        new ArgumentInfo(ArgumentSections.RealNumber),
        new ArgumentInfo(ArgumentSections.RealNumber),
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
      var T   = coolpropPlugin.GetNumberParam(args[1], ref context);
      var rho = coolpropPlugin.GetNumberParam(args[2], ref context);
      var h0  = coolpropPlugin.GetNumberParam(args[3], ref context);
      var s0  = coolpropPlugin.GetNumberParam(args[4], ref context);
      var Result = CoolPropDLLfunc(Ref, T.obj.ToDouble(), rho.obj.ToDouble(), h0.obj.ToDouble(), s0.obj.ToDouble());
      coolpropPlugin.LogInfo("[INFO ]",
        "Ref = {0} T = {1} rho = {2} h0 = {3} s0 = {4} Result = {5}",
        Ref, T.obj.ToDouble(), rho.obj.ToDouble(), h0.obj.ToDouble(), s0.obj.ToDouble(), Result);
      result = coolpropPlugin.MakeDoubleResult(Result, Unit.unitless);

      return true;
    }
  }
}
