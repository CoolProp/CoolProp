﻿using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
  class CoolProp_saturation_ancillary : IFunction
  {
    // double saturation_ancillary(const char *fluid_name, const char *output, int Q, const char *input, double value);
    [DllImport(
      "CoolProp_x86", EntryPoint = "saturation_ancillary",
      CharSet = CharSet.Ansi)]
    internal static extern double CoolPropDLLfunc_x86(
      string fluid_name,
      string output,
      int Q,
      string input,
      double value);
    [DllImport(
      "CoolProp_x64", EntryPoint = "saturation_ancillary",
      CharSet = CharSet.Ansi)]
    internal static extern double CoolPropDLLfunc_x64(
      string fluid_name,
      string output,
      int Q,
      string input,
      double value);
    internal static double CoolPropDLLfunc(
      string fluid_name,
      string output,
      int Q,
      string input,
      double value)
    {
      switch (System.IntPtr.Size) {
        case 4:
          return CoolPropDLLfunc_x86(fluid_name, output, Q, input, value);
        case 8:
          return CoolPropDLLfunc_x64(fluid_name, output, Q, input, value);
      }
      throw new System.Exception("Unknown platform!");
    }

    Term inf;
    public static int[] Arguments = new[] { 5 };

    public CoolProp_saturation_ancillary(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    TermInfo IFunction.GetTermInfo(string lang)
    {
      string funcInfo = "(fluid_name, output, Q, input, value) Extract a value from the saturation ancillary\r\n" +
        "fluid_name The name of the fluid to be used - HelmholtzEOS backend only\r\n" +
        "output The desired output variable (\"P\" for instance for pressure)\r\n" +
        "Q The quality, 0 or 1\r\n" +
        "input The input variable (\"T\")\r\n" +
        "value The input value";

      var argsInfos = new [] {
        new ArgumentInfo(ArgumentSections.String),
        new ArgumentInfo(ArgumentSections.String),
        new ArgumentInfo(ArgumentSections.RealNumber),
        new ArgumentInfo(ArgumentSections.String),
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
      var fluid_name = coolpropPlugin.GetStringParam(args[0], ref context);
      var output     = coolpropPlugin.GetStringParam(args[1], ref context);
      var Q     = (int)coolpropPlugin.GetNumberParam(args[2], ref context).obj.ToDouble();
      var input      = coolpropPlugin.GetStringParam(args[3], ref context);
      var value      = coolpropPlugin.GetNumberParam(args[4], ref context);
      Unit.match(input, value, context);
      var Result = CoolPropDLLfunc(fluid_name, output, Q, input, value.obj.ToDouble());
      coolpropPlugin.LogInfo("[INFO ]",
        "fluid_name = {0}, output = {1}, Q = {2}, input = {3}, value = {4}, Result = {5}",
        fluid_name, output, Q, input, value.obj.ToDouble(), Result);
      result = coolpropPlugin.MakeDoubleResult(Result, Unit.Find(output));

      return true;
    }
  }
}
