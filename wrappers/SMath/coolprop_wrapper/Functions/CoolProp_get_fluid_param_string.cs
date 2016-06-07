﻿using System.Runtime.InteropServices;
using SMath.Manager;

namespace coolprop_wrapper.Functions
{
  class CoolProp_get_fluid_param_string : IFunction
  {
    // long get_fluid_param_string(const char *fluid, const char *param, char *Output, int n);
    [DllImport(
      "CoolProp.x86.dll", EntryPoint = "get_fluid_param_string",
      CharSet = CharSet.Ansi)]
    internal static extern long CoolPropDLLfunc_x86(
      string fluid,
      string param,
      System.Text.StringBuilder Output,
      int n);
    [DllImport(
      "CoolProp.x64.dll", EntryPoint = "get_fluid_param_string",
      CharSet = CharSet.Ansi)]
    internal static extern long CoolPropDLLfunc_x64(
      string fluid,
      string param,
      System.Text.StringBuilder Output,
      int n);
    internal static long CoolPropDLLfunc(
      string fluid,
      string param,
      System.Text.StringBuilder Output,
      int n)
    {
      switch (System.IntPtr.Size)
      {
        case 4:
          return CoolPropDLLfunc_x86(fluid, param, Output, n);
        case 8:
          return CoolPropDLLfunc_x64(fluid, param, Output, n);
      }
      throw new System.Exception("Unknown platform!");
    }

    Term inf;
    public static int[] Arguments = new [] {2};

    public CoolProp_get_fluid_param_string(int childCount)
    {
      inf = new Term(this.GetType().Name, TermType.Function, childCount);
    }

    Term IFunction.Info { get { return inf; } }

    public TermInfo GetTermInfo(string lang)
    {
      string funcInfo = "(FluidName, ParamName) Get a string for a value from a fluid\r\n" +
        "FluidName The name of the fluid that is part of CoolProp, for instance \"n-Propane\"\r\n" +
        "ParamName A string, can be in one of \"aliases\", \"CAS\", \"CAS_number\", \"ASHRAE34\", \"REFPROPName\", \"REFPROP_name\"";

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
      var fluid = coolpropPlugin.GetStringParam(args[0], ref context);
      var param = coolpropPlugin.GetStringParam(args[1], ref context);
      var output = new System.Text.StringBuilder(10000);
      var Result = CoolPropDLLfunc(fluid, param, output, output.Capacity);
      coolpropPlugin.LogInfo("[INFO ]", "fluid = {0} param = {1} output = {2} Result = {3}", fluid, param, output.ToString(), Result);
      if (Result != 1)
          coolpropPlugin.CoolPropError();
      result = coolpropPlugin.MakeStringResult(output.ToString());

      return true;
    }
  }
}
