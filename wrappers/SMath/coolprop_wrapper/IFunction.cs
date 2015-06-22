namespace coolprop_wrapper
{
  interface IFunction
  {
    SMath.Manager.Term Info { get; }

    SMath.Manager.TermInfo GetTermInfo(string lang);

    bool ExpressionEvaluation(
      SMath.Manager.Term root,
      SMath.Manager.Term[][] args,
      ref SMath.Math.Store context,
      ref SMath.Manager.Term[] result);
  }
}
