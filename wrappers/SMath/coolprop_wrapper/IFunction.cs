namespace coolprop_wrapper
{
    interface IFunction
    {
        SMath.Manager.Term Info { get; }

        SMath.Manager.TermInfo GetTermInfo(string lang);

        bool TryEvaluateExpression(
            SMath.Manager.Entry value,
            SMath.Math.Store context,
            out SMath.Manager.Entry result);
    }
}
