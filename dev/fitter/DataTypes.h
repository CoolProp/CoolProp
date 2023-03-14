#ifndef FITTER_DATATYPES_H
#define FITTER_DATATYPES_H

class ExperimentalDataPoint
{
   public:
    EOSFitter* EOS;
    double T,     /// The temperature [K]
      rho,        /// The density [mol/m^3]
      p,          /// The pressure [Pa]
      tau,        /// The reciprocal reduced temperature [-]
      delta,      /// The reduced density [-]
      log_tau,    /// The natural logarithm of the reciprocal reduced temperature
      log_delta,  /// The natural logarithm of the reduced density
      variance;   /// The total variance of the datapoint
    virtual double residual(const std::vector<double>& n) = 0;
    double sum_squares(const std::vector<double>& n) {
        return pow(residual(n) / variance, (int)2);
    }
};

class NonlinearExperimentalDataPoint : public ExperimentalDataPoint
{
   public:
    virtual double a_0(const std::vector<double>& n) = 0;
    virtual double a_i(int i) = 0;
    double residual(const std::vector<double>& n);
};

class LinearExperimentalDataPoint : public ExperimentalDataPoint
{
   public:
    virtual double a_0(void) = 0;
    virtual double a_i(int i) = 0;
    double residual(const std::vector<double>& n);
};

class PressureDataPoint : public LinearExperimentalDataPoint
{
   public:
    PressureDataPoint(EOSFitter* EOS, double T, double rho, double p, double variance);
    /// The part that does not depend on the coefficients
    double a_0(void);
    /// The part that multiplies the coefficients
    double a_i(int i);
};

class SpecificHeatCPDataPoint : public NonlinearExperimentalDataPoint
{
   public:
    double cp, cp_over_R;
    SpecificHeatCPDataPoint(EOSFitter* EOS, double T, double rho, double cp, double variance);
    /// The part that does not depend on the coefficients
    /// Here it requires that the coefficients be passed in to calculate the precorrelation factor
    double a_0(const std::vector<double>& n);
    /// The part that multiplies the coefficients
    double a_i(int i);
};

#endif