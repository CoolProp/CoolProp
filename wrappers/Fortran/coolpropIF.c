
#include <string.h>
#include <stdio.h>

typedef int bool;
#define false 0
#define true 1

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#  define CONVENTION __declspec(dllexport)
#endif

#ifndef CONVENTION
#  define CONVENTION
#endif

double CONVENTION Props1 (char * FluidName, char * Output);
double CONVENTION Props  (char * Output,char Name1, double Prop1, char Name2, double Prop2, char * FluidName);
double CONVENTION PropsSI(char * Output,char Name1, double Prop1, char Name2, double Prop2, char * FluidName);
double CONVENTION DerivTerms(char * Term, double T, double rho, char * FluidName);
int CONVENTION set_reference_stateS(char * FluidName, char * reference_state);
bool CONVENTION enable_TTSE_LUT(char *FluidName);
bool CONVENTION isenabled_TTSE_LUT(char *FluidName);
bool CONVENTION disable_TTSE_LUT(char *FluidName);
int CONVENTION set_TTSE_mode(char *FluidName, char * Value);


void CONVENTION propsmod_(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char *FluidName, double *resu)
{
  /* return Props(std::string(Output), Name1, Prop1, Name2, Prop2, std::string(FluidName)); */
  *resu = Props(Output, *Name1, *Prop1, *Name2, *Prop2, FluidName);
}

/* Create pointer-only functions for Fortran avoiding underscores for f77 */
double CONVENTION props1_(char *FluidName, char *Output)
{
  /* return Props1(std::string(FluidName), std::string(Output)); */
  return Props1(FluidName, Output);
}

double CONVENTION props_(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char *FluidName)
{
  /* return Props(std::string(Output), Name1, Prop1, Name2, Prop2, std::string(FluidName)); */
  return Props(Output, *Name1, *Prop1, *Name2, *Prop2, FluidName);
}

double CONVENTION propssi_(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char *FluidName)
{
  /* return Props(std::string(Output), Name1, Prop1, Name2, Prop2, std::string(FluidName)); */
  return PropsSI(Output, *Name1, *Prop1, *Name2, *Prop2, FluidName);
}

double CONVENTION derivterms_(char *Term, double *T, double *rho, char *FluidName)
{
  /* return DerivTerms(std::string(Term), T, rho, std::string(FluidName)); */
  return DerivTerms(Term, *T, *rho, FluidName);
}

int CONVENTION setreferencestates_(char *FluidName, char *reference_state)
{
  /* return set_reference_stateS(std::string(FluidName), std::string(reference_state)); */
  return set_reference_stateS(FluidName, reference_state);
}
/* int set_reference_stateS(std::string FluidName, std::string reference_state); */

/*  Enable the TTSE */
int CONVENTION enablettselut_(char *FluidName)
{
  if (enable_TTSE_LUT(FluidName)) return 1;
  else return 0;
}
/* bool CONVENTION enable_TTSE_LUT(char *FluidName); */

/*  Check if TTSE is enabled */
int CONVENTION isenabledttselut_(char *FluidName)
{
  if (isenabled_TTSE_LUT(FluidName)) return 1;
  else return 0;
}
/* bool CONVENTION isenabled_TTSE_LUT(char *FluidName); */

/*  Disable the TTSE */
int CONVENTION disablettselut_(char *FluidName)
{
  if (disable_TTSE_LUT(FluidName)) return 1;
  else return 0;
}
/* bool CONVENTION disable_TTSE_LUT(char *FluidName); */

/*  Set the TTSE mode (normal or bicubic) */
int CONVENTION setttsemode_(char *FluidName, char *Value)
{
  return set_TTSE_mode(FluidName, Value);
}
/* int CONVENTION set_TTSE_mode(char *FluidName, char * Value); */



