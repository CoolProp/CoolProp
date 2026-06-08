

#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

#include <exception>
#include <string>

namespace CoolProp {

class CoolPropBaseError : public std::exception
{
   public:
    enum ErrCode
    {
        eNotImplemented,
        eSolution,
        eAttribute,
        eOutOfRange,
        eValue,
        eWrongFluid,
        eComposition,
        eInput,
        eNotAvailable,
        eHandle,
        eKey,
        eUnableToLoad,
        eDirectorySize,
        eMultipleSolutions
    };
    CoolPropBaseError(const std::string& err, ErrCode code) throw() : m_code(code), m_err(err) {}
    ~CoolPropBaseError() throw() = default;
    const char* what() const throw() override {
        return m_err.c_str();
    }
    ErrCode code() {
        return m_code;
    }

   private:
    ErrCode m_code;
    std::string m_err;
};

template <CoolPropBaseError::ErrCode errcode>
class CoolPropError : public CoolPropBaseError
{
   public:
    CoolPropError(const std::string& err = "", ErrCode ecode = errcode) throw() : CoolPropBaseError(err, ecode) {}
};

using NotImplementedError = CoolPropError<CoolPropBaseError::eNotImplemented>;
using SolutionError = CoolPropError<CoolPropBaseError::eSolution>;
using AttributeError = CoolPropError<CoolPropBaseError::eAttribute>;
using OutOfRangeError = CoolPropError<CoolPropBaseError::eOutOfRange>;
using ValueError = CoolPropError<CoolPropBaseError::eValue>;
using KeyError = CoolPropError<CoolPropBaseError::eKey>;
using HandleError = CoolPropError<CoolPropBaseError::eHandle>;
using UnableToLoadError = CoolPropError<CoolPropBaseError::eUnableToLoad>;
using DirectorySizeError = CoolPropError<CoolPropBaseError::eDirectorySize>;

// ValueError specializations
template <CoolPropBaseError::ErrCode errcode>
class ValueErrorSpec : public ValueError
{
   public:
    ValueErrorSpec(const std::string& err = "", ErrCode ecode = errcode) throw() : ValueError(err, ecode) {}
};

using WrongFluidError = ValueErrorSpec<CoolPropBaseError::eWrongFluid>;
using CompositionError = ValueErrorSpec<CoolPropBaseError::eComposition>;
using InputError = ValueErrorSpec<CoolPropBaseError::eInput>;
using NotAvailableError = ValueErrorSpec<CoolPropBaseError::eNotAvailable>;
/// Thrown when a saturation flash input maps to more than one temperature on
/// the saturation curve (e.g. water saturated-vapor enthalpy at a given h has
/// two T-roots). The caller should use update_with_guesses with a guess.T to
/// pick a branch. See GitHub #2773.
using MultipleSolutionsError = ValueErrorSpec<CoolPropBaseError::eMultipleSolutions>;

}; /* namespace CoolProp */
#endif
