

#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

#include <exception>
#include <iostream>

namespace CoolProp
{

class CoolPropBaseError: public std::exception
{
public:
    enum ErrCode { eNotImplemented, eSolution, eAttribute, eOutOfRange, eValue, eWrongFluid, eComposition, eInput, eNotAvailable };
    CoolPropBaseError(const std::string &err, ErrCode code) throw() : m_code(code), m_err(err) {}
    ~CoolPropBaseError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
    ErrCode code() { return m_code; }
private:
    ErrCode m_code;
    std::string m_err;
};

template <CoolPropBaseError::ErrCode errcode>
class CoolPropError : public CoolPropBaseError {
public:
    CoolPropError(const std::string &err = "", ErrCode ecode = errcode) throw() : CoolPropBaseError(err, ecode) {}
};

typedef CoolPropError<CoolPropBaseError::eNotImplemented> NotImplementedError;
typedef CoolPropError<CoolPropBaseError::eSolution> SolutionError;
typedef CoolPropError<CoolPropBaseError::eAttribute> AttributeError;
typedef CoolPropError<CoolPropBaseError::eOutOfRange> OutOfRangeError;
typedef CoolPropError<CoolPropBaseError::eValue> ValueError;

// ValueError specializations
template <CoolPropBaseError::ErrCode errcode>
class ValueErrorSpec : public ValueError {
public:
    ValueErrorSpec(const std::string &err = "", ErrCode ecode = errcode) throw() : ValueError(err, ecode) {}
};

typedef ValueErrorSpec<CoolPropBaseError::eWrongFluid> WrongFluidError;
typedef ValueErrorSpec<CoolPropBaseError::eComposition> CompositionError;
typedef ValueErrorSpec<CoolPropBaseError::eInput> InputError;
typedef ValueErrorSpec<CoolPropBaseError::eNotAvailable> NotAvailableError;

}; /* namespace CoolProp */
#endif
