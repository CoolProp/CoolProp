

#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

#include <exception>
#include <iostream>

namespace CoolProp
{

class CoolPropBaseError: public std::exception
{
protected:
	std::string err; // Can be accessed by subclasses since it is protected
public:
	CoolPropBaseError() throw() {};
	~CoolPropBaseError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class NotImplementedError: public CoolPropBaseError
{
public:
	NotImplementedError() throw() {};
	NotImplementedError(std::string errstring) throw(){err=errstring;};
    ~NotImplementedError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class SolutionError: public CoolPropBaseError
{
public:
	SolutionError()throw() {};
	SolutionError(std::string errstring) throw(){err=errstring;};
    ~SolutionError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class ValueError: public CoolPropBaseError
{
public:
	ValueError() throw() {};
	ValueError(std::string errstring) throw(){err=errstring;};
    ~ValueError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class AttributeError: public CoolPropBaseError
{
public:
	AttributeError() throw() {};
	AttributeError(std::string errstring) throw() {err=errstring;};
    ~AttributeError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class OutOfRangeError: public ValueError
{
public:
	OutOfRangeError() throw() {};
	OutOfRangeError(std::string errstring) throw() {err=errstring;};
    ~OutOfRangeError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class WrongFluidError: public ValueError
{
public:
	WrongFluidError() throw() {};
	WrongFluidError(std::string errstring) throw() {err=errstring;};
    ~WrongFluidError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class CompositionError: public ValueError
{
public:
	CompositionError() throw() {};
	CompositionError(std::string errstring) throw() {err=errstring;};
    ~CompositionError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class InputError: public ValueError
{
public:
	InputError() throw() {};
	InputError(std::string errstring) throw() {err=errstring;};
    ~InputError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class NotAvailableError: public ValueError
{
public:
	NotAvailableError() throw() {};
	NotAvailableError(std::string errstring) throw() {err=errstring;};
    ~NotAvailableError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

}; /* namespace CoolProp */
#endif
