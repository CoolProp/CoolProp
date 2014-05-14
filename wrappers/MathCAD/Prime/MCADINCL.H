#ifndef _MCADINCL_H_
#define _MCADINCL_H_

#include <windows.h>

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus
                           
// complex scalar type
typedef struct tagCOMPLEXSCALAR {
    double real;
    double imag;
} COMPLEXSCALAR;
    
// this is the complex scalar type received from mathcad
typedef const COMPLEXSCALAR * const LPCCOMPLEXSCALAR; 
// this is the complex scalar type that should be returned to mathcad
typedef COMPLEXSCALAR * const LPCOMPLEXSCALAR;
    

// complex array type
typedef struct tagCOMPLEXARRAY  {
    unsigned int rows;
    unsigned int cols;
    double **hReal; // hReal[cols][rows],  == NULL when the real part is zero 
    double **hImag; // hImag[cols][rows],  == NULL when the imaginary part is zero
} COMPLEXARRAY;
    
// readonly complex array type
typedef struct tagReadOnlyCOMPLEXARRAY  {
    const unsigned int rows;
    const unsigned int cols;
    const double *const *const hReal; // hReal[cols][rows],  == NULL when the real part is zero 
    const double *const *const hImag; // hImag[cols][rows],  == NULL when the imaginary part is zero
} ReadOnlyCOMPLEXARRAY;
    
// this is the complex array type received from mathcad
typedef const ReadOnlyCOMPLEXARRAY * const LPCCOMPLEXARRAY;    
// this is the complex array type that should be returned to mathcad
typedef COMPLEXARRAY * const LPCOMPLEXARRAY;


typedef struct tagMCSTRING {
	char *str;
}MCSTRING;

typedef const MCSTRING * const LPCMCSTRING;
typedef MCSTRING * const LPMCSTRING;


// types to be used in declaration of the function's
// arguments and of the return value
#define COMPLEX_SCALAR  1
#define COMPLEX_ARRAY   2
#define STRING			8


//
// File name variables. These are passed as const char *pointers
// if the string doesn't look like it has a path in it then 
// the current working directory will be prepended to the string 
// before it is passed to the user function.
// 
//

// your function will be passed a const char * pointer
#define INFILE			13
// an OUTFILE is like an INFILE except it allows you
// to put your function on the left side of a := like
// the WRITEPRN() builtin
#define OUTFILE			14

// use this structure to create a function
#define MAX_ARGS        10

typedef LRESULT (* LPCFUNCTION ) ( void * const, const void * const, ... );    

typedef struct tagFUNCTIONINFO {
    char *  lpstrName;
    char *  lpstrParameters; 
    char *  lpstrDescription;
    LPCFUNCTION lpfnMyCFunction;
    long unsigned int returnType;
    unsigned int nArgs;
    long unsigned int argType[MAX_ARGS];
} FUNCTIONINFO;

const void * CreateUserFunction( HINSTANCE, FUNCTIONINFO * );
    
BOOL CreateUserErrorMessageTable(   HINSTANCE, 
                                    unsigned int nErrorMessages,
                                    char  * ErrorMessageTable[] );
                                       
// memory management routines
char * MathcadAllocate( unsigned int size );
void MathcadFree( char * address );


// array allocation -- should be used to allocate
// return array
BOOL    MathcadArrayAllocate(   COMPLEXARRAY * const, 
                                unsigned int rows,  
                                unsigned int cols,
                                BOOL allocateReal, 
                                BOOL allocateImag );
// should be used to free ( in case of an error )
// Mathcad allocated return array
void MathcadArrayFree( COMPLEXARRAY * const );
    
    
// this routine can be used to find out 
// whether the user has attempted to interrupt
// Mathcad 
// this routine slows down the execution -- so use judiciously
BOOL isUserInterrupted( void );
    
#ifdef __cplusplus
}
#endif // __cplusplus
    
#endif // _MCADINCL_H_  
