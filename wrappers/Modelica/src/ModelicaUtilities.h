#ifndef MODELICA_UTILITIES_H
#define MODELICA_UTILITIES_H

/* Utility functions which can be called by external Modelica functions.

  These are defined in the Modelica Language specification and implemented
  by all Modelica tools
*/

extern void ModelicaMessage(const char* string);
/*
Output the message string (no format control).
*/

extern void ModelicaFormatMessage(const char* string, ...);
/*
Output the message under the same format control as the C-function printf.
  */

extern void ModelicaError(const char* string);
/*
Output the error message string (no format control). This function
never returns to the calling function, but handles the error
similarly to an assert in the Modelica code.
*/

extern void ModelicaFormatError(const char* string, ...);
/*
Output the error message under the same format control as the C-function
printf. This function never returns to the calling function,
but handles the error similarly to an assert in the Modelica code.
*/

extern char* ModelicaAllocateString(size_t len);
/*
Allocate memory for a Modelica string which is used as return
argument of an external Modelica function. Note, that the storage
for string arrays (= pointer to string array) is still provided by the
calling program, as for any other array. If an error occurs, this
function does not return, but calls "ModelicaError".
*/

extern char* ModelicaAllocateStringWithErrorReturn(size_t len);
/*
Same as ModelicaAllocateString, except that in case of error, the
function returns 0. This allows the external function to close files
and free other open resources in case of error. After cleaning up
resources use ModelicaError or ModelicaFormatError to signal
the error.
*/

#endif
