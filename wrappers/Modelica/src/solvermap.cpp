#include "solvermap.h"
#include "basesolver.h"
#include "testsolver.h"

#if (FLUIDPROP == 1)
#    include "fluidpropsolver.h"
#endif  // FLUIDPROP == 1

#if (COOLPROP == 1)
#    include "coolpropsolver.h"
#endif  // COOLPROP == 1

//! Get a specific solver
/*!
  This function returns the solver for the specified library name, substance name
  and possibly medium name. It creates a new solver if the solver does not already
  exist. When implementing new solvers, one has to add the newly created solvers to
  this function. An error message is generated if the specific library is not supported
  by the interface library.
  @param mediumName Medium name
  @param libraryName Library name
  @param substanceName Substance name
*/
BaseSolver* SolverMap::getSolver(const string& mediumName, const string& libraryName, const string& substanceName) {
    // Get solver key from library and substance name
    string solverKeyString(solverKey(libraryName, substanceName));
    // Check whether solver already exists
    if (_solvers.find(solverKeyString) != _solvers.end()) return _solvers[solverKeyString];
    // Create new solver if it doesn't exist
    // Test solver for compiler setup debugging
    if (libraryName.compare("TestMedium") == 0) _solvers[solverKeyString] = new TestSolver(mediumName, libraryName, substanceName);
#if (FLUIDPROP == 1)
    // FluidProp solver
    else if (libraryName.find("FluidProp") == 0)
        _solvers[solverKeyString] = new FluidPropSolver(mediumName, libraryName, substanceName);
#endif  // FLUIDPROP == 1

#if (COOLPROP == 1)
    // CoolProp solver
    else if (libraryName.find("CoolProp") == 0)
        _solvers[solverKeyString] = new CoolPropSolver(mediumName, libraryName, substanceName);
#endif  // COOLPROP == 1

    else {
        // Generate error message
        char error[100];
        sprintf(error, "Error: libraryName = %s is not supported by any external solver\n", libraryName.c_str());
        errorMessage(error);
    }
    // Return pointer to solver
    return _solvers[solverKeyString];
};

//! Generate a unique solver key
/*!
  This function generates a unique solver key based on the library name and
  substance name.
*/
string SolverMap::solverKey(const string& libraryName, const string& substanceName) {
    // This function returns the solver key and may be changed by advanced users
    return libraryName + "." + substanceName;
}

map<string, BaseSolver*> SolverMap::_solvers;
