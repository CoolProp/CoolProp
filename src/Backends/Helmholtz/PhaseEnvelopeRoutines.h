#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

class PhaseEnvelopeRoutines{
    public:
    /** \brief Build the phase envelope
     *
     * @param HEOS The HelmholtzEOSMixtureBackend instance to be used
     */
    static void build(HelmholtzEOSMixtureBackend &HEOS);
    
    /** \brief Refine the phase envelope, adding points in places that are sparse
     *
     * @param HEOS The HelmholtzEOSMixtureBackend instance to be used
     */
    static void refine(HelmholtzEOSMixtureBackend &HEOS);
    
    /** \brief Finalize the phase envelope and calculate maxima values, critical point, etc.
     * 
     * @param HEOS The HelmholtzEOSMixtureBackend instance to be used
     */
    static void finalize(HelmholtzEOSMixtureBackend &HEOS);
    
    /** \brief Determine which indices bound a given value
     * 
     * If you provide pressure for instance, it will return each of the indices 
     * that bound crossings in the pressure versus rhov curve.  Thus this information
     * can be used to determine whether another input is "inside" or "outside" the phase
     * boundary.
     * 
     * @param HEOS The HelmholtzEOSMixtureBackend instance to be used
     * @param iInput The key for the variable type that is to be checked
     * @param value The value associated with iInput
     */
    static std::vector<std::pair<std::size_t, std::size_t> > find_intersections(HelmholtzEOSMixtureBackend &HEOS, parameters iInput, CoolPropDbl value);
    
    /** \brief Determine whether a pair of inputs is inside or outside the phase envelope
     * 
     * @param HEOS The HelmholtzEOSMixtureBackend instance to be used
     * @param iInput1 The key for the first input
     * @param value1 The value of the first input
     * @param iInput2 The key for the second input
     * @param value2 The value of the second input
     * @param iclosest The index of the phase envelope for the closest point
     * @param closest_state A SimpleState corresponding to the closest point found on the phase envelope
     */
    static bool is_inside(HelmholtzEOSMixtureBackend &HEOS, parameters iInput1, CoolPropDbl value1, parameters iInput2, CoolPropDbl value2, std::size_t &iclosest, SimpleState &closest_state);
};
    
} /* namespace CoolProp */