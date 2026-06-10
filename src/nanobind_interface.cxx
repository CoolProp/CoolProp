#ifdef NANOBIND

#    include "CoolProp/CoolProp.h"
#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/HumidAirProp.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/numerics/numerics.h"  // ValidNumber: scalar PropsSI/HAPropsSI raise on non-finite (Cython parity)
#    include "CoolProp/detail/state_capi.h"
#    include "Backends/Helmholtz/MixtureParameters.h"

#    include <nanobind/nanobind.h>
#    include <nanobind/stl/string.h>
#    include <nanobind/stl/string_view.h>  // get_parameter_information takes std::string_view
#    include <nanobind/stl/vector.h>
#    include <nanobind/stl/tuple.h>
#    include <nanobind/stl/pair.h>
#    include <nanobind/stl/map.h>
#    include <nanobind/ndarray.h>  // vectorized PropsSI returns a numpy array
#    include <algorithm>
#    include <array>
namespace nb = nanobind;

CoolProp::AbstractState* factory(const std::string& backend, const std::string& fluid_names) {
    return CoolProp::AbstractState::factory(backend, fluid_names);
}

// Split a delimited param string into a list, mirroring the str.split() used by
// the Cython high-level convenience wrappers below (FluidsList, get_aliases).
static std::vector<std::string> _split_str(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == delim) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    // Keep Python str.split() semantics for non-empty input (trailing/empty
    // tokens preserved), but return [] -- not [""] -- for an empty input, so
    // e.g. an empty FluidsList/aliases string yields an empty list.
    if (!cur.empty() || !s.empty()) {
        out.push_back(cur);
    }
    return out;
}

// Scalar-or-sequence -> vector<double> for the vectorized PropsSI/HAPropsSI
// overloads; is_seq reports whether the input was array-like (vs a scalar).
//
// Mirrors the legacy Cython `iterable()` test: a 1-D list/tuple/ndarray is a
// sequence; a 0-D array, a numpy scalar (np.float64/np.int64), a Python int or
// float -- anything with __float__ -- is a SCALAR.  Crucially this keeps the
// canonical int-literal call PropsSI("T","P",101325,...) on the scalar path so
// it returns a float, not a 1-element ndarray (bd CoolProp-r9sq.21).
static std::vector<double> _to_vec(nb::handle o, bool& is_seq) {
    if (nb::hasattr(o, "ndim")) {  // numpy array / scalar: dispatch on ndim
        long ndim = nb::cast<long>(nb::getattr(o, "ndim"));
        if (ndim > 1) {
            // Reject multi-dimensional input rather than silently flattening it.
            PyErr_SetString(PyExc_ValueError, "vectorized PropsSI/HAPropsSI input is not one-dimensional");
            throw nb::python_error();
        }
        if (ndim == 0) {  // 0-D array is a scalar
            is_seq = false;
            return {nb::cast<double>(nb::float_(o))};
        }
        is_seq = true;  // 1-D ndarray
        return nb::cast<std::vector<double>>(o);
    }
    if (nb::isinstance<nb::list>(o) || nb::isinstance<nb::tuple>(o)) {
        is_seq = true;
        return nb::cast<std::vector<double>>(o);
    }
    // Scalar: Python int/float, numpy scalar, or anything with __float__.  Coerce
    // through nb::float_ (PyNumber_Float) so int/np.float64 are accepted as scalars.
    is_seq = false;
    return {nb::cast<double>(nb::float_(o))};
}

// Raise a Python ValueError carrying CoolProp's global error string if x is not
// finite -- the legacy scalar PropsSI/HAPropsSI raised ValueError rather than
// returning inf/nan (bd CoolProp-r9sq.21/.22).
static void _raise_if_invalid(double x) {
    // Qualified: this helper is at file scope, before init_CoolProp's `using namespace CoolProp`.
    // (ValidNumber lives in the global namespace; get_global_param_string in CoolProp.)
    if (!ValidNumber(x)) {
        PyErr_SetString(PyExc_ValueError, CoolProp::get_global_param_string("errstring").c_str());
        throw nb::python_error();
    }
}

// Broadcast a size-1 vector up to length n; lengths other than 1 or n are an error.
static void _broadcast_to(std::vector<double>& v, std::size_t n, const char* which) {
    if (v.size() == n) {
        return;
    }
    if (v.size() == 1) {
        double val = v[0];  // copy first: v[0] aliases the vector that assign() clears
        v.assign(n, val);
        return;
    }
    // Match the legacy Cython wrapper, which raises TypeError on a length mismatch.
    PyErr_SetString(PyExc_TypeError, (std::string("vectorized input ") + which + " has an incompatible length").c_str());
    throw nb::python_error();
}

// ---- State C-ABI capsule (bridge for the frozen Cython `State` compat shim) --
// The opaque handle carries a shared_ptr<AbstractState> directly.
namespace {
using _CAPI_SP = std::shared_ptr<CoolProp::AbstractState>;
thread_local std::string g_capi_error;  // message for the most recent failure; "" == none

// Stash an error message without ever throwing out of the extern "C" frame (the
// assignment allocates, hence the catch).  On OOM the message is left empty and
// the shim falls back to its own generic text.
void capi_set_error(const char* msg) noexcept {
    try {
        g_capi_error = (msg != nullptr) ? msg : "unknown C++ exception";
    } catch (...) {
        g_capi_error.clear();  // clear() is noexcept
    }
}
void* capi_make(const char* backend, const char* fluids) {
    try {
        auto* p = new _CAPI_SP(CoolProp::AbstractState::factory(backend, fluids));
        g_capi_error.clear();
        return p;
    } catch (const std::exception& e) {
        capi_set_error(e.what());
        return nullptr;
    } catch (...) {
        capi_set_error(nullptr);
        return nullptr;
    }
}
void capi_destroy(void* h) {
    try {
        delete static_cast<_CAPI_SP*>(h);
    } catch (...) {  // ~AbstractState is effectively noexcept; defensive against future changes
    }
}
void capi_update(void* h, long input_pair, double v1, double v2) {
    try {
        (*static_cast<_CAPI_SP*>(h))->update(static_cast<CoolProp::input_pairs>(input_pair), v1, v2);
        g_capi_error.clear();
    } catch (const std::exception& e) {
        capi_set_error(e.what());
    } catch (...) {
        capi_set_error(nullptr);
    }
}
double capi_keyed_output(void* h, long key) {
    try {
        double v = (*static_cast<_CAPI_SP*>(h))->keyed_output(static_cast<CoolProp::parameters>(key));
        g_capi_error.clear();
        return v;
    } catch (const std::exception& e) {
        capi_set_error(e.what());
        return NAN;
    } catch (...) {
        capi_set_error(nullptr);
        return NAN;
    }
}
double capi_first_partial_deriv(void* h, long Of, long Wrt, long Constant) {
    try {
        double v = (*static_cast<_CAPI_SP*>(h))
                     ->first_partial_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt),
                                           static_cast<CoolProp::parameters>(Constant));
        g_capi_error.clear();
        return v;
    } catch (const std::exception& e) {
        capi_set_error(e.what());
        return NAN;
    } catch (...) {
        capi_set_error(nullptr);
        return NAN;
    }
}
const char* capi_last_error() {
    return g_capi_error.empty() ? nullptr : g_capi_error.c_str();
}
void capi_set_mole_fractions(void* h, const double* z, long n) {
    // Validate the C-ABI pointers/length before constructing or dereferencing
    // (a null handle/fractions or negative length from any consumer must not crash).
    if (h == nullptr) {
        capi_set_error("set_mole_fractions: null handle");
        return;
    }
    if (n < 0) {
        capi_set_error("set_mole_fractions: negative length");
        return;
    }
    if (n > 0 && z == nullptr) {
        capi_set_error("set_mole_fractions: null fractions pointer");
        return;
    }
    try {
        std::vector<double> fractions(z, z + n);  // binds the vector<double> set_mole_fractions overload
        (*static_cast<_CAPI_SP*>(h))->set_mole_fractions(fractions);
        g_capi_error.clear();
    } catch (const std::exception& e) {
        capi_set_error(e.what());
    } catch (...) {
        capi_set_error(nullptr);
    }
}
void capi_specify_phase(void* h, long phase) {
    if (h == nullptr) {
        capi_set_error("specify_phase: null handle");
        return;
    }
    try {
        (*static_cast<_CAPI_SP*>(h))->specify_phase(static_cast<CoolProp::phases>(phase));
        g_capi_error.clear();
    } catch (const std::exception& e) {
        capi_set_error(e.what());
    } catch (...) {
        capi_set_error(nullptr);
    }
}
const CoolProp_StateCAPI g_state_capi = {
  capi_make, capi_destroy, capi_update, capi_keyed_output, capi_first_partial_deriv, capi_last_error, capi_set_mole_fractions, capi_specify_phase};
}  // namespace

void init_CoolProp(nb::module_& m) {
    using namespace CoolProp;

    nb::class_<SimpleState>(m, "SimpleState")
      .def(nb::init<>())
      .def_rw("T", &SimpleState::T)
      .def_rw("p", &SimpleState::p)
      .def_rw("rhomolar", &SimpleState::rhomolar);

    nb::class_<GuessesStructure>(m, "GuessesStructure")
      .def(nb::init<>())
      .def_rw("T", &GuessesStructure::T)
      .def_rw("p", &GuessesStructure::p)
      .def_rw("rhomolar", &GuessesStructure::rhomolar)
      .def_rw("hmolar", &GuessesStructure::hmolar)
      .def_rw("smolar", &GuessesStructure::smolar)
      .def_rw("rhomolar_liq", &GuessesStructure::rhomolar_liq)
      .def_rw("rhomolar_vap", &GuessesStructure::rhomolar_vap)
      .def_rw("x", &GuessesStructure::x)
      .def_rw("y", &GuessesStructure::y)
      .def("clear", &GuessesStructure::clear);

    nb::class_<CriticalState, SimpleState>(m, "CriticalState").def(nb::init<>()).def_rw("stable", &CriticalState::stable);

    nb::class_<PhaseEnvelopeData>(m, "PhaseEnvelopeData")
      .def(nb::init<>())
      .def_rw("K", &PhaseEnvelopeData::K)
      .def_rw("lnK", &PhaseEnvelopeData::lnK)
      .def_rw("x", &PhaseEnvelopeData::x)
      .def_rw("y", &PhaseEnvelopeData::y)
      .def_rw("T", &PhaseEnvelopeData::T)
      .def_rw("p", &PhaseEnvelopeData::p)
      .def_rw("lnT", &PhaseEnvelopeData::lnT)
      .def_rw("lnp", &PhaseEnvelopeData::lnp)
      .def_rw("rhomolar_liq", &PhaseEnvelopeData::rhomolar_liq)
      .def_rw("rhomolar_vap", &PhaseEnvelopeData::rhomolar_vap)
      .def_rw("lnrhomolar_liq", &PhaseEnvelopeData::lnrhomolar_liq)
      .def_rw("lnrhomolar_vap", &PhaseEnvelopeData::lnrhomolar_vap)
      .def_rw("hmolar_liq", &PhaseEnvelopeData::hmolar_liq)
      .def_rw("hmolar_vap", &PhaseEnvelopeData::hmolar_vap)
      .def_rw("smolar_liq", &PhaseEnvelopeData::smolar_liq)
      .def_rw("smolar_vap", &PhaseEnvelopeData::smolar_vap)
      .def_rw("Q", &PhaseEnvelopeData::Q)
      .def_rw("cpmolar_liq", &PhaseEnvelopeData::cpmolar_liq)
      .def_rw("cpmolar_vap", &PhaseEnvelopeData::cpmolar_vap)
      .def_rw("cvmolar_liq", &PhaseEnvelopeData::cvmolar_liq)
      .def_rw("cvmolar_vap", &PhaseEnvelopeData::cvmolar_vap)
      .def_rw("viscosity_liq", &PhaseEnvelopeData::viscosity_liq)
      .def_rw("viscosity_vap", &PhaseEnvelopeData::viscosity_vap)
      .def_rw("conductivity_liq", &PhaseEnvelopeData::conductivity_liq)
      .def_rw("conductivity_vap", &PhaseEnvelopeData::conductivity_vap)
      .def_rw("speed_sound_vap", &PhaseEnvelopeData::speed_sound_vap);

    // See http://stackoverflow.com/a/148610 and http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
    nb::enum_<configuration_keys>(m, "configuration_keys", nb::is_arithmetic())
#    define X(Enum, String, Default, Desc) .value(String, configuration_keys::Enum)
      CONFIGURATION_KEYS_ENUM
#    undef X
        .export_values();

    nb::enum_<parameters>(m, "parameters", nb::is_arithmetic())
      .value("igas_constant", parameters::igas_constant)
      .value("imolar_mass", parameters::imolar_mass)
      .value("iacentric_factor", parameters::iacentric_factor)
      .value("irhomolar_reducing", parameters::irhomolar_reducing)
      .value("irhomolar_critical", parameters::irhomolar_critical)
      .value("iT_reducing", parameters::iT_reducing)
      .value("iT_critical", parameters::iT_critical)
      .value("irhomass_reducing", parameters::irhomass_reducing)
      .value("irhomass_critical", parameters::irhomass_critical)
      .value("iP_critical", parameters::iP_critical)
      .value("iP_reducing", parameters::iP_reducing)
      .value("iT_triple", parameters::iT_triple)
      .value("iP_triple", parameters::iP_triple)
      .value("iT_min", parameters::iT_min)
      .value("iT_max", parameters::iT_max)
      .value("iP_max", parameters::iP_max)
      .value("iP_min", parameters::iP_min)
      .value("idipole_moment", parameters::idipole_moment)
      .value("iT", parameters::iT)
      .value("iP", parameters::iP)
      .value("iQ", parameters::iQ)
      .value("iTau", parameters::iTau)
      .value("iDelta", parameters::iDelta)
      .value("iDmolar", parameters::iDmolar)
      .value("iHmolar", parameters::iHmolar)
      .value("iSmolar", parameters::iSmolar)
      .value("iCpmolar", parameters::iCpmolar)
      .value("iCp0molar", parameters::iCp0molar)
      .value("iCvmolar", parameters::iCvmolar)
      .value("iUmolar", parameters::iUmolar)
      .value("iGmolar", parameters::iGmolar)
      .value("iHelmholtzmolar", parameters::iHelmholtzmolar)
      .value("iSmolar_residual", parameters::iSmolar_residual)
      .value("iHmolar_residual", parameters::iHmolar_residual)
      .value("iGmolar_residual", parameters::iGmolar_residual)
      .value("iDmass", parameters::iDmass)
      .value("iHmass", parameters::iHmass)
      .value("iSmass", parameters::iSmass)
      .value("iCpmass", parameters::iCpmass)
      .value("iCp0mass", parameters::iCp0mass)
      .value("iCvmass", parameters::iCvmass)
      .value("iUmass", parameters::iUmass)
      .value("iGmass", parameters::iGmass)
      .value("iHelmholtzmass", parameters::iHelmholtzmass)
      .value("iviscosity", parameters::iviscosity)
      .value("iconductivity", parameters::iconductivity)
      .value("isurface_tension", parameters::isurface_tension)
      .value("iPrandtl", parameters::iPrandtl)
      .value("ispeed_sound", parameters::ispeed_sound)
      .value("iisothermal_compressibility", parameters::iisothermal_compressibility)
      .value("iisobaric_expansion_coefficient", parameters::iisobaric_expansion_coefficient)
      .value("ifundamental_derivative_of_gas_dynamics", parameters::ifundamental_derivative_of_gas_dynamics)
      .value("ialphar", parameters::ialphar)
      .value("idalphar_ddelta_consttau", parameters::idalphar_ddelta_consttau)
      .value("idalpha0_dtau_constdelta", parameters::idalpha0_dtau_constdelta)
      .value("iBvirial", parameters::iBvirial)
      .value("iCvirial", parameters::iCvirial)
      .value("idBvirial_dT", parameters::idBvirial_dT)
      .value("idCvirial_dT", parameters::idCvirial_dT)
      .value("iZ", parameters::iZ)
      .value("iPIP", parameters::iPIP)
      .value("ifraction_min", parameters::ifraction_min)
      .value("ifraction_max", parameters::ifraction_max)
      .value("iT_freeze", parameters::iT_freeze)
      .value("iGWP20", parameters::iGWP20)
      .value("iGWP100", parameters::iGWP100)
      .value("iGWP500", parameters::iGWP500)
      .value("iFH", parameters::iFH)
      .value("iHH", parameters::iHH)
      .value("iPH", parameters::iPH)
      .value("iODP", parameters::iODP)
      .value("iPhase", parameters::iPhase)
      // bd CoolProp-r9sq.23: values the hand-written enum had drifted from
      // DataStructures.h (mass-basis quality, ideal-gas decompositions, the
      // isentropic-expansion + alpha0 derivatives).
      .value("INVALID_PARAMETER", parameters::INVALID_PARAMETER)
      .value("iQmass", parameters::iQmass)
      .value("iHmolar_idealgas", parameters::iHmolar_idealgas)
      .value("iSmolar_idealgas", parameters::iSmolar_idealgas)
      .value("iUmolar_idealgas", parameters::iUmolar_idealgas)
      .value("iHmass_idealgas", parameters::iHmass_idealgas)
      .value("iSmass_idealgas", parameters::iSmass_idealgas)
      .value("iUmass_idealgas", parameters::iUmass_idealgas)
      .value("iisentropic_expansion_coefficient", parameters::iisentropic_expansion_coefficient)
      .value("idalphar_dtau_constdelta", parameters::idalphar_dtau_constdelta)
      .value("ialpha0", parameters::ialpha0)
      .value("idalpha0_ddelta_consttau", parameters::idalpha0_ddelta_consttau)
      .value("id2alpha0_ddelta2_consttau", parameters::id2alpha0_ddelta2_consttau)
      .value("id3alpha0_ddelta3_consttau", parameters::id3alpha0_ddelta3_consttau)
      .value("iundefined_parameter", parameters::iundefined_parameter)
      .export_values();

    nb::enum_<input_pairs>(m, "input_pairs", nb::is_arithmetic())
      .value("QT_INPUTS", input_pairs::QT_INPUTS)
      .value("PQ_INPUTS", input_pairs::PQ_INPUTS)
      .value("QSmolar_INPUTS", input_pairs::QSmolar_INPUTS)
      .value("QSmass_INPUTS", input_pairs::QSmass_INPUTS)
      .value("HmolarQ_INPUTS", input_pairs::HmolarQ_INPUTS)
      .value("HmassQ_INPUTS", input_pairs::HmassQ_INPUTS)
      .value("DmolarQ_INPUTS", input_pairs::DmolarQ_INPUTS)
      .value("DmassQ_INPUTS", input_pairs::DmassQ_INPUTS)
      .value("PT_INPUTS", input_pairs::PT_INPUTS)
      .value("DmassT_INPUTS", input_pairs::DmassT_INPUTS)
      .value("DmolarT_INPUTS", input_pairs::DmolarT_INPUTS)
      .value("HmolarT_INPUTS", input_pairs::HmolarT_INPUTS)
      .value("HmassT_INPUTS", input_pairs::HmassT_INPUTS)
      .value("SmolarT_INPUTS", input_pairs::SmolarT_INPUTS)
      .value("SmassT_INPUTS", input_pairs::SmassT_INPUTS)
      .value("TUmolar_INPUTS", input_pairs::TUmolar_INPUTS)
      .value("TUmass_INPUTS", input_pairs::TUmass_INPUTS)
      .value("DmassP_INPUTS", input_pairs::DmassP_INPUTS)
      .value("DmolarP_INPUTS", input_pairs::DmolarP_INPUTS)
      .value("HmassP_INPUTS", input_pairs::HmassP_INPUTS)
      .value("HmolarP_INPUTS", input_pairs::HmolarP_INPUTS)
      .value("PSmass_INPUTS", input_pairs::PSmass_INPUTS)
      .value("PSmolar_INPUTS", input_pairs::PSmolar_INPUTS)
      .value("PUmass_INPUTS", input_pairs::PUmass_INPUTS)
      .value("PUmolar_INPUTS", input_pairs::PUmolar_INPUTS)
      .value("HmassSmass_INPUTS", input_pairs::HmassSmass_INPUTS)
      .value("HmolarSmolar_INPUTS", input_pairs::HmolarSmolar_INPUTS)
      .value("SmassUmass_INPUTS", input_pairs::SmassUmass_INPUTS)
      .value("SmolarUmolar_INPUTS", input_pairs::SmolarUmolar_INPUTS)
      .value("DmassHmass_INPUTS", input_pairs::DmassHmass_INPUTS)
      .value("DmolarHmolar_INPUTS", input_pairs::DmolarHmolar_INPUTS)
      .value("DmassSmass_INPUTS", input_pairs::DmassSmass_INPUTS)
      .value("DmolarSmolar_INPUTS", input_pairs::DmolarSmolar_INPUTS)
      .value("DmassUmass_INPUTS", input_pairs::DmassUmass_INPUTS)
      .value("DmolarUmolar_INPUTS", input_pairs::DmolarUmolar_INPUTS)
      // bd CoolProp-r9sq.23: the mass-basis quality pairs (needed for
      // AbstractState.update(QmassT_INPUTS, ...)) + INVALID were missing.
      .value("INPUT_PAIR_INVALID", input_pairs::INPUT_PAIR_INVALID)
      .value("QmassT_INPUTS", input_pairs::QmassT_INPUTS)
      .value("PQmass_INPUTS", input_pairs::PQmass_INPUTS)
      .value("QmassSmolar_INPUTS", input_pairs::QmassSmolar_INPUTS)
      .value("QmassSmass_INPUTS", input_pairs::QmassSmass_INPUTS)
      .value("HmolarQmass_INPUTS", input_pairs::HmolarQmass_INPUTS)
      .value("HmassQmass_INPUTS", input_pairs::HmassQmass_INPUTS)
      .value("DmolarQmass_INPUTS", input_pairs::DmolarQmass_INPUTS)
      .value("DmassQmass_INPUTS", input_pairs::DmassQmass_INPUTS)
      .export_values();

    nb::enum_<phases>(m, "phases", nb::is_arithmetic())
      .value("iphase_liquid", phases::iphase_liquid)
      .value("iphase_supercritical", phases::iphase_supercritical)
      .value("iphase_supercritical_gas", phases::iphase_supercritical_gas)
      .value("iphase_supercritical_liquid", phases::iphase_supercritical_liquid)
      .value("iphase_critical_point", phases::iphase_critical_point)
      .value("iphase_gas", phases::iphase_gas)
      .value("iphase_twophase", phases::iphase_twophase)
      .value("iphase_unknown", phases::iphase_unknown)
      .value("iphase_not_imposed", phases::iphase_not_imposed)
      .export_values();

    // bd CoolProp-r9sq.23: enums the legacy interface exposed but the nanobind
    // interface did not bind at all.
    nb::enum_<fluid_types>(m, "fluid_types", nb::is_arithmetic())
      .value("FLUID_TYPE_PURE", fluid_types::FLUID_TYPE_PURE)
      .value("FLUID_TYPE_PSEUDOPURE", fluid_types::FLUID_TYPE_PSEUDOPURE)
      .value("FLUID_TYPE_REFPROP", fluid_types::FLUID_TYPE_REFPROP)
      .value("FLUID_TYPE_INCOMPRESSIBLE_LIQUID", fluid_types::FLUID_TYPE_INCOMPRESSIBLE_LIQUID)
      .value("FLUID_TYPE_INCOMPRESSIBLE_SOLUTION", fluid_types::FLUID_TYPE_INCOMPRESSIBLE_SOLUTION)
      .value("FLUID_TYPE_UNDEFINED", fluid_types::FLUID_TYPE_UNDEFINED)
      .export_values();

    nb::enum_<fast_evaluate_status>(m, "fast_evaluate_status", nb::is_arithmetic())
      .value("fast_evaluate_ok", fast_evaluate_status::fast_evaluate_ok)
      .value("fast_evaluate_out_of_range", fast_evaluate_status::fast_evaluate_out_of_range)
      .value("fast_evaluate_two_phase_disallowed", fast_evaluate_status::fast_evaluate_two_phase_disallowed)
      .value("fast_evaluate_unsupported_input", fast_evaluate_status::fast_evaluate_unsupported_input)
      .value("fast_evaluate_unsupported_output", fast_evaluate_status::fast_evaluate_unsupported_output)
      .value("fast_evaluate_internal_error", fast_evaluate_status::fast_evaluate_internal_error)
      .export_values();

    // bd CoolProp-r9sq.24: register SpinodalData so get_spinodal_data() (below)
    // can convert its return value instead of raising a cast error.
    nb::class_<SpinodalData>(m, "SpinodalData")
      .def(nb::init<>())
      .def_ro("tau", &SpinodalData::tau)
      .def_ro("delta", &SpinodalData::delta)
      .def_ro("M1", &SpinodalData::M1);

    // bd CoolProp-r9sq.18/.24: the bound state structs use the C++ names
    // (CriticalState, GuessesStructure, PhaseEnvelopeData, SpinodalData); retain
    // the legacy Cython Py*-prefixed names as backwards-compatibility aliases so
    // downstream code (e.g. CoolProp.Plots.Common imports + constructs
    // PyCriticalState / PyGuessesStructure) keeps working against the v8 core.
    m.attr("PyCriticalState") = m.attr("CriticalState");
    m.attr("PyGuessesStructure") = m.attr("GuessesStructure");
    m.attr("PyPhaseEnvelopeData") = m.attr("PhaseEnvelopeData");
    m.attr("PySpinodalData") = m.attr("SpinodalData");

    // bd CoolProp-r9sq.28: register AbstractState as a real TYPE with a factory
    // constructor (nb::new_), not a module-level factory FUNCTION returning a
    // private _AbstractState.  So `AbstractState("HEOS","Water")` still builds a
    // state AND `isinstance(x, AbstractState)` works (it previously raised
    // "arg 2 must be a type", forcing callers like Plots/Common.py to reach into
    // the private class).  nb::new_ binds __new__ to `factory` + a no-op __init__.
    nb::class_<AbstractState>(m, "AbstractState")
      .def(nb::new_(&factory))
      .def("set_T", &AbstractState::set_T)
      .def("backend_name", &AbstractState::backend_name)
      .def("using_mole_fractions", &AbstractState::using_mole_fractions)
      .def("using_mass_fractions", &AbstractState::using_mass_fractions)
      .def("using_volu_fractions", &AbstractState::using_volu_fractions)
      .def("set_mole_fractions", &AbstractState::set_mole_fractions)
      .def("set_mass_fractions", &AbstractState::set_mass_fractions)
      .def("set_volu_fractions", &AbstractState::set_volu_fractions)
      .def("mole_fractions_liquid", &AbstractState::mole_fractions_liquid)
      .def("mole_fractions_liquid_double", &AbstractState::mole_fractions_liquid_double)
      .def("mole_fractions_vapor", &AbstractState::mole_fractions_vapor)
      .def("mole_fractions_vapor_double", &AbstractState::mole_fractions_vapor_double)
      .def("get_mole_fractions", &AbstractState::get_mole_fractions)
      .def("get_mass_fractions", &AbstractState::get_mass_fractions)
      .def("update", &AbstractState::update)
      .def("update_with_guesses", &AbstractState::update_with_guesses)
      .def("available_in_high_level", &AbstractState::available_in_high_level)
      .def("build_options_json", &AbstractState::build_options_json)  // bd CoolProp-r9sq.25
      .def("fluid_param_string", &AbstractState::fluid_param_string)
      .def("fluid_names", &AbstractState::fluid_names)
      .def("set_binary_interaction_double", (void(AbstractState::*)(const std::string&, const std::string&, const std::string&, const double))
                                              & AbstractState::set_binary_interaction_double)
      .def("set_binary_interaction_double", (void(AbstractState::*)(const std::size_t, const std::size_t, const std::string&, const double))
                                              & AbstractState::set_binary_interaction_double)
      .def("set_binary_interaction_string", (void(AbstractState::*)(const std::string&, const std::string&, const std::string&, const std::string&))
                                              & AbstractState::set_binary_interaction_string)
      .def("set_binary_interaction_string", (void(AbstractState::*)(const std::size_t, const std::size_t, const std::string&, const std::string&))
                                              & AbstractState::set_binary_interaction_string)
      .def("get_binary_interaction_double",
           (double(AbstractState::*)(const std::string&, const std::string&, const std::string&)) & AbstractState::get_binary_interaction_double)
      .def("get_binary_interaction_double",
           (double(AbstractState::*)(const std::size_t, const std::size_t, const std::string&)) & AbstractState::get_binary_interaction_double)
      .def("get_binary_interaction_string", &AbstractState::get_binary_interaction_string)
      .def("apply_simple_mixing_rule", &AbstractState::apply_simple_mixing_rule)
      .def("set_fluid_parameter_double", &AbstractState::set_fluid_parameter_double)
      .def("clear", &AbstractState::clear)
      .def("get_reducing_state", &AbstractState::get_reducing_state)
      .def("get_state", &AbstractState::get_state)
      .def("Tmin", &AbstractState::Tmin)
      .def("Tmax", &AbstractState::Tmax)
      .def("pmax", &AbstractState::pmax)
      .def("Ttriple", &AbstractState::Ttriple)
      .def("phase", &AbstractState::phase)
      .def("specify_phase", &AbstractState::specify_phase)
      .def("unspecify_phase", &AbstractState::unspecify_phase)
      .def("T_critical", &AbstractState::T_critical)
      .def("p_critical", &AbstractState::p_critical)
      .def("rhomolar_critical", &AbstractState::rhomolar_critical)
      .def("rhomass_critical", &AbstractState::rhomass_critical)
      .def("all_critical_points", &AbstractState::all_critical_points)
      .def("build_spinodal", &AbstractState::build_spinodal)
      .def("get_spinodal_data", &AbstractState::get_spinodal_data)
      .def("criticality_contour_values",
           [](AbstractState& AS) {
               double L, M;
               AS.criticality_contour_values(L, M);
               return nb::make_tuple(L, M);
           })
      .def("tangent_plane_distance", &AbstractState::tangent_plane_distance)
      .def("T_reducing", &AbstractState::T_reducing)
      .def("rhomolar_reducing", &AbstractState::rhomolar_reducing)
      .def("rhomass_reducing", &AbstractState::rhomass_reducing)
      .def("p_triple", &AbstractState::p_triple)
      .def("name", &AbstractState::name)
      .def("dipole_moment", &AbstractState::dipole_moment)
      .def("keyed_output", &AbstractState::keyed_output)
      .def("trivial_keyed_output", &AbstractState::trivial_keyed_output)
      // bd CoolProp-r9sq.25: zero-allocation vectorized batch path (tabular / IF97
      // backends).  Caller-allocated, shape-validated contiguous buffers, mirroring
      // the legacy Cython fast_evaluate; out is (N_inputs, N_outputs), filled in place.
      .def(
        "fast_evaluate",
        [](AbstractState& AS, input_pairs input_pair, nb::ndarray<const double, nb::ndim<1>, nb::c_contig> val1,
           nb::ndarray<const double, nb::ndim<1>, nb::c_contig> val2, nb::ndarray<const int, nb::ndim<1>, nb::c_contig> outputs,
           nb::ndarray<double, nb::ndim<2>, nb::c_contig> out, nb::ndarray<int, nb::ndim<1>, nb::c_contig> status, phases imposed_phase) {
            std::size_t N = val1.shape(0);
            std::size_t M = outputs.shape(0);
            if (val2.shape(0) != N) {
                throw nb::value_error("val1 and val2 must have the same length");
            }
            if (out.shape(0) != N || out.shape(1) != M) {
                throw nb::value_error("out must have shape (N_inputs, N_outputs)");
            }
            if (status.shape(0) != N) {
                throw nb::value_error("status must have length N_inputs");
            }
            if (N == 0) {
                return;
            }
            if (M == 0) {  // nothing to compute per point; C++ contract: status=ok, no output writes
                for (std::size_t k = 0; k < N; ++k) {
                    status.data()[k] = 0;
                }
                return;
            }
            AS.fast_evaluate(input_pair, val1.data(), val2.data(), N, reinterpret_cast<const parameters*>(outputs.data()), M, out.data(), N * M,
                             status.data(), N, imposed_phase);
        },
        nb::arg("input_pair"), nb::arg("val1"), nb::arg("val2"), nb::arg("outputs"), nb::arg("out"), nb::arg("status"),
        nb::arg("imposed_phase") = phases::iphase_not_imposed)
      .def("saturated_liquid_keyed_output", &AbstractState::saturated_liquid_keyed_output)
      .def("saturated_vapor_keyed_output", &AbstractState::saturated_vapor_keyed_output)
      .def("T", &AbstractState::T)
      .def("rhomolar", &AbstractState::rhomolar)
      .def("rhomass", &AbstractState::rhomass)
      .def("p", &AbstractState::p)
      .def("Q", &AbstractState::Q)
      .def("Qmass", &AbstractState::Qmass)
      .def("tau", &AbstractState::tau)
      .def("delta", &AbstractState::delta)
      .def("molar_mass", &AbstractState::molar_mass)
      .def("acentric_factor", &AbstractState::acentric_factor)
      .def("gas_constant", &AbstractState::gas_constant)
      .def("Bvirial", &AbstractState::Bvirial)
      .def("dBvirial_dT", &AbstractState::dBvirial_dT)
      .def("Cvirial", &AbstractState::Cvirial)
      .def("dCvirial_dT", &AbstractState::dCvirial_dT)
      .def("compressibility_factor", &AbstractState::compressibility_factor)
      .def("hmolar", &AbstractState::hmolar)
      .def("hmass", &AbstractState::hmass)
      .def("hmolar_excess", &AbstractState::hmolar_excess)
      .def("hmass_excess", &AbstractState::hmass_excess)
      .def("smolar", &AbstractState::smolar)
      .def("smass", &AbstractState::smass)
      .def("smolar_excess", &AbstractState::smolar_excess)
      .def("smass_excess", &AbstractState::smass_excess)
      .def("umolar", &AbstractState::umolar)
      .def("umass", &AbstractState::umass)
      .def("umolar_excess", &AbstractState::umolar_excess)
      .def("umass_excess", &AbstractState::umass_excess)
      .def("cpmolar", &AbstractState::cpmolar)
      .def("cpmass", &AbstractState::cpmass)
      .def("cp0molar", &AbstractState::cp0molar)
      .def("cp0mass", &AbstractState::cp0mass)
      .def("cvmolar", &AbstractState::cvmolar)
      .def("cvmass", &AbstractState::cvmass)
      .def("gibbsmolar", &AbstractState::gibbsmolar)
      .def("gibbsmass", &AbstractState::gibbsmass)
      .def("gibbsmolar_excess", &AbstractState::gibbsmolar_excess)
      .def("gibbsmass_excess", &AbstractState::gibbsmass_excess)
      .def("helmholtzmolar", &AbstractState::helmholtzmolar)
      .def("helmholtzmass", &AbstractState::helmholtzmass)
      .def("helmholtzmolar_excess", &AbstractState::helmholtzmolar_excess)
      .def("helmholtzmass_excess", &AbstractState::helmholtzmass_excess)
      .def("volumemolar_excess", &AbstractState::volumemolar_excess)
      .def("volumemass_excess", &AbstractState::volumemass_excess)
      // q2sh: Cython-parity additions -- idealgas/residual decompositions
      .def("hmolar_idealgas", &AbstractState::hmolar_idealgas)
      .def("hmass_idealgas", &AbstractState::hmass_idealgas)
      .def("hmolar_residual", &AbstractState::hmolar_residual)
      .def("smolar_idealgas", &AbstractState::smolar_idealgas)
      .def("smass_idealgas", &AbstractState::smass_idealgas)
      .def("smolar_residual", &AbstractState::smolar_residual)
      .def("umolar_idealgas", &AbstractState::umolar_idealgas)
      .def("umass_idealgas", &AbstractState::umass_idealgas)
      .def("gibbsmolar_residual", &AbstractState::gibbsmolar_residual)
      .def("neff", &AbstractState::neff)
      // q2sh: fluid-constant / cubic-alpha / superancillary accessors
      .def("get_fluid_constant", &AbstractState::get_fluid_constant)
      .def("get_fluid_parameter_double", &AbstractState::get_fluid_parameter_double)
      .def("set_cubic_alpha_C", &AbstractState::set_cubic_alpha_C)
      .def("update_QT_pure_superanc", &AbstractState::update_QT_pure_superanc)
      .def("speed_sound", &AbstractState::speed_sound)
      .def("isothermal_compressibility", &AbstractState::isothermal_compressibility)
      .def("isobaric_expansion_coefficient", &AbstractState::isobaric_expansion_coefficient)
      .def("fugacity_coefficient", &AbstractState::fugacity_coefficient)
      .def("fugacity", &AbstractState::fugacity)
      .def("chemical_potential", &AbstractState::chemical_potential)
      .def("fundamental_derivative_of_gas_dynamics", &AbstractState::fundamental_derivative_of_gas_dynamics)
      .def("PIP", &AbstractState::PIP)
      // bd CoolProp-r9sq.24: out-reference params -> return a (T, rho) tuple.
      .def("true_critical_point",
           [](AbstractState& AS) {
               double T = 0, rho = 0;
               AS.true_critical_point(T, rho);
               return nb::make_tuple(T, rho);
           })
      .def("ideal_curve",
           [](AbstractState& AS, const std::string& name) {
               std::vector<double> T, p;
               AS.ideal_curve(name, T, p);
               return nb::make_tuple(T, p);
           })
      .def("first_partial_deriv", &AbstractState::first_partial_deriv)
      .def("second_partial_deriv", &AbstractState::second_partial_deriv)
      .def("first_saturation_deriv", &AbstractState::first_saturation_deriv)
      .def("second_saturation_deriv", &AbstractState::second_saturation_deriv)
      .def("first_two_phase_deriv", &AbstractState::first_two_phase_deriv)
      .def("second_two_phase_deriv", &AbstractState::second_two_phase_deriv)
      .def("first_two_phase_deriv_splined", &AbstractState::first_two_phase_deriv_splined)
      .def("build_phase_envelope", &AbstractState::build_phase_envelope)
      .def("get_phase_envelope_data", &AbstractState::get_phase_envelope_data)
      .def("has_melting_line", &AbstractState::has_melting_line)
      .def("melting_line", &AbstractState::melting_line)
      .def("saturation_ancillary", &AbstractState::saturation_ancillary)
      .def("viscosity", &AbstractState::viscosity)
      // bd CoolProp-r9sq.24: out-reference params -> return the legacy dict.
      .def("viscosity_contributions",
           [](AbstractState& AS) {
               CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0;
               AS.viscosity_contributions(dilute, initial_density, residual, critical);
               nb::dict d;
               d["dilute"] = static_cast<double>(dilute);
               d["initial_density"] = static_cast<double>(initial_density);
               d["residual"] = static_cast<double>(residual);
               d["critical"] = static_cast<double>(critical);
               return d;
           })
      .def("conductivity", &AbstractState::conductivity)
      .def("conductivity_contributions",
           [](AbstractState& AS) {
               CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0;
               AS.conductivity_contributions(dilute, initial_density, residual, critical);
               nb::dict d;
               d["dilute"] = static_cast<double>(dilute);
               d["initial_density"] = static_cast<double>(initial_density);
               d["residual"] = static_cast<double>(residual);
               d["critical"] = static_cast<double>(critical);
               return d;
           })
      .def("surface_tension", &AbstractState::surface_tension)
      .def("Prandtl", &AbstractState::Prandtl)
      // bd CoolProp-r9sq.24: T/rho are in/out -> return dict(T, rhomolar) like legacy.
      .def("conformal_state",
           [](AbstractState& AS, const std::string& reference_fluid, double T, double rho) {
               CoolPropDbl T0 = T, rho0 = rho;
               AS.conformal_state(reference_fluid, T0, rho0);
               nb::dict d;
               d["T"] = static_cast<double>(T0);
               d["rhomolar"] = static_cast<double>(rho0);
               return d;
           })
      .def("change_EOS", &AbstractState::change_EOS)
      .def("alpha0", &AbstractState::alpha0)
      .def("dalpha0_dDelta", &AbstractState::dalpha0_dDelta)
      .def("dalpha0_dTau", &AbstractState::dalpha0_dTau)
      .def("d2alpha0_dDelta2", &AbstractState::d2alpha0_dDelta2)
      .def("d2alpha0_dDelta_dTau", &AbstractState::d2alpha0_dDelta_dTau)
      .def("d2alpha0_dTau2", &AbstractState::d2alpha0_dTau2)
      .def("d3alpha0_dTau3", &AbstractState::d3alpha0_dTau3)
      .def("d3alpha0_dDelta_dTau2", &AbstractState::d3alpha0_dDelta_dTau2)
      .def("d3alpha0_dDelta2_dTau", &AbstractState::d3alpha0_dDelta2_dTau)
      .def("d3alpha0_dDelta3", &AbstractState::d3alpha0_dDelta3)
      .def("alphar", &AbstractState::alphar)
      .def("dalphar_dDelta", &AbstractState::dalphar_dDelta)
      .def("dalphar_dTau", &AbstractState::dalphar_dTau)
      .def("d2alphar_dDelta2", &AbstractState::d2alphar_dDelta2)
      .def("d2alphar_dDelta_dTau", &AbstractState::d2alphar_dDelta_dTau)
      .def("d2alphar_dTau2", &AbstractState::d2alphar_dTau2)
      .def("d3alphar_dDelta3", &AbstractState::d3alphar_dDelta3)
      .def("d3alphar_dDelta2_dTau", &AbstractState::d3alphar_dDelta2_dTau)
      .def("d3alphar_dDelta_dTau2", &AbstractState::d3alphar_dDelta_dTau2)
      .def("d3alphar_dTau3", &AbstractState::d3alphar_dTau3)
      .def("d4alphar_dDelta4", &AbstractState::d4alphar_dDelta4)
      .def("d4alphar_dDelta3_dTau", &AbstractState::d4alphar_dDelta3_dTau)
      .def("d4alphar_dDelta2_dTau2", &AbstractState::d4alphar_dDelta2_dTau2)
      .def("d4alphar_dDelta_dTau3", &AbstractState::d4alphar_dDelta_dTau3)
      .def("d4alphar_dTau4", &AbstractState::d4alphar_dTau4);

    // NOTE: `AbstractState(backend, fluids)` is the class constructor registered
    // above via nb::new_(&factory) -- no separate module-level factory function
    // (bd CoolProp-r9sq.28), so isinstance(x, AbstractState) works.

    m.def("get_config_as_json_string", &get_config_as_json_string);
    m.def("set_config_as_json_string", &set_config_as_json_string);
    m.def("config_key_description", (std::string(*)(configuration_keys)) & config_key_description);
    m.def("config_key_description", (std::string(*)(const std::string&)) & config_key_description);
    m.def("set_config_string", &set_config_string);
    m.def("set_config_double", &set_config_double);
    m.def("set_departure_functions", &set_departure_functions);
    m.def("set_config_bool", &set_config_bool);
    m.def("get_config_string", &get_config_string);
    m.def("get_config_double", &get_config_double);
    m.def("get_config_bool", &get_config_bool);
    m.def("get_config_int", &get_config_int);
    m.def("set_config_int", &set_config_int);
    m.def("get_parameter_information", &get_parameter_information);
    m.def("get_parameter_index", &get_parameter_index);
    m.def("get_phase_index", &get_phase_index);
    m.def("is_trivial_parameter", &is_trivial_parameter);
    // generate_update_pair has two out-reference params and returns the input_pair;
    // wrap it to return (input_pair, out1, out2) like the legacy Cython wrapper
    // (binding &generate_update_pair<double> directly exposes an unusable signature).
    m.def("generate_update_pair", [](parameters key1, double value1, parameters key2, double value2) {
        double out1 = 0.0, out2 = 0.0;
        input_pairs pair = generate_update_pair<double>(key1, value1, key2, value2, out1, out2);
        return nb::make_tuple(pair, out1, out2);
    });
    m.def("Props1SI", &Props1SI);
    // Scalar PropsSI (all-float args): raise ValueError on a non-finite result
    // rather than silently returning inf/nan, matching the legacy Cython wrapper
    // (bd CoolProp-r9sq.21).
    m.def("PropsSI", [](const std::string& Output, const std::string& Name1, double Prop1, const std::string& Name2, double Prop2,
                        const std::string& FluidName) {
        double val = PropsSI(Output, Name1, Prop1, Name2, Prop2, FluidName);
        _raise_if_invalid(val);
        return val;
    });
    // Vectorized PropsSI: list/ndarray Prop1/Prop2 (with scalar broadcast) dispatch
    // to the C++ PropsSImulti path.  Returns a numpy array for array inputs, but a
    // plain float when every input was scalar -- so the canonical int-literal call
    // PropsSI("T","P",101325,"Q",0,"Water") (whose int args do NOT match the all-double
    // scalar overload above and so bind here) returns a float, not a 1-element ndarray
    // (bd CoolProp-r9sq.21).
    m.def("PropsSI",
          [](const std::string& Output, const std::string& Name1, nb::object Prop1, const std::string& Name2, nb::object Prop2,
             const std::string& FluidName) -> nb::object {
              bool s1 = false, s2 = false;
              std::vector<double> v1 = _to_vec(Prop1, s1);
              std::vector<double> v2 = _to_vec(Prop2, s2);
              bool any_seq = s1 || s2;
              // Empty state inputs -> empty array (mirrors the legacy #2417 short-circuit).
              if ((s1 && v1.empty()) || (s2 && v2.empty())) {
                  auto* empty = new double[1];
                  nb::capsule owner(empty, [](void* p) noexcept { delete[] static_cast<double*>(p); });
                  return nb::cast(nb::ndarray<nb::numpy, double>(empty, {static_cast<std::size_t>(0)}, owner));
              }
              std::size_t n = std::max(v1.size(), v2.size());
              _broadcast_to(v1, n, "Prop1");
              _broadcast_to(v2, n, "Prop2");
              std::string backend, fluid;
              extract_backend(FluidName, backend, fluid);
              std::vector<double> fractions{1.0};
              std::string delimited = extract_fractions(fluid, fractions);
              std::vector<std::string> fluids = _split_str(delimited, '&');
              std::vector<std::vector<double>> out = PropsSImulti({Output}, Name1, v1, Name2, v2, backend, fluids, fractions);
              if (out.empty() || out[0].empty()) {
                  // Legacy raised ValueError (not RuntimeError) on a failed evaluation.
                  PyErr_SetString(PyExc_ValueError, get_global_param_string("errstring").c_str());
                  throw nb::python_error();
              }
              // All-scalar inputs -> a scalar float, guarded for finiteness like the
              // all-float scalar overload above.
              if (!any_seq) {
                  _raise_if_invalid(out[0][0]);
                  return nb::float_(out[0][0]);
              }
              // PropsSImulti returns a matrix indexed [input][output]; a single Output
              // means one column, so take out[i][0] across the inputs.
              std::size_t n_out = out.size();
              auto* data = new double[n_out];
              for (std::size_t i = 0; i < n_out; ++i) {
                  data[i] = out[i][0];
              }
              nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
              return nb::cast(nb::ndarray<nb::numpy, double>(data, {n_out}, owner));
          });
    // Legacy compatibility: the Cython PropsSI also accepted the 2-arg trivial
    // form PropsSI("Tcrit", "Water") (order-lenient).  Restore it as an overload
    // dispatching to Props1SI; raise ValueError on a non-finite result (parity).
    m.def("PropsSI", [](const std::string& Output, const std::string& FluidName) {
        double val = Props1SI(Output, FluidName);
        _raise_if_invalid(val);
        return val;
    });
    m.def("PhaseSI", &PhaseSI);
    m.def("PropsSImulti", &PropsSImulti);
    m.def("get_global_param_string", &get_global_param_string);
    m.def("get_debug_level", &get_debug_level);
    m.def("set_debug_level", &set_debug_level);
    m.def("get_fluid_param_string", &get_fluid_param_string);
    // bd CoolProp-r9sq.24: these have C++ out-reference params; bind the legacy
    // one-string-in, tuple-out forms (extract_backend("REFPROP::Water") -> (backend,
    // fluid); extract_fractions("R32&R125[0.7]") -> (fluids_list, fractions)).
    m.def("extract_backend", [](const std::string& in_str) {
        std::string backend, fluid;
        extract_backend(in_str, backend, fluid);
        return nb::make_tuple(backend, fluid);
    });
    m.def("extract_fractions", [](const std::string& in_str) {
        std::vector<double> fractions;
        std::string delimited = extract_fractions(in_str, fractions);
        return nb::make_tuple(_split_str(delimited, '&'), fractions);
    });
    m.def("set_reference_stateS", &set_reference_stateS);
    m.def("set_reference_stateD", &set_reference_stateD);
    m.def("saturation_ancillary", &saturation_ancillary);
    m.def("add_fluids_as_JSON", &add_fluids_as_JSON);
    // Scalar HAPropsSI (all-float args): raise ValueError on a non-finite result,
    // matching the legacy Cython wrapper (bd CoolProp-r9sq.22).
    m.def("HAPropsSI",
          [](const std::string& Output, const std::string& N1, double V1, const std::string& N2, double V2, const std::string& N3, double V3) {
              double val = HumidAir::HAPropsSI(Output, N1, V1, N2, V2, N3, V3);
              _raise_if_invalid(val);
              return val;
          });
    // Vectorized HAPropsSI: array inputs loop over the scalar C++ HAPropsSI.  Returns
    // a float when every input was scalar (so HAPropsSI('H','T',298.15,'P',101325,'R',0.5)
    // with an int pressure returns a float, not a 1-element list), a numpy array when any
    // input was an ndarray (preserving the array type), else a list -- matching the legacy
    // Cython wrapper (bd CoolProp-r9sq.22).
    m.def("HAPropsSI",
          [](const std::string& Output, const std::string& N1, nb::object V1, const std::string& N2, nb::object V2, const std::string& N3,
             nb::object V3) -> nb::object {
              bool s1 = false, s2 = false, s3 = false;
              std::vector<double> v1 = _to_vec(V1, s1), v2 = _to_vec(V2, s2), v3 = _to_vec(V3, s3);
              bool any_seq = s1 || s2 || s3;
              bool any_ndarray = (s1 && nb::hasattr(V1, "ndim")) || (s2 && nb::hasattr(V2, "ndim")) || (s3 && nb::hasattr(V3, "ndim"));
              std::size_t n = std::max({v1.size(), v2.size(), v3.size()});
              _broadcast_to(v1, n, "Input1");
              _broadcast_to(v2, n, "Input2");
              _broadcast_to(v3, n, "Input3");
              std::vector<double> out(n);
              for (std::size_t i = 0; i < n; ++i) {
                  out[i] = HumidAir::HAPropsSI(Output, N1, v1[i], N2, v2[i], N3, v3[i]);
              }
              if (!any_seq) {
                  _raise_if_invalid(out[0]);
                  return nb::float_(out[0]);
              }
              if (any_ndarray) {
                  auto* data = new double[n != 0u ? n : 1u];
                  for (std::size_t i = 0; i < n; ++i) {
                      data[i] = out[i];
                  }
                  nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
                  return nb::cast(nb::ndarray<nb::numpy, double>(data, {n}, owner));
              }
              return nb::cast(out);  // list input -> list output
          });
    // HAProps (non-SI humid air) intentionally NOT bound -- removed for v8 (SI-only); use HAPropsSI.
    m.def("HAProps_Aux", [](std::string out_string, double T, double p, double psi_w) {
        std::array<char, 1000> units{};
        double out = HumidAir::HAProps_Aux(out_string.c_str(), T, p, psi_w, units.data());
        return nb::make_tuple(out, std::string(units.data()));
    });
    m.def("cair_sat", &HumidAir::cair_sat);
    m.def("get_mixture_binary_pair_data", &get_mixture_binary_pair_data);
    m.def("set_mixture_binary_pair_data", &set_mixture_binary_pair_data);
    m.def("apply_simple_mixing_rule", &apply_simple_mixing_rule);

    // ---- q2sh: remaining Cython-parity module functions ---------------------
    // Lower-level C++-backed wrappers
    m.def("set_interaction_parameters", &set_interaction_parameters);
    m.def("set_predefined_mixtures", &set_predefined_mixtures);
    m.def("get_mixture_binary_pair_pcsaft", &get_mixture_binary_pair_pcsaft);
    m.def("set_mixture_binary_pair_pcsaft", &set_mixture_binary_pair_pcsaft);
    // Higher-level convenience wrappers (mirror the Cython module surface; the
    // lower-level get_global_param_string / get_fluid_param_string stay above).
    m.def("FluidsList", []() { return _split_str(get_global_param_string("FluidsList"), ','); });
    m.def("get_aliases", [](const std::string& Fluid) { return _split_str(get_fluid_param_string(Fluid, "aliases_bar"), '|'); });
    m.def("get_REFPROPname", [](const std::string& Fluid) { return get_fluid_param_string(Fluid, "REFPROP_name"); });
    m.def("get_BibTeXKey", [](const std::string& Fluid, const std::string& key) { return get_fluid_param_string(Fluid, "BibTeX-" + key); });
    m.def("get_errstr", []() { return get_global_param_string("errstring"); });
    // set_reference_state: keep the low-level D/S binds above and add the unified
    // (FluidName, *args) dispatcher that the Cython module exposes.
    m.def("set_reference_state", [](const std::string& FluidName, nb::args args) {
        if (args.size() == 1) {
            set_reference_stateS(FluidName, nb::cast<std::string>(args[0]));
        } else if (args.size() == 4) {
            // Coerce via the Python __float__ protocol (nb::float_ -> PyNumber_Float)
            // so numpy scalars / 1-element arrays (e.g. a density from PropsSI) are
            // accepted, matching the legacy Cython `<double>obj` coercion. nb::cast<double>
            // would raise std::bad_cast on a non-float object (bd CoolProp-r9sq.16).
            // Take an nb::handle explicitly: MSVC won't function-style-cast an
            // args[] accessor straight to nb::float_, but accessor->handle is implicit.
            auto _as_double = [](nb::handle h) { return nb::cast<double>(nb::float_(h)); };
            set_reference_stateD(FluidName, _as_double(args[0]), _as_double(args[1]), _as_double(args[2]), _as_double(args[3]));
        } else {
            throw std::invalid_argument("Invalid number of inputs to set_reference_state");
        }
    });
}

#    if defined(COOLPROP_NANOBIND_MODULE)
NB_MODULE(CoolProp, m) {
    init_CoolProp(m);
    // Export the State C-ABI table so the frozen Cython `State` shim can forward
    // through it (PDSim cimport compatibility without a Cython AbstractState).
    m.attr("_capi") = nb::capsule(&g_state_capi, "CoolProp._capi");
}
#    endif

#endif
