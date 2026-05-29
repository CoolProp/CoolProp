// P+X flash preconditioner — standalone unit (CoolProp-cdj prep, Phase 3a)
//
// Two facilities:
//   1. ensure_melt_curve_table(H) / get_melt_curve_table(H) — lazily builds a
//      per-fluid melt-curve table indexed by pressure (single-valued in p even
//      for water, whose melt curve is multi-valued in T).  Cached per-backend
//      via a static map + mutex inside the .cpp, so the public header for
//      HelmholtzEOSMixtureBackend does not need to change.
//   2. regime_classified_rho_seed(H, T, p) — returns a positive finite rho_seed
//      (mol/m^3) chosen by region (supercritical / subcritical liquid /
//      subcritical vapor).  Builds the melt curve table on first call (lazy).
//      Cost target: <= ~200 ns per call in Release builds.  Always finite + > 0.
//
// This unit is NOT yet wired into the production flash path; Phase 3b will
// integrate it into HSU_P_flash_singlephase_Brent.  This phase ships the
// preconditioner as a separable, unit-tested module.

#ifndef COOLPROP_PX_PRECONDITIONER_H
#define COOLPROP_PX_PRECONDITIONER_H

#include <vector>

namespace CoolProp {

class HelmholtzEOSMixtureBackend;

namespace pxprecond {

/// Per-fluid melt-curve density table, p-indexed (single-valued in p for all
/// fluids, including water whose melt curve is multi-valued in T).  Built
/// lazily on first use.
struct MeltCurveTable
{
    std::vector<double> p;         ///< log-spaced from p_triple to a useful upper bound; monotone increasing
    std::vector<double> T_melt;    ///< melting temperature at each p_i
    std::vector<double> rho_melt;  ///< liquid-side molar density at each (T_melt(p_i), p_i)
    bool has_data = false;         ///< false if the fluid has no melting line / build failed
    bool initialized = false;      ///< guard against repeated init

    /// Linear-interpolate rho_melt at the query pressure.  Clamps p to the
    /// table's p range; returns NaN if !has_data.  Cost: O(log N) binary search
    /// + O(1) interp.
    double rho_at(double p_query) const;
};

/// Lazily build (or no-op if already built) the per-fluid melt-curve table for
/// H.  The table is stored in a static map keyed by H's address, protected by a
/// mutex.  If H.has_melting_line() is false, sets has_data=false and returns —
/// no entries.
void ensure_melt_curve_table(HelmholtzEOSMixtureBackend& H);

/// Accessor: returns a const reference to the melt-curve table on H (calling
/// ensure_melt_curve_table first if needed).
const MeltCurveTable& get_melt_curve_table(HelmholtzEOSMixtureBackend& H);

/// Regime-classified seed for the inner density Newton at (T, p), single-phase.
/// Returns a positive rho_seed (mol/m^3) chosen by region:
///   - supercritical (T > Tc): ideal-gas density p/(RT) for high T or low p,
///     blended toward rho_c near the critical point.
///   - subcritical liquid (T < T_sat(p)): log-p interpolation between
///     rho_sat,L(T) and rho_melt(p) when the melt table has data, else
///     compressibility extrapolation.
///   - subcritical vapor (T > T_sat(p)): rho_sat,V(T) blended toward p/(RT) on
///     p/p_sat.
/// Cost target: <= 200 ns per call (Release).  Always finite + positive — if
/// any internal computation yields NaN/Inf/<=0, falls back to p/(RT).  Builds
/// the melt curve table on first call (lazy).
double regime_classified_rho_seed(HelmholtzEOSMixtureBackend& H, double T, double p);

}  // namespace pxprecond
}  // namespace CoolProp

#endif  // COOLPROP_PX_PRECONDITIONER_H
