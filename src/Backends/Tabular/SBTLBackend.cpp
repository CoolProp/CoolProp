#if !defined(NO_TABULAR_BACKENDS)

#include "SBTLBackend.h"
#include "DataStructures.h"

namespace CoolProp {

// ---------------------------------------------------------------------------
// Helper: get pointer to a table field by parameter key
// ---------------------------------------------------------------------------
static std::vector<std::vector<double>>* sbtl_get_field(SinglePhaseGriddedTableData& table, parameters param) {
    switch (param) {
        case iT:      return &table.T;
        case iP:      return &table.p;
        case iDmolar: return &table.rhomolar;
        case iSmolar: return &table.smolar;
        case iHmolar: return &table.hmolar;
        case iUmolar: return &table.umolar;
        default:      throw ValueError("Invalid param in sbtl_get_field");
    }
}

// ---------------------------------------------------------------------------
// DU table set_limits() implementations
// ---------------------------------------------------------------------------

/**
 * Liquid DU table bounds:
 *   x (D): D_crit  →  D_L(T_triple)@P_max   [compressed liquid at coldest+densest]
 *   y (U): U_L(T_triple)  →  U(D_crit, T_max_ext)
 *
 * Many cells in this rectangle will lie in the two-phase or supercritical region;
 * build() leaves them as _HUGE holes, and build_sbtl_coeffs marks them invalid.
 */
void DULiquidTable::set_limits() {
    if (!AS) throw ValueError("DULiquidTable: AS is not set");
    const double T_triple = static_cast<double>(std::max(AS->Ttriple(), AS->Tmin()));

    // Minimum density: D_crit (lower bound of the liquid-like region)
    xmin = static_cast<double>(AS->rhomolar_critical());

    // Maximum density: saturated liquid at triple point + 10% margin.
    // Using PT_INPUTS at (T_triple, P_max) would cross the ice-liquid boundary for water
    // (melting temperature rises above T_triple at high pressure), so use the saturation
    // density at triple point as a safe upper bound; the 10% margin catches moderately
    // compressed liquid states.
    AS->update(QT_INPUTS, 0, T_triple);
    xmax = static_cast<double>(AS->rhomolar()) * 1.10;

    // Minimum U: saturated liquid at triple point
    ymin = static_cast<double>(AS->umolar()) - 50.0;

    // Maximum U: hot state at D_crit and the extended-temperature upper limit
    // (matches the pT table's T_max * 1.499 ceiling)
    double T_max_ext = static_cast<double>(AS->Tmax()) * 1.499;
    CoolPropDbl v1, v2;
    input_pairs pair = generate_update_pair(iDmolar, AS->rhomolar_critical(), iT, T_max_ext, v1, v2);
    AS->update(pair, v1, v2);
    ymax = static_cast<double>(AS->umolar()) * 1.01;
}

/**
 * Gas DU table bounds:
 *   x (D): D_V(T_triple)  →  D_crit   [log-spaced: ratio ~66000 for water]
 *   y (U): U_V(T_triple)  →  U(D_V(T_triple), T_max_ext)
 */
void DUGasTable::set_limits() {
    if (!AS) throw ValueError("DUGasTable: AS is not set");
    const double T_triple = static_cast<double>(std::max(AS->Ttriple(), AS->Tmin()));

    // Minimum density: saturated vapor at triple point
    AS->update(QT_INPUTS, 1, T_triple);
    xmin = static_cast<double>(AS->rhomolar()) * 0.999;

    // Maximum density: D_crit
    xmax = static_cast<double>(AS->rhomolar_critical());

    // Minimum U: saturated vapor at triple point (also computes xmin state)
    ymin = static_cast<double>(AS->umolar()) - 50.0;

    // Maximum U: hot dilute gas at the minimum density and extended T_max
    double T_max_ext = static_cast<double>(AS->Tmax()) * 1.499;
    CoolPropDbl v1, v2;
    input_pairs pair = generate_update_pair(iDmolar, xmin, iT, T_max_ext, v1, v2);
    AS->update(pair, v1, v2);
    ymax = static_cast<double>(AS->umolar()) * 1.01;
}

// ---------------------------------------------------------------------------
// Helper: cubic Lagrange basis polynomial coefficients for 4 non-uniform nodes
// ---------------------------------------------------------------------------
//
// Given 4 node positions px[0..3], computes cx[k][m] = coefficient of ξ^m in
// L_k(ξ), for k=0..3 and m=0..3.  L_k is the unique cubic polynomial satisfying
// L_k(px[k])=1 and L_k(px[m])=0 for m≠k.
//
// Algorithm: for each k, build the product ∏_{m≠k}(ξ - px[m]) by polynomial
// multiplication (working from high degree to low to avoid in-place aliasing),
// then divide by the scalar denominator ∏_{m≠k}(px[k] - px[m]).
static void lagrange_cubic_basis(const double px[4], double cx[4][4]) {
    for (int k = 0; k < 4; ++k) {
        double poly[4] = {1.0, 0.0, 0.0, 0.0};  // running polynomial (starts as 1)
        double denom = 1.0;
        for (int m = 0; m < 4; ++m) {
            if (m == k) continue;
            denom *= (px[k] - px[m]);
            // Multiply poly by (ξ - px[m]): new_poly[d] = old_poly[d-1] - px[m]*old_poly[d]
            // Iterate high-to-low so that poly[d-1] is still the old value when poly[d] is updated.
            for (int d = 3; d >= 1; --d)
                poly[d] = poly[d - 1] - px[m] * poly[d];
            poly[0] *= -px[m];
        }
        for (int d = 0; d < 4; ++d) cx[k][d] = poly[d] / denom;
    }
}

/**
 * @brief Build SBTL polynomial coefficients from node values.
 *
 * @param table    The grid table to process.
 * @param coeffs   Output coefficient array (sized (Nx-1)×(Ny-1)).
 * @param nstencil Stencil width: 3 (degree-2, 9 coefficients per cell, default)
 *                 or 4 (degree-3, 16 coefficients per cell).
 *
 * For nstencil=3 (bi-quadratic):
 *   Stencil nodes in x: i-1, i, i+1  (normalized positions: -rx, 0, 1)
 *   Coefficient layout: alpha[m*3+n]  for m,n ∈ {0,1,2}
 *   Boundary cells: i=0 or j=0 are invalid.
 *
 * For nstencil=4 (bi-cubic):
 *   Stencil nodes in x: i-1, i, i+1, i+2  (normalized: -rx, 0, 1, 1+sx)
 *   Coefficient layout: alpha[m*4+n]  for m,n ∈ {0,1,2,3}
 *   Boundary cells: i=0, i≥Nx-2, j=0, or j≥Ny-2 are invalid.
 *
 * Passes:
 *   1. Validity: boundary cells, non-finite stencil nodes, phase-boundary check.
 *   2. Coefficients: Lagrange basis products.
 *   3. Subdomain decomposition: assign liquid/vapor alternates to invalid cells.
 */
void SBTLBackend::build_sbtl_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs, int nstencil) {
    if (!coeffs.empty()) {
        return;
    }

    const int param_count = 6;
    parameters param_list[param_count] = {iDmolar, iT, iSmolar, iHmolar, iP, iUmolar};

    // (Nx-1) × (Ny-1) cells
    coeffs.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));

    // For nstencil=4, tracks whether each valid cell uses degree-3 (true) or
    // fell back to degree-2 (false).  Allocated regardless of nstencil to
    // keep the second-pass indexing uniform.
    std::vector<std::vector<bool>> use_cubic(table.Nx - 1, std::vector<bool>(table.Ny - 1, false));

    // Returns true when the (hi_off+2)×(hi_off+2) stencil rooted at cell (i,j)
    // is entirely finite and does not cross a phase boundary.
    // hi_off=1 → 3×3, hi_off=2 → 4×4.
    // Caller must guarantee i≥1, j≥1, i+hi_off < Nx, j+hi_off < Ny.
    auto stencil_ok = [&](std::size_t i, std::size_t j, int hi_off) -> bool {
        for (int k = 0; k < param_count; ++k) {
            parameters param = param_list[k];
            if (param == table.xkey || param == table.ykey) continue;
            std::vector<std::vector<double>>* f = sbtl_get_field(table, param);
            for (int di = -1; di <= hi_off; ++di)
                for (int dj = -1; dj <= hi_off; ++dj)
                    if (!ValidNumber((*f)[i + di][j + dj])) return false;
        }
        // Phase-boundary check: the stencil must not straddle the saturation dome.
        // Skip for DU tables (density is a grid axis — no cross-phase stencil issue).
        //
        // Near the critical point, the liquid/gas density ratio can be as small as
        // 1.5-2 (vs the old 100× threshold), so ratio-based detection fails there.
        // Instead, for SUBCRITICAL stencils, check that all nodes are on the SAME
        // side of rho_crit (liquid vs vapor).
        //
        // For SUPERCRITICAL stencils (any T-node ≥ T_crit for pT tables, or any
        // P-node ≥ P_crit for pH tables), the saturation curve doesn't exist —
        // fall back to a 100× ratio safety check instead.
        if (table.xkey != iDmolar && table.ykey != iDmolar) {
            if (table.AS) {
                const double rho_crit = static_cast<double>(table.AS->rhomolar_critical());

                // Determine whether the stencil is in the subcritical region.
                // Subcritical for pT tables:  all T-nodes < T_crit
                // Subcritical for pH tables:  all P-nodes < P_crit
                bool use_phase_check = false;
                if (table.xkey == iT) {
                    // pT table: subcritical if every T-stencil node < T_crit
                    const double T_crit = static_cast<double>(table.AS->T_critical());
                    use_phase_check = true;
                    for (int di = -1; di <= hi_off && use_phase_check; ++di)
                        if (table.xvec[i + di] >= T_crit) use_phase_check = false;
                } else if (table.ykey == iP) {
                    // pH (or similar) table: subcritical if every P-stencil node < P_crit
                    const double P_crit = static_cast<double>(table.AS->p_critical());
                    use_phase_check = true;
                    for (int dj = -1; dj <= hi_off && use_phase_check; ++dj)
                        if (table.yvec[j + dj] >= P_crit) use_phase_check = false;
                }

                if (use_phase_check) {
                    // Subcritical stencil: all nodes must be in the same phase.
                    bool has_liquid = false, has_vapor = false;
                    for (int di = -1; di <= hi_off; ++di) {
                        for (int dj = -1; dj <= hi_off; ++dj) {
                            const double v = table.rhomolar[i + di][j + dj];
                            if (v > rho_crit) has_liquid = true;
                            else has_vapor = true;
                        }
                    }
                    if (has_liquid && has_vapor) return false;
                } else {
                    // Supercritical stencil: no phase boundary — density ratio safety check
                    double rho_min = 1e300, rho_max = 0.0;
                    for (int di = -1; di <= hi_off; ++di) {
                        for (int dj = -1; dj <= hi_off; ++dj) {
                            double v = table.rhomolar[i + di][j + dj];
                            if (v < rho_min) rho_min = v;
                            if (v > rho_max) rho_max = v;
                        }
                    }
                    if (rho_max > 100.0 * rho_min) return false;
                }
            } else {
                // AS unavailable: use density ratio check
                double rho_min = 1e300, rho_max = 0.0;
                for (int di = -1; di <= hi_off; ++di) {
                    for (int dj = -1; dj <= hi_off; ++dj) {
                        double v = table.rhomolar[i + di][j + dj];
                        if (v < rho_min) rho_min = v;
                        if (v > rho_max) rho_max = v;
                    }
                }
                if (rho_max > 100.0 * rho_min) return false;
            }
        }
        return true;
    };

    // -----------------------------------------------------------------------
    // First pass: determine cell validity
    //
    // For nstencil=4: try 4×4 stencil first (degree-3); if it fails due to a
    // phase-boundary node, fall back to 3×3 (degree-2) before marking invalid.
    // This prevents the saturation-dome cells from becoming fully invalid and
    // being redirected to distant alternates that cause large extrapolation errors.
    // -----------------------------------------------------------------------
    for (std::size_t i = 0; i < table.Nx - 1; ++i) {
        for (std::size_t j = 0; j < table.Ny - 1; ++j) {
            // Minimum boundary for a 3×3 stencil: need node i-1 (i≥1) and j-1 (j≥1)
            if (i == 0 || j == 0) {
                coeffs[i][j].set_invalid();
                continue;
            }

            bool cell_valid = false;

            if (nstencil == 4) {
                // Try 4×4 (degree-3) first — requires nodes up to i+2, j+2
                bool can_4x4 = (i + 2 < table.Nx) && (j + 2 < table.Ny);
                if (can_4x4 && stencil_ok(i, j, 2)) {
                    use_cubic[i][j] = true;
                    cell_valid = true;
                } else if (stencil_ok(i, j, 1)) {
                    // 4×4 failed or out-of-bounds — fall back to 3×3 (degree-2)
                    use_cubic[i][j] = false;
                    cell_valid = true;
                }
            } else {
                // nstencil=3: 3×3 stencil only
                cell_valid = stencil_ok(i, j, 1);
            }

            if (cell_valid) {
                coeffs[i][j].set_valid();
                coeffs[i][j].dx_dxhat = table.xvec[i + 1] - table.xvec[i];
                coeffs[i][j].dy_dyhat = table.yvec[j + 1] - table.yvec[j];
            } else {
                coeffs[i][j].set_invalid();
            }
        }
    }

    // -----------------------------------------------------------------------
    // Second pass: compute polynomial coefficients for valid cells
    // -----------------------------------------------------------------------
    for (int k = 0; k < param_count; ++k) {
        parameters param = param_list[k];
        if (param == table.xkey || param == table.ykey) {
            continue;
        }

        std::vector<std::vector<double>>* f = sbtl_get_field(table, param);

        for (std::size_t i = 0; i < table.Nx - 1; ++i) {
            for (std::size_t j = 0; j < table.Ny - 1; ++j) {
                if (!coeffs[i][j].valid()) continue;

                if (nstencil == 4 && use_cubic[i][j]) {
                    // ---- Degree-3 (bi-cubic): 4×4 stencil, 16 coefficients ----
                    // Normalized positions of the 4 stencil nodes in x: {-rx, 0, 1, 1+sx}
                    double rx = (table.xvec[i]     - table.xvec[i - 1]) / (table.xvec[i + 1] - table.xvec[i]);
                    double sx = (table.xvec[i + 2] - table.xvec[i + 1]) / (table.xvec[i + 1] - table.xvec[i]);
                    double ry = (table.yvec[j]     - table.yvec[j - 1]) / (table.yvec[j + 1] - table.yvec[j]);
                    double sy = (table.yvec[j + 2] - table.yvec[j + 1]) / (table.yvec[j + 1] - table.yvec[j]);

                    double px[4] = {-rx, 0.0, 1.0, 1.0 + sx};
                    double py[4] = {-ry, 0.0, 1.0, 1.0 + sy};

                    // z[k][l]: k=0 → node i-1, k=1 → node i, k=2 → node i+1, k=3 → node i+2
                    double z[4][4];
                    for (int di = -1; di <= 2; ++di)
                        for (int dj = -1; dj <= 2; ++dj)
                            z[di + 1][dj + 1] = (*f)[i + di][j + dj];

                    // Cubic Lagrange basis: cx[k][m] = coeff of ξ^m in L_k(ξ)
                    double cx[4][4], cy[4][4];
                    lagrange_cubic_basis(px, cx);
                    lagrange_cubic_basis(py, cy);

                    // a_{mn} = Σ_{k,l} z[k][l] · cx[k][m] · cy[l][n]
                    std::vector<double> alpha(16, 0.0);
                    for (int m = 0; m < 4; ++m) {
                        for (int n = 0; n < 4; ++n) {
                            double val = 0.0;
                            for (int kk = 0; kk < 4; ++kk) {
                                if (cx[kk][m] == 0.0) continue;
                                for (int ll = 0; ll < 4; ++ll)
                                    val += z[kk][ll] * cx[kk][m] * cy[ll][n];
                            }
                            alpha[m * 4 + n] = val;
                        }
                    }
                    coeffs[i][j].set(param, alpha);

                } else {
                    // ---- Degree-2 (bi-quadratic): 3×3 stencil, 9 coefficients ----
                    // Used for nstencil=3 tables (DU tables) and as the degree-3
                    // fallback for nstencil=4 cells near the saturation dome.
                    double rx = (table.xvec[i] - table.xvec[i - 1]) / (table.xvec[i + 1] - table.xvec[i]);
                    double ry = (table.yvec[j] - table.yvec[j - 1]) / (table.yvec[j + 1] - table.yvec[j]);

                    // z[k][l]: k=0 → node i-1, k=1 → node i, k=2 → node i+1
                    double z[3][3];
                    for (int di = -1; di <= 1; ++di)
                        for (int dj = -1; dj <= 1; ++dj)
                            z[di + 1][dj + 1] = (*f)[i + di][j + dj];

                    // Lagrange basis coefficients: cx[k][m] = coeff of ξ^m in L_k(ξ)
                    double cx[3][3], cy[3][3];
                    cx[0][0] = 0.0;  cx[0][1] = -1.0 / (rx * (rx + 1.0));  cx[0][2] = 1.0 / (rx * (rx + 1.0));
                    cx[1][0] = 1.0;  cx[1][1] = (1.0 - rx) / rx;           cx[1][2] = -1.0 / rx;
                    cx[2][0] = 0.0;  cx[2][1] = rx / (1.0 + rx);           cx[2][2] = 1.0 / (1.0 + rx);

                    cy[0][0] = 0.0;  cy[0][1] = -1.0 / (ry * (ry + 1.0));  cy[0][2] = 1.0 / (ry * (ry + 1.0));
                    cy[1][0] = 1.0;  cy[1][1] = (1.0 - ry) / ry;           cy[1][2] = -1.0 / ry;
                    cy[2][0] = 0.0;  cy[2][1] = ry / (1.0 + ry);           cy[2][2] = 1.0 / (1.0 + ry);

                    // a_{mn} = Σ_{k,l} z[k][l] · cx[k][m] · cy[l][n]
                    std::vector<double> alpha(9, 0.0);
                    for (int m = 0; m < 3; ++m) {
                        for (int n = 0; n < 3; ++n) {
                            double val = 0.0;
                            for (int kk = 0; kk < 3; ++kk) {
                                if (cx[kk][m] == 0.0) continue;
                                for (int ll = 0; ll < 3; ++ll)
                                    val += z[kk][ll] * cx[kk][m] * cy[ll][n];
                            }
                            alpha[m * 3 + n] = val;
                        }
                    }
                    coeffs[i][j].set(param, alpha);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Third pass: subdomain decomposition — give each invalid cell two alternates.
    //
    // Primary alternate   (alt_i/alt_j):   nearest valid LIQUID cell
    //   → found by scanning toward lower-T (i.e., i decreasing) then all directions
    //
    // Secondary alternate (alt_i2/alt_j2): nearest valid VAPOR cell
    //   → found by scanning toward higher-T (i.e., i increasing) then all directions
    //
    // A valid cell is "liquid" if its center-node density > rho_crit, "vapor" otherwise.
    // Since valid cells passed the rho_max/rho_min < 100 test and the liquid/vapor
    // density ratio for real fluids is >> 100, every valid cell is entirely in one phase.
    // -----------------------------------------------------------------------
    const double rho_crit_val = table.AS ? static_cast<double>(table.AS->rhomolar_critical()) : 1e20;

    auto is_liquid_valid = [&](std::size_t ci, std::size_t cj) -> bool {
        return coeffs[ci][cj].valid()
               && ValidNumber(table.rhomolar[ci][cj])
               && table.rhomolar[ci][cj] > rho_crit_val;
    };
    auto is_vapor_valid = [&](std::size_t ci, std::size_t cj) -> bool {
        return coeffs[ci][cj].valid()
               && ValidNumber(table.rhomolar[ci][cj])
               && table.rhomolar[ci][cj] < rho_crit_val;
    };

    // Expanded search patterns: 8 directions, each represented as (di, dj).
    // For liquid (lower-T): search -x first; for vapor (higher-T): search +x first.
    static const int liq_di[8] = {-1,  1,  0,  0, -1,  1, -1,  1};
    static const int liq_dj[8] = { 0,  0, -1,  1, -1, -1,  1,  1};
    static const int vap_di[8] = { 1, -1,  0,  0,  1, -1,  1, -1};
    static const int vap_dj[8] = { 0,  0, -1,  1, -1, -1,  1,  1};

    for (std::size_t i = 0; i < table.Nx - 1; ++i) {
        for (std::size_t j = 0; j < table.Ny - 1; ++j) {
            if (coeffs[i][j].valid()) continue;

            // Scan increasing radii until we find valid cells in both phases.
            // For efficiency, cap the search at Nx (should always find within a few steps).
            const int maxR = static_cast<int>(std::max(table.Nx, table.Ny));
            bool found_liq = false, found_vap = false;

            for (int r = 1; r < maxR && !(found_liq && found_vap); ++r) {
                // 8-direction scan at this radius
                for (int d = 0; d < 8; ++d) {
                    int il = static_cast<int>(i) + liq_di[d] * r;
                    int jl = static_cast<int>(j) + liq_dj[d] * r;
                    if (!found_liq && il >= 0 && il < static_cast<int>(table.Nx) - 1
                        && jl >= 0 && jl < static_cast<int>(table.Ny) - 1
                        && is_liquid_valid(static_cast<std::size_t>(il), static_cast<std::size_t>(jl))) {
                        coeffs[i][j].set_alternate(static_cast<std::size_t>(il), static_cast<std::size_t>(jl));
                        found_liq = true;
                    }

                    int iv = static_cast<int>(i) + vap_di[d] * r;
                    int jv = static_cast<int>(j) + vap_dj[d] * r;
                    if (!found_vap && iv >= 0 && iv < static_cast<int>(table.Nx) - 1
                        && jv >= 0 && jv < static_cast<int>(table.Ny) - 1
                        && is_vapor_valid(static_cast<std::size_t>(iv), static_cast<std::size_t>(jv))) {
                        coeffs[i][j].set_alternate2(static_cast<std::size_t>(iv), static_cast<std::size_t>(jv));
                        found_vap = true;
                    }

                    if (found_liq && found_vap) break;
                }
            }

            // If no phase-specific alternate found (e.g. supercritical-only table),
            // fall back to any valid neighbour for the missing phase.
            if (!found_liq && !found_vap) {
                throw ValueError(format("SBTLBackend: cell (%zu,%zu) has no valid neighbors", i, j));
            }
            if (!found_liq) {
                // Use vapor alternate as fallback for liquid too
                std::size_t ia = i, ja = j;
                coeffs[i][j].get_alternate2(ia, ja);
                coeffs[i][j].set_alternate(ia, ja);
            }
            if (!found_vap) {
                // Use liquid alternate as fallback for vapor too
                std::size_t ia = i, ja = j;
                coeffs[i][j].get_alternate(ia, ja);
                coeffs[i][j].set_alternate2(ia, ja);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// DU table builder
// ---------------------------------------------------------------------------

void SBTLBackend::build_du_tables() {
    if (du_tables_built) return;

    du_liquid.AS = this->AS;
    du_liquid.set_limits();
    du_liquid.build(this->AS);
    build_sbtl_coeffs(du_liquid, coeffs_du_liquid);

    du_gas.AS = this->AS;
    du_gas.set_limits();
    du_gas.build(this->AS);
    build_sbtl_coeffs(du_gas, coeffs_du_gas);

    build_sat_cache();
    du_tables_built = true;
}

// ---------------------------------------------------------------------------
// Saturation property cache
// ---------------------------------------------------------------------------

/**
 * Build a precomputed saturation cache for flash_DmolarUmolar.
 *
 * Resamples the existing PureFluidSaturationTableData onto N_SAT log-spaced
 * pressure nodes, enabling O(1) index computation during flash (no binary search).
 * Reduces per-saturation-lookup cost from ~300 ns (EOS call) to ~5–10 ns.
 *
 * Per SBTL PDF Appendix A4, dedicated saturation splines v'(u), v''(p), u'(T) are
 * the standard approach for fast two-phase calculations from (v,u).
 */
void SBTLBackend::build_sat_cache() {
    if (sat_cache_built) return;

    PureFluidSaturationTableData& sat = dataset->pure_saturation;

    const int N_SAT = 1000;
    // Use the inner range of the saturation pressure array to avoid edge effects
    const double p_lo = sat.pL.front() * 1.0001;
    const double p_hi = sat.pL.back()  * (1.0 - 1e-5);

    _sat_log_p_min = std::log(p_lo);
    const double log_p_max = std::log(p_hi);
    const double dlogp = (log_p_max - _sat_log_p_min) / (N_SAT - 1);
    _sat_inv_dlogp = 1.0 / dlogp;

    _sat_logp.resize(N_SAT);
    _sat_DL.resize(N_SAT);
    _sat_DV.resize(N_SAT);
    _sat_UL.resize(N_SAT);
    _sat_UV.resize(N_SAT);
    _sat_T.resize(N_SAT);

    const std::size_t Nsat = sat.N;

    for (int k = 0; k < N_SAT; ++k) {
        const double log_p = _sat_log_p_min + k * dlogp;
        const double p = std::exp(log_p);
        _sat_logp[k] = log_p;

        // Binary search for bracket indices in the saturation table pressure vectors
        std::size_t iL = 0, iV = 0;
        bisect_vector(sat.pL, p, iL);
        bisect_vector(sat.pV, p, iV);
        // Clamp so CubicInterp has a valid 4-point stencil {idx-2, idx-1, idx, idx+1}
        iL = std::max(std::size_t(2), std::min(Nsat - 2, iL));
        iV = std::max(std::size_t(2), std::min(Nsat - 2, iV));

        // Same cubic interpolation as PureFluidSaturationTableData::evaluate
        _sat_DL[k] = std::exp(CubicInterp(sat.logpL, sat.logrhomolarL, iL-2, iL-1, iL, iL+1, log_p));
        _sat_DV[k] = std::exp(CubicInterp(sat.logpV, sat.logrhomolarV, iV-2, iV-1, iV, iV+1, log_p));
        _sat_UL[k] = CubicInterp(sat.logpL, sat.umolarL, iL-2, iL-1, iL, iL+1, log_p);
        _sat_UV[k] = CubicInterp(sat.logpV, sat.umolarV, iV-2, iV-1, iV, iV+1, log_p);
        _sat_T[k]  = CubicInterp(sat.logpL, sat.TL,      iL-2, iL-1, iL, iL+1, log_p);
    }

    // Cache triple-point saturation densities for the always-executed phase bounds check
    _sat_DL_triple = _sat_DL[0];
    _sat_DV_triple = _sat_DV[0];

    sat_cache_built = true;
}

/**
 * O(1) saturation liquid property lookup at pressure p.
 *
 * Log-spaced pressure indexing lets us compute the array slot analytically — no
 * binary search.  Linear interpolation gives the value.
 */
double SBTLBackend::fast_sat_lookup(parameters param, double p) const {
    const double log_p = std::log(p);
    const double idx_f = (log_p - _sat_log_p_min) * _sat_inv_dlogp;
    const int N = static_cast<int>(_sat_logp.size());
    const int idx = std::max(0, std::min(N - 2, static_cast<int>(idx_f)));
    const double t = idx_f - static_cast<double>(idx);

    switch (param) {
        case iDmolar: return (1.0 - t) * _sat_DL[idx] + t * _sat_DL[idx + 1];
        case iUmolar: return (1.0 - t) * _sat_UL[idx] + t * _sat_UL[idx + 1];
        case iT:      return (1.0 - t) * _sat_T[idx]  + t * _sat_T[idx + 1];
        default:      throw ValueError("fast_sat_lookup: unsupported param");
    }
}

/// O(1) saturation vapor property lookup at pressure p.
double SBTLBackend::fast_sat_lookup_vapor(parameters param, double p) const {
    const double log_p = std::log(p);
    const double idx_f = (log_p - _sat_log_p_min) * _sat_inv_dlogp;
    const int N = static_cast<int>(_sat_logp.size());
    const int idx = std::max(0, std::min(N - 2, static_cast<int>(idx_f)));
    const double t = idx_f - static_cast<double>(idx);

    switch (param) {
        case iDmolar: return (1.0 - t) * _sat_DV[idx] + t * _sat_DV[idx + 1];
        case iUmolar: return (1.0 - t) * _sat_UV[idx] + t * _sat_UV[idx + 1];
        case iT:      return (1.0 - t) * _sat_T[idx]  + t * _sat_T[idx + 1];
        default:      throw ValueError("fast_sat_lookup_vapor: unsupported param");
    }
}

// ---------------------------------------------------------------------------
// Cell finding: phase-aware subdomain selection
// ---------------------------------------------------------------------------

/**
 * Find the correct cell for (x, y), using the phase-appropriate subdomain alternate
 * for cells near the saturation boundary.
 *
 * The saturation table is queried at query time to determine whether (x, y) is on
 * the liquid side (x < x_sat_liquid) or vapor side (x > x_sat_vapor).  The result
 * selects which of the two pre-computed alternates to use:
 *   - Liquid query → primary alternate (set_alternate)   → nearest all-liquid cell
 *   - Vapor  query → secondary alternate (set_alternate2) → nearest all-vapor cell
 *
 * This implements subdomain decomposition at lookup time without needing separate
 * grids: liquid and vapor each have their own "territory" of valid cells, and
 * invalid cells (near the saturation boundary) are redirected to the correct territory.
 *
 * NOTE: Phase detection is only performed when table.ykey == iP (i.e. P is the
 * second coordinate).  For DU tables (ykey == iUmolar) the detection is skipped
 * and the primary alternate is used unconditionally.
 */
void SBTLBackend::find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                                   double x, double y, std::size_t& i, std::size_t& j) {
    // Use SBTL's own coefficient arrays for pT/pH tables, falling back to the
    // passed argument for DU tables (which are already private to SBTL).
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    table.find_native_nearest_good_cell(x, y, i, j);
    if (!coeffs[i][j].valid()) {
        // Determine query phase from the saturation table (for pure fluids at subcritical P).
        // Only applicable when y is pressure (ykey == iP); skip for DU tables where y is energy.
        bool query_is_vapor = false;
        bool phase_known = false;

        if (!is_mixture && table.AS && table.ykey == iP && y < table.AS->p_critical()) {
            try {
                // Use iQ lookup: always sets T_sat_L/T_sat_V and returns true within the
                // saturation table pressure range, regardless of where x lies.
                std::size_t iL = 0, iV = 0;
                CoolPropDbl T_sat_L = 0, T_sat_V = 0;
                bool in_range = dataset->pure_saturation.is_inside(iP, y, iQ, 0.5, iL, iV, T_sat_L, T_sat_V);
                if (in_range) {
                    if (table.xkey == iT) {
                        // x is temperature; for pure fluid T_sat_L ≈ T_sat_V = T_sat(P)
                        if (x < T_sat_L - 1e-6) {
                            query_is_vapor = false;
                            phase_known = true;
                        } else if (x > T_sat_V + 1e-6) {
                            query_is_vapor = true;
                            phase_known = true;
                        }
                    } else if (table.xkey == iHmolar) {
                        // x is enthalpy; need hL_sat and hV_sat at pressure y.
                        // iL and iV are bracket indices from is_inside (bisect on pL/pV).
                        auto& sat = dataset->pure_saturation;
                        std::size_t iVp = std::min(iV + 1, sat.N - 1);
                        std::size_t iLp = std::min(iL + 1, sat.N - 1);
                        if (iVp < 3) iVp = 3;
                        if (iLp < 3) iLp = 3;
                        double logp = std::log(y);
                        double hV_sat = CubicInterp(sat.logpV, sat.hmolarV, iVp - 3, iVp - 2, iVp - 1, iVp, logp);
                        double hL_sat = CubicInterp(sat.logpL, sat.hmolarL, iLp - 3, iLp - 2, iLp - 1, iLp, logp);
                        if (x < hL_sat - 1.0) {
                            query_is_vapor = false;
                            phase_known = true;
                        } else if (x > hV_sat + 1.0) {
                            query_is_vapor = true;
                            phase_known = true;
                        }
                    }
                }
            } catch (...) {
                // Saturation lookup failed; fall through to generic fallback.
            }
        }

        if (phase_known) {
            if (!query_is_vapor) {
                // Liquid phase: use primary alternate (nearest all-liquid cell)
                if (coeffs[i][j].has_valid_neighbor()) {
                    std::size_t ia = i, ja = j;
                    coeffs[ia][ja].get_alternate(ia, ja);
                    i = ia; j = ja;
                    return;
                }
            } else {
                // Vapor phase: use secondary alternate (nearest all-vapor cell)
                if (coeffs[i][j].has_valid_neighbor2()) {
                    std::size_t ia = i, ja = j;
                    coeffs[ia][ja].get_alternate2(ia, ja);
                    i = ia; j = ja;
                    return;
                }
                // If no vapor alternate, fall through to primary
                if (coeffs[i][j].has_valid_neighbor()) {
                    coeffs[i][j].get_alternate(i, j);
                    return;
                }
            }
        } else {
            // Phase unknown (supercritical, DU table, or saturation lookup failed):
            // use primary alternate (which for DU tables always points to a same-phase cell)
            if (coeffs[i][j].has_valid_neighbor()) {
                coeffs[i][j].get_alternate(i, j);
                return;
            }
        }

        throw ValueError(format("SBTLBackend: cell is invalid and has no good neighbors for x = %g, y = %g", x, y));
    }
}

void SBTLBackend::find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                        const parameters variable1, const double value1, const parameters otherkey, const double otherval,
                                        std::size_t& i, std::size_t& j) {
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    table.find_nearest_neighbor(variable1, value1, otherkey, otherval, i, j);
    const CellCoeffs& cell = coeffs[i][j];
    if (!cell.valid()) {
        if (cell.has_valid_neighbor()) {
            cell.get_alternate(i, j);
        } else {
            throw ValueError(format("Cell is invalid and has no good neighbors for x = %g, y = %g", value1, otherval));
        }
    }
}

// ---------------------------------------------------------------------------
// Transport properties: bilinear interpolation
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y, std::size_t i,
                                                    std::size_t j) {
    std::vector<std::vector<double>>* f = NULL;
    switch (output) {
        case iconductivity:
            f = &table.cond;
            break;
        case iviscosity:
            f = &table.visc;
            break;
        default:
            throw ValueError(format("invalid output variable to SBTLBackend::evaluate_single_phase_transport"));
    }
    double x1 = table.xvec[i], x2 = table.xvec[i + 1], y1 = table.yvec[j], y2 = table.yvec[j + 1];
    double f11 = (*f)[i][j], f12 = (*f)[i][j + 1], f21 = (*f)[i + 1][j], f22 = (*f)[i + 1][j + 1];
    double val =
      1.0 / ((x2 - x1) * (y2 - y1)) * (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y) + f12 * (x2 - x) * (y - y1) + f22 * (x - x1) * (y - y1));
    switch (output) {
        case iconductivity:
            _conductivity = val;
            break;
        case iviscosity:
            _viscosity = val;
            break;
        default:
            throw ValueError("Invalid output variable in SBTLBackend::evaluate_single_phase_transport");
    }
    return val;
}

// ---------------------------------------------------------------------------
// Polynomial evaluation — degree-2 (9 coefficients) or degree-3 (16 coefficients)
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                          const parameters output, const double x, const double y, const std::size_t i, const std::size_t j) {
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    const CellCoeffs& cell = coeffs[i][j];
    const std::vector<double>& alpha = cell.get(output);

    double xi  = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    double eta = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);

    double val;
    if (alpha.size() == 9) {
        // Degree-2 (bi-quadratic): z = Σ_{m,n=0}^{2} alpha[m*3+n] xi^m eta^n
        double eta2 = eta * eta;
        double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2;
        double B1 = alpha[3] + alpha[4] * eta + alpha[5] * eta2;
        double B2 = alpha[6] + alpha[7] * eta + alpha[8] * eta2;
        val = B0 + B1 * xi + B2 * xi * xi;
    } else {
        // Degree-3 (bi-cubic): z = Σ_{m,n=0}^{3} alpha[m*4+n] xi^m eta^n
        double xi2 = xi * xi, xi3 = xi2 * xi;
        double eta2 = eta * eta, eta3 = eta2 * eta;
        double B0 = alpha[ 0] + alpha[ 1] * eta + alpha[ 2] * eta2 + alpha[ 3] * eta3;
        double B1 = alpha[ 4] + alpha[ 5] * eta + alpha[ 6] * eta2 + alpha[ 7] * eta3;
        double B2 = alpha[ 8] + alpha[ 9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
        double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
        val = B0 + B1 * xi + B2 * xi2 + B3 * xi3;
    }

    switch (output) {
        case iT:      _T = val;        break;
        case iP:      _p = val;        break;
        case iDmolar: _rhomolar = val; break;
        case iSmolar: _smolar = val;   break;
        case iHmolar: _hmolar = val;   break;
        case iUmolar: _umolar = val;   break;
        default:
            throw ValueError("Invalid output variable in SBTLBackend::evaluate_single_phase");
    }
    return val;
}

// ---------------------------------------------------------------------------
// Analytical derivatives dz/dx and dz/dy via chain rule
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                                     parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx,
                                                     std::size_t Ny) {
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    const CellCoeffs& cell = coeffs[i][j];
    const std::vector<double>& alpha = cell.get(output);

    double xi  = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    double eta = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);

    if (Nx == 1 && Ny == 0) {
        if (output == table.xkey) return 1.0;
        if (output == table.ykey) return 0.0;
        double dzdxi;
        if (alpha.size() == 9) {
            double eta2 = eta * eta;
            double B1 = alpha[3] + alpha[4] * eta + alpha[5] * eta2;
            double B2 = alpha[6] + alpha[7] * eta + alpha[8] * eta2;
            dzdxi = B1 + 2.0 * B2 * xi;
        } else {
            double xi2 = xi * xi;
            double eta2 = eta * eta, eta3 = eta2 * eta;
            double B1 = alpha[ 4] + alpha[ 5] * eta + alpha[ 6] * eta2 + alpha[ 7] * eta3;
            double B2 = alpha[ 8] + alpha[ 9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
            double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
            dzdxi = B1 + 2.0 * B2 * xi + 3.0 * B3 * xi2;
        }
        return dzdxi / (table.xvec[i + 1] - table.xvec[i]);
    } else if (Ny == 1 && Nx == 0) {
        if (output == table.ykey) return 1.0;
        if (output == table.xkey) return 0.0;
        double dzdeta;
        if (alpha.size() == 9) {
            double eta2 = eta * eta;
            dzdeta = (alpha[1] + 2.0 * alpha[2] * eta) + (alpha[4] + 2.0 * alpha[5] * eta) * xi
                     + (alpha[7] + 2.0 * alpha[8] * eta) * xi * xi;
        } else {
            double xi2 = xi * xi, xi3 = xi2 * xi;
            double eta2 = eta * eta;
            dzdeta = (alpha[ 1] + 2.0 * alpha[ 2] * eta + 3.0 * alpha[ 3] * eta2)
                   + (alpha[ 5] + 2.0 * alpha[ 6] * eta + 3.0 * alpha[ 7] * eta2) * xi
                   + (alpha[ 9] + 2.0 * alpha[10] * eta + 3.0 * alpha[11] * eta2) * xi2
                   + (alpha[13] + 2.0 * alpha[14] * eta + 3.0 * alpha[15] * eta2) * xi3;
        }
        return dzdeta / (table.yvec[j + 1] - table.yvec[j]);
    } else {
        throw ValueError("Invalid Nx/Ny in SBTLBackend::evaluate_single_phase_derivative");
    }
}

// ---------------------------------------------------------------------------
// Inversion: quadratic formula (degree-2) or Newton from quadratic seed (degree-3)
// ---------------------------------------------------------------------------

void SBTLBackend::invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                        parameters other_key, double other, double y, std::size_t i, std::size_t j) {
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    const CellCoeffs& cell = coeffs[i][j];
    const std::vector<double>& alpha = cell.get(other_key);

    double eta  = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);
    double eta2 = eta * eta;

    double xi;
    if (alpha.size() == 9) {
        // Degree-2: solve A*xi^2 + B*xi + C = 0 analytically
        double A = alpha[6] + alpha[7] * eta + alpha[8] * eta2;
        double B = alpha[3] + alpha[4] * eta + alpha[5] * eta2;
        double C = alpha[0] + alpha[1] * eta + alpha[2] * eta2 - other;
        if (std::abs(A) < 1e-14 * std::abs(B)) {
            if (std::abs(B) > 1e-14)
                xi = -C / B;
            else
                throw ValueError("Degenerate quadratic in SBTLBackend::invert_single_phase_x");
        } else {
            double disc = B * B - 4.0 * A * C;
            if (disc < 0.0) disc = 0.0;
            double sqrtD = std::sqrt(disc);
            double xi1 = (-B + sqrtD) / (2.0 * A);
            double xi2 = (-B - sqrtD) / (2.0 * A);
            xi = (std::min(std::abs(xi1), std::abs(xi1 - 1.0)) <= std::min(std::abs(xi2), std::abs(xi2 - 1.0))) ? xi1 : xi2;
        }
    } else {
        // Degree-3: solve D3*xi^3 + A*xi^2 + B*xi + C = 0 via Newton from quadratic seed
        double eta3 = eta2 * eta;
        double D3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
        double A  = alpha[ 8] + alpha[ 9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
        double B  = alpha[ 4] + alpha[ 5] * eta + alpha[ 6] * eta2 + alpha[ 7] * eta3;
        double C  = alpha[ 0] + alpha[ 1] * eta + alpha[ 2] * eta2 + alpha[ 3] * eta3 - other;

        // Quadratic seed (ignore D3 term)
        double xi_seed;
        if (std::abs(A) < 1e-14 * std::abs(B)) {
            xi_seed = std::abs(B) > 1e-14 ? -C / B : 0.5;
        } else {
            double disc = B * B - 4.0 * A * C;
            if (disc < 0.0) disc = 0.0;
            double sq = std::sqrt(disc);
            double x1 = (-B + sq) / (2.0 * A);
            double x2 = (-B - sq) / (2.0 * A);
            xi_seed = (std::min(std::abs(x1), std::abs(x1 - 1.0)) <= std::min(std::abs(x2), std::abs(x2 - 1.0))) ? x1 : x2;
        }

        // Newton refinement (converges in 2–3 iterations for smooth thermodynamic data)
        xi = xi_seed;
        for (int iter = 0; iter < 6; ++iter) {
            double fval  = ((D3 * xi + A) * xi + B) * xi + C;
            double dfval = (3.0 * D3 * xi + 2.0 * A) * xi + B;
            if (std::abs(dfval) < 1e-30) break;
            double dxi = fval / dfval;
            xi -= dxi;
            if (std::abs(dxi) < 1e-13 * (1.0 + std::abs(xi))) break;
        }
    }

    double val = xi * (table.xvec[i + 1] - table.xvec[i]) + table.xvec[i];

    switch (table.xkey) {
        case iHmolar: _hmolar = val; break;
        case iT:      _T = val;      break;
        default:
            throw ValueError("Invalid x-key in SBTLBackend::invert_single_phase_x");
    }
}

void SBTLBackend::invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs_arg,
                                        parameters other_key, double other, double x, std::size_t i, std::size_t j) {
    const auto& coeffs = sbtl_coeffs_for(table, coeffs_arg);
    const CellCoeffs& cell = coeffs[i][j];
    const std::vector<double>& alpha = cell.get(other_key);

    double xi  = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    double xi2 = xi * xi;

    double eta;
    if (alpha.size() == 9) {
        // Degree-2: solve A*eta^2 + B*eta + C = 0 analytically
        double A = alpha[2] + alpha[5] * xi + alpha[8] * xi2;
        double B = alpha[1] + alpha[4] * xi + alpha[7] * xi2;
        double C = alpha[0] + alpha[3] * xi + alpha[6] * xi2 - other;
        if (std::abs(A) < 1e-14 * std::abs(B)) {
            if (std::abs(B) > 1e-14)
                eta = -C / B;
            else
                throw ValueError("Degenerate quadratic in SBTLBackend::invert_single_phase_y");
        } else {
            double disc = B * B - 4.0 * A * C;
            if (disc < 0.0) disc = 0.0;
            double sqrtD = std::sqrt(disc);
            double e1 = (-B + sqrtD) / (2.0 * A);
            double e2 = (-B - sqrtD) / (2.0 * A);
            eta = (std::min(std::abs(e1), std::abs(e1 - 1.0)) <= std::min(std::abs(e2), std::abs(e2 - 1.0))) ? e1 : e2;
        }
    } else {
        // Degree-3: solve D3*eta^3 + A*eta^2 + B*eta + C = 0 via Newton from quadratic seed
        double xi3 = xi2 * xi;
        double D3 = alpha[ 3] + alpha[ 7] * xi + alpha[11] * xi2 + alpha[15] * xi3;
        double A  = alpha[ 2] + alpha[ 6] * xi + alpha[10] * xi2 + alpha[14] * xi3;
        double B  = alpha[ 1] + alpha[ 5] * xi + alpha[ 9] * xi2 + alpha[13] * xi3;
        double C  = alpha[ 0] + alpha[ 4] * xi + alpha[ 8] * xi2 + alpha[12] * xi3 - other;

        // Quadratic seed (ignore D3 term)
        double eta_seed;
        if (std::abs(A) < 1e-14 * std::abs(B)) {
            eta_seed = std::abs(B) > 1e-14 ? -C / B : 0.5;
        } else {
            double disc = B * B - 4.0 * A * C;
            if (disc < 0.0) disc = 0.0;
            double sq = std::sqrt(disc);
            double e1 = (-B + sq) / (2.0 * A);
            double e2 = (-B - sq) / (2.0 * A);
            eta_seed = (std::min(std::abs(e1), std::abs(e1 - 1.0)) <= std::min(std::abs(e2), std::abs(e2 - 1.0))) ? e1 : e2;
        }

        // Newton refinement
        eta = eta_seed;
        for (int iter = 0; iter < 6; ++iter) {
            double fval  = ((D3 * eta + A) * eta + B) * eta + C;
            double dfval = (3.0 * D3 * eta + 2.0 * A) * eta + B;
            if (std::abs(dfval) < 1e-30) break;
            double deta = fval / dfval;
            eta -= deta;
            if (std::abs(deta) < 1e-13 * (1.0 + std::abs(eta))) break;
        }
    }

    double val = eta * (table.yvec[j + 1] - table.yvec[j]) + table.yvec[j];

    switch (table.ykey) {
        case iP: _p = val; break;
        default:
            throw ValueError("Invalid y-key in SBTLBackend::invert_single_phase_y");
    }
}

// ---------------------------------------------------------------------------
// D,U flash: direct DU-table lookup — no Newton iteration
//
// Given molar density D [mol/m³] and molar internal energy U [J/mol],
// find T and P using dedicated DU-space polynomial tables where T(D,U)
// and P(D,U) are stored as forward evaluations.
//
// Algorithm:
//  1. Phase detection via saturation table: widest saturation density range
//     (at the triple point) → if D is between D_V and D_L at triple point,
//     the state is potentially two-phase.
//
//  2. Two-phase Newton (SBTL Appendix A4): if possibly two-phase, run Newton
//     to find p_sat such that x_v(p) = x_u(p).  If the quality estimate
//     Q ∈ [0.01, 0.99] and both estimates agree, treat as two-phase.
//
//  3. Single-phase: direct DU-table forward evaluation.
//     - Select liquid table (D > D_crit) or gas table (D ≤ D_crit).
//     - find_native_nearest_good_indices → cell (i, j)
//     - T = evaluate_single_phase(..., iT, D, U, i, j)
//     - P = evaluate_single_phase(..., iP, D, U, i, j)
//     - Map (T, P) to the pT table cell for subsequent property lookups.
// ---------------------------------------------------------------------------
void SBTLBackend::flash_DmolarUmolar(CoolPropDbl D_in, CoolPropDbl U_in) {
    const double D_target = static_cast<double>(D_in);
    const double U_target = static_cast<double>(U_in);
    _rhomolar = D_target;
    _umolar = U_target;

    const double D_crit   = static_cast<double>(this->AS->rhomolar_critical());

    // ----- Step 1: Cheap phase detection from cached triple-point densities -----
    // _sat_DL_triple and _sat_DV_triple bound the entire two-phase dome.
    // States outside this density range are guaranteed single-phase.
    bool possibly_two_phase = (D_target > _sat_DV_triple && D_target < _sat_DL_triple);
    bool is_confirmed_two_phase = false;
    double p_sat_found = 0.0;

    if (possibly_two_phase) {
        // ----- Step 2: Two-phase Newton (SBTL Appendix A4) -----
        // Solve  f(p) = x_v(p) - x_u(p) = 0  using fast saturation cache.
        // Each function evaluation costs ~30 ns instead of ~300 ns with sat_eval.
        const double v_target = 1.0 / D_target;
        const double p_lo_sat = std::exp(_sat_log_p_min);  // triple point pressure
        const double p_hi_sat = std::exp(_sat_log_p_min + (_sat_logp.size()-1) / _sat_inv_dlogp);

        // Inline fast f(p): avoids repeated function-call overhead
        auto fast_f = [&](double p_test) -> double {
            const double DL = fast_sat_lookup(iDmolar, p_test);
            const double DV = fast_sat_lookup_vapor(iDmolar, p_test);
            const double UL = fast_sat_lookup(iUmolar, p_test);
            const double UV = fast_sat_lookup_vapor(iUmolar, p_test);
            const double vL = 1.0 / DL, vV = 1.0 / DV;
            const double dv = vV - vL, du = UV - UL;
            if (std::abs(dv) < 1e-20 || std::abs(du) < 1e-20) return _HUGE;
            return (v_target - vL) / dv - (U_target - UL) / du;
        };

        double f_lo = fast_f(p_lo_sat);
        double f_hi = fast_f(p_hi_sat);
        // Use geometric mean (log-midpoint) as initial guess — much better than
        // arithmetic mean for the ~36000× pressure range of water.
        double p_sat = std::sqrt(p_lo_sat * p_hi_sat);

        if (f_lo * f_hi < 0.0) {
            double p_bracket_lo = p_lo_sat, p_bracket_hi = p_hi_sat;
            double f_bracket_lo = f_lo;

            for (int iter = 0; iter < 60; ++iter) {
                double f0 = fast_f(p_sat);
                if (std::abs(f0) < 1e-10) break;
                double dp = p_sat * 1e-5;
                double fp = fast_f(p_sat + dp);
                double dfdp = (fp - f0) / dp;
                double p_new = (std::abs(dfdp) > 1e-30) ? p_sat - f0 / dfdp : _HUGE;
                if (p_new < p_bracket_lo || p_new > p_bracket_hi || p_new != p_new) {
                    if ((f0 > 0) == (f_bracket_lo > 0)) { p_bracket_lo = p_sat; f_bracket_lo = f0; }
                    else p_bracket_hi = p_sat;
                    p_sat = 0.5 * (p_bracket_lo + p_bracket_hi);
                } else {
                    p_sat = p_new;
                }
            }

            const double DL = fast_sat_lookup(iDmolar, p_sat);
            const double DV = fast_sat_lookup_vapor(iDmolar, p_sat);
            const double UL = fast_sat_lookup(iUmolar, p_sat);
            const double UV = fast_sat_lookup_vapor(iUmolar, p_sat);
            const double vL = 1.0 / DL, vV = 1.0 / DV;
            const double Q_v = (std::abs(vV - vL) > 1e-20) ? (v_target - vL) / (vV - vL) : -1.0;
            const double Q_u = (std::abs(UV - UL) > 1e-20) ? (U_target - UL) / (UV - UL) : -1.0;

            if (Q_v >= 0.01 && Q_v <= 0.99 && std::abs(Q_v - Q_u) < 0.05) {
                is_confirmed_two_phase = true;
                p_sat_found = p_sat;
            }
        }
    }

    if (is_confirmed_two_phase) {
        const double DL = fast_sat_lookup(iDmolar, p_sat_found);
        const double DV = fast_sat_lookup_vapor(iDmolar, p_sat_found);
        const double vL = 1.0 / DL, vV = 1.0 / DV;
        const double v_target = 1.0 / D_target;
        const double Q_val = std::max(0.0, std::min(1.0, (v_target - vL) / (vV - vL)));

        _p = p_sat_found;
        _Q = Q_val;
        _T = fast_sat_lookup(iT, p_sat_found);  // T_sat is independent of quality
        using_single_phase_table = false;
        _phase = iphase_twophase;

        if (!is_mixture) {
            // Cache saturation bracket indices from the original saturation table
            // for subsequent property evaluations (h, s, cp, etc.)
            PureFluidSaturationTableData& sat = dataset->pure_saturation;
            CoolPropDbl yL, yV;
            double T_dummy = 0.5 * (sat.TL.front() + sat.TL.back());
            sat.is_inside(iP, p_sat_found, iT, T_dummy, cached_saturation_iL, cached_saturation_iV, yL, yV);
        }
        return;
    }

    // ----- Step 3: Single-phase via DU-table forward evaluation -----
    // Select liquid (D > D_crit) or gas (D ≤ D_crit) table.
    // Use explicit pointer assignment because DULiquidTable and DUGasTable are
    // different derived types and the ternary operator cannot deduce the base reference.
    bool use_liquid_table = (D_target > D_crit);
    SinglePhaseGriddedTableData* du_tbl_ptr;
    std::vector<std::vector<CellCoeffs>>* du_coeff_ptr;
    if (use_liquid_table) {
        du_tbl_ptr   = &du_liquid;
        du_coeff_ptr = &coeffs_du_liquid;
    } else {
        du_tbl_ptr   = &du_gas;
        du_coeff_ptr = &coeffs_du_gas;
    }
    SinglePhaseGriddedTableData& du_tbl  = *du_tbl_ptr;
    std::vector<std::vector<CellCoeffs>>& du_coeff = *du_coeff_ptr;

    std::size_t du_i = 0, du_j = 0;
    find_native_nearest_good_indices(du_tbl, du_coeff, D_target, U_target, du_i, du_j);

    // Direct forward evaluation: T(D, U) and P(D, U)
    double T_found = evaluate_single_phase(du_tbl, du_coeff, iT, D_target, U_target, du_i, du_j);
    double P_found = evaluate_single_phase(du_tbl, du_coeff, iP, D_target, U_target, du_i, du_j);
    // NOTE: P accuracy is limited for nearly-incompressible liquid states at low pressure.
    // For water near saturation at P < ~1 MPa, ΔD/D < 0.01% per bar of pressure change,
    // so adjacent (D, U) states at different P are indistinguishable at 200-point table
    // resolution.  T accuracy is always good (≤0.001%).  The P error is physically
    // inconsequential because downstream property differences (Δh, Δs) are negligible for
    // incompressible fluid when ΔP << bulk modulus.

    _T = T_found;
    _p = P_found;

    // Map (T, P) to the pT table for subsequent property lookups (h, s, cp, etc.)
    auto& pT_tbl   = dataset->single_phase_logpT;
    auto& pT_coeff = dataset->coeffs_pT;
    std::size_t pT_i = 0, pT_j = 0;
    find_native_nearest_good_indices(pT_tbl, pT_coeff, _T, _p, pT_i, pT_j);

    using_single_phase_table = true;
    selected_table = SELECTED_PT_TABLE;
    cached_single_phase_i = pT_i;
    cached_single_phase_j = pT_j;
    recalculate_singlephase_phase();
}

}  // namespace CoolProp

#endif  // !defined(NO_TABULAR_BACKENDS)
