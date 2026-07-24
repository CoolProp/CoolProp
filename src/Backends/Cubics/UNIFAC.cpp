#include "UNIFAC.h"

#include <algorithm>

void UNIFAC::UNIFACMixture::set_interaction_parameters() {
    for (const auto& isgi : unique_groups) {
        for (const auto& jsgi : unique_groups) {
            if (jsgi >= isgi) {
                continue;
            }
            int mgi1 = static_cast<int>(m_sgi_to_mgi.find(isgi)->second);
            int mgi2 = static_cast<int>(m_sgi_to_mgi.find(jsgi)->second);
            // Insert in normal order
            std::pair<std::pair<int, int>, UNIFACLibrary::InteractionParameters> m_pair(std::pair<int, int>(mgi1, mgi2),
                                                                                        library.get_interaction_parameters(mgi1, mgi2));
            interaction.insert(m_pair);
            // Insert in backwards order
            if (mgi1 != mgi2) {
                std::pair<std::pair<int, int>, UNIFACLibrary::InteractionParameters> m_pair(std::pair<int, int>(mgi2, mgi1),
                                                                                            library.get_interaction_parameters(mgi2, mgi1));
                interaction.insert(m_pair);
            }
        }
    }
}

void UNIFAC::UNIFACMixture::set_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter,
                                                      const double value) {
    if (parameter == "aij") {
        this->interaction[std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2))].a_ij = value;
    } else if (parameter == "bij") {
        this->interaction[std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2))].b_ij = value;
    } else if (parameter == "cij") {
        this->interaction[std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2))].c_ij = value;
    } else {
        throw CoolProp::ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
    }
}
double UNIFAC::UNIFACMixture::get_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter) {
    auto it = this->interaction.find(std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2)));
    if (it == this->interaction.end()) {
        throw CoolProp::ValueError(format("Unable to match mgi-mgi pair: [%d,%d]", static_cast<int>(mgi1), static_cast<int>(mgi1)));
    } else {
        if (parameter == "aij") {
            return it->second.a_ij;
        } else if (parameter == "bij") {
            return it->second.b_ij;
        } else if (parameter == "cij") {
            return it->second.c_ij;
        } else {
            throw CoolProp::ValueError(format("I don't know what to do with parameter [%s]", parameter.c_str()));
        }
    }
}

/// Set the mole fractions of the components in the mixtures (not the groups)
void UNIFAC::UNIFACMixture::set_mole_fractions(const std::vector<double>& z) {
    //    // If the vector fractions are the same as last ones, don't do anything and return
    //    if (!mole_fractions.empty() && maxvectordiff(z, mole_fractions) < 1e-15){
    //        return;
    //    }
    this->mole_fractions = z;
    if (this->N != z.size()) {
        throw CoolProp::ValueError("Size of molar fraction do not match number of components.");
    }
    // Composition changed -> the temperature-keyed group-table cache (whose mixture group residual
    // depends on the group surface fractions m_thetag computed below) is now stale.
    m_T_cache.clear();
    m_tables_valid = false;

    m_Xg.assign(m_G, 0.0);
    m_thetag.assign(m_G, 0.0);

    // Iterate over the fluids
    double X_summer = 0;
    for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
        X_summer += this->mole_fractions[i] * pure_data[i].group_count;
    }
    /// Calculations for each group in the total mixture (ascending compact-group order == ascending sgi)
    for (std::size_t g = 0; g < m_G; ++g) {
        double X = 0;
        // Iterate over the fluids
        for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
            X += this->mole_fractions[i] * m_group_count[i][g];
        }
        m_Xg[g] = X;
    }
    /// Now come back through and divide by the sum(z_i*count) for this fluid
    for (std::size_t g = 0; g < m_G; ++g) {
        m_Xg[g] /= X_summer;
        //printf("X_{%d}: %g\n", m_groups[g], m_Xg[g]);
    }
    double theta_summer = 0;
    for (std::size_t g = 0; g < m_G; ++g) {
        double cont = m_Xg[g] * m_Q[g];
        theta_summer += cont;
        m_thetag[g] = cont;
    }
    /// Now come back through and divide by the sum(X*Q) for this fluid
    for (std::size_t g = 0; g < m_G; ++g) {
        m_thetag[g] /= theta_summer;
        //printf("theta_{%d}: %g\n", m_groups[g], m_thetag[g]);
    }
}

double UNIFAC::UNIFACMixture::Psi(std::size_t sgi1, std::size_t sgi2) const {

    if (this->interaction.size() == 0) {
        throw CoolProp::ValueError("interaction parameters for UNIFAC not yet set");
    }
    std::size_t mgi1 = m_sgi_to_mgi.find(sgi1)->second;
    std::size_t mgi2 = m_sgi_to_mgi.find(sgi2)->second;
    if (mgi1 == mgi2) {
        return 1;
    } else {
        auto it = this->interaction.find(std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2)));
        if (it != this->interaction.end()) {
            return exp(-(it->second.a_ij / this->m_T + it->second.b_ij + it->second.c_ij * this->m_T));
        } else {
            throw CoolProp::ValueError(
              format("Could not match mgi[%d]-mgi[%d] interaction in UNIFAC", static_cast<int>(mgi1), static_cast<int>(mgi2)));
        }
    }
}

std::size_t UNIFAC::UNIFACMixture::group_count(std::size_t i, std::size_t sgi) const {
    auto it = m_gidx.find(sgi);
    if (it == m_gidx.end()) {
        return 0;
    }
    return static_cast<std::size_t>(m_group_count[i][it->second]);
}

double UNIFAC::UNIFACMixture::theta_pure(std::size_t i, std::size_t sgi) const {
    return pure_data[i].theta[m_gidx.find(sgi)->second];
}

void UNIFAC::UNIFACMixture::set_temperature(const double T) {
    if (this->mole_fractions.empty()) {
        throw CoolProp::ValueError("mole fractions must be set before calling set_temperature");
    }

    // Fast path: the working tables already hold this T at the current composition.  The exact
    // floating-point comparison is intentional: T is a cache key, not a physical tolerance -- the
    // caller re-passes the identical T value (T_r/tau at fixed tau) across a flash, and the perturbed
    // temperatures of the finite-difference tau-derivatives are likewise exact repeats.  A tolerance
    // would wrongly return tables computed at a *different* temperature.  (Same rationale for the
    // std::map<double, ...> lookup below.)
    if (m_tables_valid && m_T == T) {
        return;
    }
    m_T = T;

    // Cache hit: restore the fully-computed tables for this T (at the current composition).  The
    // finite-difference tau-derivatives in ln_gamma_R evaluate at perturbed temperatures that repeat
    // across the per-component loop and across density iterations, so this restore (a few small
    // vector copies) replaces a full Psi + pure + group-residual rebuild.
    auto it = m_T_cache.find(T);
    if (it != m_T_cache.end()) {
        const GroupTables& t = it->second;
        Psi_ = t.Psi;
        for (std::size_t i = 0; i < N; ++i) {
            std::copy(t.pure_lnGamma.begin() + i * m_G, t.pure_lnGamma.begin() + (i + 1) * m_G, pure_data[i].lnGamma.begin());
        }
        m_lnGammag = t.lnGammag;
        m_tables_valid = true;
        return;
    }

    // Miss: compute the dense Psi table (row-major, [gk*m_G + gm]), the pure-component reference
    // ln(Gamma), and the mixture group residual ln(Gamma), then store them keyed by T.
    Psi_.assign(m_G * m_G, 0.0);
    for (std::size_t gk = 0; gk < m_G; ++gk) {
        for (std::size_t gm = 0; gm < m_G; ++gm) {
            Psi_[gk * m_G + gm] = Psi(m_groups[gk], m_groups[gm]);
        }
    }

    for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
        const UNIFACLibrary::Component& c = components[i];
        const std::vector<double>& theta_i = pure_data[i].theta;
        for (std::size_t k = 0; k < c.groups.size(); ++k) {
            double Q = c.groups[k].group.Q_k;
            std::size_t gk = m_gidx[c.groups[k].group.sgi];
            double sum1 = 0;
            for (const auto& group : c.groups) {
                std::size_t gm = m_gidx[group.group.sgi];
                sum1 += theta_i[gm] * Psi_[gm * m_G + gk];
            }
            // sum1 > 0: c.groups is the set of groups present in component i, theta_pure
            // returns the (positive) surface-area fraction of each, and Psi = exp(-a/T) > 0.
            // cppcheck cannot prove this symbolically.
            // cppcheck-suppress invalidFunctionArg
            double s = 1 - log(sum1);
            for (std::size_t m = 0; m < c.groups.size(); ++m) {
                std::size_t gm = m_gidx[c.groups[m].group.sgi];
                double sum2 = 0;
                for (const auto& group : c.groups) {
                    std::size_t gn = m_gidx[group.group.sgi];
                    sum2 += theta_i[gn] * Psi_[gn * m_G + gm];
                }
                s -= theta_i[gm] * Psi_[gk * m_G + gm] / sum2;
            }
            pure_data[i].lnGamma[gk] = Q * s;
            //printf("ln(Gamma)^(%d)_{%d}: %g\n", static_cast<int>(i + 1), m_groups[gk], Q*s);
        }
    }

    m_lnGammag.assign(m_G, 0.0);
    for (std::size_t gk = 0; gk < m_G; ++gk) {
        double sum1 = 0;
        for (std::size_t gm = 0; gm < m_G; ++gm) {
            sum1 += m_thetag[gm] * Psi_[gm * m_G + gk];
        }
        // sum1 > 0: m_thetag are the (positive) group surface fractions and Psi = exp(-a/T) > 0,
        // so the sum over the (non-empty) group set is strictly positive.  cppcheck cannot prove
        // this symbolically.
        // cppcheck-suppress invalidFunctionArg
        double s = 1 - log(sum1);
        for (std::size_t gm = 0; gm < m_G; ++gm) {
            double sum3 = 0;
            for (std::size_t gn = 0; gn < m_G; ++gn) {
                sum3 += m_thetag[gn] * Psi_[gn * m_G + gm];
            }
            s -= m_thetag[gm] * Psi_[gk * m_G + gm] / sum3;
        }
        m_lnGammag[gk] = m_Q[gk] * s;
        //printf("log(Gamma)_{%d}: %g\n", m_groups[gk], m_Q[gk]*s);
    }
    m_tables_valid = true;

    // Store a snapshot for this T (at the current composition).  The cap bounds memory if a caller
    // sweeps many distinct temperatures at fixed composition (e.g. a temperature-search flash).
    if (m_T_cache.size() >= m_T_cache_cap) {
        m_T_cache.clear();
    }
    GroupTables t;
    t.Psi = Psi_;
    t.pure_lnGamma.resize(N * m_G);
    for (std::size_t i = 0; i < N; ++i) {
        std::copy(pure_data[i].lnGamma.begin(), pure_data[i].lnGamma.end(), t.pure_lnGamma.begin() + i * m_G);
    }
    t.lnGammag = m_lnGammag;
    m_T_cache.emplace(T, std::move(t));
}
double UNIFAC::UNIFACMixture::ln_gamma_R(const double tau, std::size_t i, std::size_t itau) {
    if (itau == 0) {
        set_temperature(T_r / tau);
        double summer = 0;
        const std::vector<int>& gc_i = m_group_count[i];
        const std::vector<double>& lnGamma_i = pure_data[i].lnGamma;
        for (std::size_t g = 0; g < m_G; ++g) {
            int count = gc_i[g];
            if (count > 0) {
                summer += count * (m_lnGammag[g] - lnGamma_i[g]);
            }
        }
        //printf("log(gamma)_{%d}: %g\n", i+1, summer);
        return summer;
    } else {
        double dtau = 0.01 * tau;
        return (ln_gamma_R(tau + dtau, i, itau - 1) - ln_gamma_R(tau - dtau, i, itau - 1)) / (2 * dtau);
    }
}
void UNIFAC::UNIFACMixture::activity_coefficients(double tau, const std::vector<double>& z, std::vector<double>& gamma) {
    if (this->N != z.size()) {
        throw CoolProp::ValueError("Size of molar fraction do not match number of components.");
    }
    std::vector<double> r(N), q(N), l(N), phi(N), theta(N), ln_Gamma_C(N);
    double summerzr = 0, summerzq = 0, summerzl = 0;
    for (std::size_t i = 0; i < N; ++i) {
        double summerr = 0, summerq = 0;
        const UNIFACLibrary::Component& c = components[i];
        for (const auto& cg : c.groups) {
            summerr += cg.count * cg.group.R_k;
            summerq += cg.count * cg.group.Q_k;
        }
        r[i] = summerr;
        q[i] = summerq;
        summerzr += z[i] * r[i];
        summerzq += z[i] * q[i];
    }
    for (std::size_t i = 0; i < N; ++i) {
        phi[i] = z[i] * r[i] / summerzr;
        theta[i] = z[i] * q[i] / summerzq;
        l[i] = 10.0 / 2.0 * (r[i] - q[i]) - (r[i] - 1);
        summerzl += z[i] * l[i];
    }
    for (std::size_t i = 0; i < N; ++i) {
        ln_Gamma_C[i] = log(phi[i] / z[i]) + 10.0 / 2.0 * q[i] * log(theta[i] / phi[i]) + l[i] - phi[i] / z[i] * summerzl;
        gamma[i] = exp(ln_gamma_R(tau, i, 0) + ln_Gamma_C[i]);
    }
}

/// Add a component with the defined groups defined by (count, sgi) pairs
void UNIFAC::UNIFACMixture::add_component(const UNIFACLibrary::Component& comp) {
    components.push_back(comp);
    for (const auto& cg : comp.groups) {
        m_sgi_to_mgi.emplace(cg.group.sgi, cg.group.mgi);
    }
}

void UNIFAC::UNIFACMixture::set_components(const std::string& identifier_type, const std::vector<std::string>& identifiers) {
    components.clear();
    N = identifiers.size();
    if (identifier_type == "name") {
        // Iterate over the provided names
        for (const auto& id : identifiers) {
            // Get and add the component
            UNIFACLibrary::Component c = library.get_component("name", id);
            add_component(c);
        }
    } else {
        throw CoolProp::ValueError("Cannot understand identifier_type");
    }
    /// Calculate the parameters X and theta for the pure components, which does not depend on temperature nor molar fraction
    set_pure_data();
}

/// Calculate the parameters X and theta for the pure components, which does not depend on temperature nor molar fraction
void UNIFAC::UNIFACMixture::set_pure_data() {
    // The group set / Q values are changing, so every cached temperature table is stale.
    m_T_cache.clear();
    m_tables_valid = false;

    pure_data.clear();
    unique_groups.clear();
    m_gidx.clear();
    m_groups.clear();

    // Pass 1: collect the unique subgroups across all components.  std::set keeps them ascending,
    // which fixes the compact-group ordering so that summation order matches the previous
    // std::map-based implementation (results stay bit-for-bit identical).
    for (std::size_t i = 0; i < N; ++i) {
        for (const auto& cg : components[i].groups) {
            unique_groups.insert(cg.group.sgi);
        }
    }
    m_G = unique_groups.size();
    m_groups.assign(unique_groups.begin(), unique_groups.end());
    for (std::size_t g = 0; g < m_G; ++g) {
        m_gidx[m_groups[g]] = g;
    }

    // Pass 2: build the dense per-group Q, per-component group counts, and per-component X / theta.
    m_Q.assign(m_G, 0.0);
    m_lnGammag.assign(m_G, 0.0);
    m_group_count.assign(N, std::vector<int>(m_G, 0));
    for (std::size_t i = 0; i < N; ++i) {
        const UNIFACLibrary::Component& c = components[i];
        ComponentData cd;
        cd.X.assign(m_G, 0.0);
        cd.theta.assign(m_G, 0.0);
        cd.lnGamma.assign(m_G, 0.0);
        cd.group_count = 0;
        double summerxq = 0;
        for (const auto& cg : c.groups) {
            std::size_t g = m_gidx[cg.group.sgi];
            auto x = static_cast<double>(cg.count);
            auto theta = static_cast<double>(cg.count * cg.group.Q_k);
            cd.X[g] = x;
            cd.theta[g] = theta;
            cd.group_count += cg.count;
            m_group_count[i][g] = cg.count;
            summerxq += x * cg.group.Q_k;
            m_Q[g] = cg.group.Q_k;  // Q_k is a property of the subgroup, identical across components
        }
        /// Now come back through and divide by the total # groups / the sum(X*Q) for this fluid.
        /// Entries for groups absent from this component stay zero (0 / c = 0) and are never read.
        for (std::size_t g = 0; g < m_G; ++g) {
            cd.X[g] /= cd.group_count;
            cd.theta[g] /= summerxq;
        }
        pure_data.push_back(cd);
    }

    // If a composition was already set (e.g. set_Q_k called after set_mole_fractions), the mixture
    // group surface fractions m_Xg/m_thetag depend on the group layout and Q that just changed, so
    // refresh them from the stored composition.  (When called from set_components the composition is
    // not yet set, so this is skipped.)
    if (!mole_fractions.empty() && mole_fractions.size() == N) {
        set_mole_fractions(std::vector<double>(mole_fractions));
    }
}

/// Modify the surface parameter Q_k of the sub group sgi
void UNIFAC::UNIFACMixture::set_Q_k(const size_t sgi, const double value) {
    for (std::size_t i = 0; i < N; ++i) {
        for (auto& group : components[i].groups) {
            if (group.group.sgi == sgi) {
                group.group.Q_k = value;
            }
        }
    }

    /// Re-calculate the parameters X and theta for the pure components, which does not depend on temperature nor molar fraction
    set_pure_data();
}

/// Modify the surface parameter Q_k of the sub group sgi
double UNIFAC::UNIFACMixture::get_Q_k(const size_t sgi) const {
    for (std::size_t i = 0; i < N; ++i) {
        for (const auto& group : components[i].groups) {
            if (group.group.sgi == sgi) {
                return group.group.Q_k;
            }
        }
    }
    throw CoolProp::ValueError("Could not get Q_k");
}
