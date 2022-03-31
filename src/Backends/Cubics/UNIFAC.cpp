#include "UNIFAC.h"

void UNIFAC::UNIFACMixture::set_interaction_parameters() {
    for (std::set<std::size_t>::const_iterator itisgi = unique_groups.begin(); itisgi != unique_groups.end(); ++itisgi) {
        for (std::set<std::size_t>::const_iterator itjsgi = unique_groups.begin(); itjsgi != unique_groups.end(); ++itjsgi) {
            if (*itjsgi >= *itisgi) {
                continue;
            }
            std::size_t mgi1 = m_sgi_to_mgi.find(*itisgi)->second;
            std::size_t mgi2 = m_sgi_to_mgi.find(*itjsgi)->second;
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
    std::map<std::pair<int, int>, UNIFACLibrary::InteractionParameters>::iterator it = this->interaction.find(std::pair<int, int>(mgi1, mgi2));
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

    std::map<std::size_t, double>&Xg = m_Xg, &thetag = m_thetag;
    Xg.clear();
    thetag.clear();

    // Iterate over the fluids
    double X_summer = 0;
    for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
        X_summer += this->mole_fractions[i] * pure_data[i].group_count;
    }
    /// Calculations for each group in the total mixture
    for (std::set<std::size_t>::iterator itsgi = unique_groups.begin(); itsgi != unique_groups.end(); ++itsgi) {
        double X = 0;
        // Iterate over the fluids
        for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
            X += this->mole_fractions[i] * group_count(i, *itsgi);
        }
        Xg.insert(std::pair<std::size_t, double>(*itsgi, X));
    }
    /// Now come back through and divide by the sum(z_i*count) for this fluid
    for (std::map<std::size_t, double>::iterator it = Xg.begin(); it != Xg.end(); ++it) {
        it->second /= X_summer;
        //printf("X_{%d}: %g\n", it->first, it->second);
    }
    double theta_summer = 0;
    for (std::set<std::size_t>::iterator itsgi = unique_groups.begin(); itsgi != unique_groups.end(); ++itsgi) {
        double cont = Xg.find(*itsgi)->second * m_Q.find(*itsgi)->second;
        theta_summer += cont;
        thetag.insert(std::pair<std::size_t, double>(*itsgi, cont));
    }
    /// Now come back through and divide by the sum(X*Q) for this fluid
    for (std::map<std::size_t, double>::iterator it = thetag.begin(); it != thetag.end(); ++it) {
        it->second /= theta_summer;
        //printf("theta_{%d}: %g\n", it->first, it->second);
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
        std::map<std::pair<int, int>, UNIFACLibrary::InteractionParameters>::const_iterator it =
          this->interaction.find(std::pair<int, int>(static_cast<int>(mgi1), static_cast<int>(mgi2)));
        if (it != this->interaction.end()) {
            return exp(-(it->second.a_ij / this->m_T + it->second.b_ij + it->second.c_ij * this->m_T));
        } else {
            throw CoolProp::ValueError(
              format("Could not match mgi[%d]-mgi[%d] interaction in UNIFAC", static_cast<int>(mgi1), static_cast<int>(mgi2)));
        }
    }
}

std::size_t UNIFAC::UNIFACMixture::group_count(std::size_t i, std::size_t sgi) const {
    const UNIFACLibrary::Component& c = components[i];
    for (std::vector<UNIFACLibrary::ComponentGroup>::const_iterator it = c.groups.begin(); it != c.groups.end(); ++it) {
        if (it->group.sgi == sgi) {
            return it->count;
        }
    }
    return 0;
}

double UNIFAC::UNIFACMixture::theta_pure(std::size_t i, std::size_t sgi) const {
    return pure_data[i].theta.find(sgi)->second;
}

void UNIFAC::UNIFACMixture::set_temperature(const double T) {
    //    // Check whether you are using exactly the same temperature as last time
    //    if (static_cast<bool>(_T) && std::abs(static_cast<double>(_T) - T) < 1e-15 {
    //        //
    //        return;
    //    }
    this->m_T = T;

    if (this->mole_fractions.empty()) {
        throw CoolProp::ValueError("mole fractions must be set before calling set_temperature");
    }

    // Compute Psi once for the different calls
    for (std::set<std::size_t>::iterator itk = unique_groups.begin(); itk != unique_groups.end(); ++itk) {
        for (std::set<std::size_t>::iterator itm = unique_groups.begin(); itm != unique_groups.end(); ++itm) {
            Psi_[std::pair<std::size_t, std::size_t>(*itk, *itm)] = Psi(*itk, *itm);
        }
    }

    for (std::size_t i = 0; i < this->mole_fractions.size(); ++i) {
        const UNIFACLibrary::Component& c = components[i];
        for (std::size_t k = 0; k < c.groups.size(); ++k) {
            double Q = c.groups[k].group.Q_k;
            int sgik = c.groups[k].group.sgi;
            double sum1 = 0;
            for (std::size_t m = 0; m < c.groups.size(); ++m) {
                int sgim = c.groups[m].group.sgi;
                sum1 += theta_pure(i, sgim) * Psi_.find(std::pair<std::size_t, std::size_t>(sgim, sgik))->second;
            }
            double s = 1 - log(sum1);
            for (std::size_t m = 0; m < c.groups.size(); ++m) {
                int sgim = c.groups[m].group.sgi;
                double sum2 = 0;
                for (std::size_t n = 0; n < c.groups.size(); ++n) {
                    int sgin = c.groups[n].group.sgi;
                    sum2 += theta_pure(i, sgin) * Psi_.find(std::pair<std::size_t, std::size_t>(sgin, sgim))->second;
                }
                s -= theta_pure(i, sgim) * Psi_.find(std::pair<std::size_t, std::size_t>(sgik, sgim))->second / sum2;
            }
            pure_data[i].lnGamma[sgik] = Q * s;
            //printf("ln(Gamma)^(%d)_{%d}: %g\n", static_cast<int>(i + 1), sgik, Q*s);
        }
    }

    std::map<std::size_t, double>&thetag = m_thetag, &lnGammag = m_lnGammag;
    lnGammag.clear();

    for (std::set<std::size_t>::iterator itksgi = unique_groups.begin(); itksgi != unique_groups.end(); ++itksgi) {
        double sum1 = 0;
        for (std::set<std::size_t>::iterator itmsgi = unique_groups.begin(); itmsgi != unique_groups.end(); ++itmsgi) {
            sum1 += thetag.find(*itmsgi)->second * Psi_.find(std::pair<std::size_t, std::size_t>(*itmsgi, *itksgi))->second;
        }
        double s = 1 - log(sum1);
        for (std::set<std::size_t>::iterator itmsgi = unique_groups.begin(); itmsgi != unique_groups.end(); ++itmsgi) {
            double sum3 = 0;
            for (std::set<std::size_t>::iterator itnsgi = unique_groups.begin(); itnsgi != unique_groups.end(); ++itnsgi) {
                sum3 += thetag.find(*itnsgi)->second * Psi_.find(std::pair<std::size_t, std::size_t>(*itnsgi, *itmsgi))->second;
            }
            s -= thetag.find(*itmsgi)->second * Psi_.find(std::pair<std::size_t, std::size_t>(*itksgi, *itmsgi))->second / sum3;
        }
        lnGammag.insert(std::pair<std::size_t, double>(*itksgi, m_Q.find(*itksgi)->second * s));
        //printf("log(Gamma)_{%d}: %g\n", itk->sgi, itk->Q_k*s);
    }
    _T = m_T;
}
double UNIFAC::UNIFACMixture::ln_gamma_R(const double tau, std::size_t i, std::size_t itau) {
    if (itau == 0) {
        set_temperature(T_r / tau);
        double summer = 0;
        for (std::set<std::size_t>::const_iterator itsgi = unique_groups.begin(); itsgi != unique_groups.end(); ++itsgi) {
            std::size_t count = group_count(i, *itsgi);
            if (count > 0) {
                summer += count * (m_lnGammag.find(*itsgi)->second - pure_data[i].lnGamma.find(*itsgi)->second);
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
        for (std::size_t j = 0; j < c.groups.size(); ++j) {
            const UNIFACLibrary::ComponentGroup& cg = c.groups[j];
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
    for (std::vector<UNIFACLibrary::ComponentGroup>::const_iterator it = comp.groups.begin(); it != comp.groups.end(); ++it) {
        m_sgi_to_mgi.insert(std::pair<std::size_t, std::size_t>(it->group.sgi, it->group.mgi));
    }
}

void UNIFAC::UNIFACMixture::set_components(const std::string& identifier_type, std::vector<std::string> identifiers) {
    components.clear();
    N = identifiers.size();
    if (identifier_type == "name") {
        // Iterate over the provided names
        for (std::vector<std::string>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it) {
            // Get and add the component
            UNIFACLibrary::Component c = library.get_component("name", *it);
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
    pure_data.clear();
    unique_groups.clear();
    m_Q.clear();
    for (std::size_t i = 0; i < N; ++i) {
        const UNIFACLibrary::Component& c = components[i];
        ComponentData cd;
        double summerxq = 0;
        cd.group_count = 0;
        for (std::size_t j = 0; j < c.groups.size(); ++j) {
            const UNIFACLibrary::ComponentGroup& cg = c.groups[j];
            double x = static_cast<double>(cg.count);
            double theta = static_cast<double>(cg.count * cg.group.Q_k);
            cd.X.insert(std::pair<int, double>(cg.group.sgi, x));
            cd.theta.insert(std::pair<int, double>(cg.group.sgi, theta));
            cd.group_count += cg.count;
            summerxq += x * cg.group.Q_k;
            unique_groups.insert(cg.group.sgi);
            m_Q.insert(std::pair<std::size_t, double>(cg.group.sgi, cg.group.Q_k));
        }
        /// Now come back through and divide by the total # groups for this fluid
        for (std::map<std::size_t, double>::iterator it = cd.X.begin(); it != cd.X.end(); ++it) {
            it->second /= cd.group_count;
            //printf("X^(%d)_{%d}: %g\n", static_cast<int>(i + 1), static_cast<int>(it->first), it->second);
        }
        /// Now come back through and divide by the sum(X*Q) for this fluid
        for (std::map<std::size_t, double>::iterator it = cd.theta.begin(); it != cd.theta.end(); ++it) {
            it->second /= summerxq;
            //printf("theta^(%d)_{%d}: %g\n", static_cast<int>(i+1), static_cast<int>(it->first), it->second);
        }
        pure_data.push_back(cd);
    }
}

/// Modify the surface parameter Q_k of the sub group sgi
void UNIFAC::UNIFACMixture::set_Q_k(const size_t sgi, const double value) {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < components[i].groups.size(); ++j) {
            if (components[i].groups[j].group.sgi == sgi) {
                components[i].groups[j].group.Q_k = value;
            }
        }
    }

    /// Re-calculate the parameters X and theta for the pure components, which does not depend on temperature nor molar fraction
    set_pure_data();
}

/// Modify the surface parameter Q_k of the sub group sgi
double UNIFAC::UNIFACMixture::get_Q_k(const size_t sgi) const {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < components[i].groups.size(); ++j) {
            if (components[i].groups[j].group.sgi == sgi) {
                return components[i].groups[j].group.Q_k;
            }
        }
    }
    throw CoolProp::ValueError("Could not get Q_k");
}
