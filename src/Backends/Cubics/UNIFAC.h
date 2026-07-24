#ifndef UNIFAC_H_
#define UNIFAC_H_

#include <map>
#include <set>
#include <vector>

#include "UNIFACLibrary.h"
#include "CoolProp/Exceptions.h"

/// Structure containing data for the pure fluid in the mixture.
/// X/theta/lnGamma are dense vectors indexed by the mixture's compact group index
/// (see UNIFACMixture::m_groups); entries for groups absent from this component are zero.
struct ComponentData
{
    std::vector<double> X, theta, lnGamma;
    int group_count;  ///< The total number of groups in the pure fluid
};

namespace UNIFAC {
class UNIFACMixture
{
   private:
    /// A const reference to the library of group and interaction parameters
    const UNIFACLibrary::UNIFACParameterLibrary& library;

    std::size_t N = 0;  ///< Number of components

    double m_T = _HUGE;  ///< The temperature in K
    double T_r;          ///< Reducing temperature

    // Group-indexed quantities are stored in DENSE arrays indexed by a compact group index
    // (0..m_G-1), built once in set_pure_data from the sorted set of unique subgroups.  This
    // replaces the std::map<...> lookups that dominated the hot UNIFAC evaluation loops.  The
    // compact index preserves the ascending-sgi ordering of the original std::map iteration, so
    // summation order -- and therefore results -- are bit-for-bit unchanged.
    std::size_t m_G = 0;                        ///< Number of unique subgroups in the mixture
    std::vector<std::size_t> m_groups;          ///< Compact index -> subgroup index (sgi), ascending
    std::map<std::size_t, std::size_t> m_gidx;  ///< sgi -> compact group index (cold-path lookups only)

    std::vector<double> Psi_;  ///< Dense m_G x m_G Psi table, row-major: Psi_[gi*m_G + gj]

    std::vector<double> m_Xg,  ///< Compact-group-indexed mole fraction of group in the mixture
      m_thetag,                ///< Compact-group-indexed theta for the group in the mixture
      m_lnGammag,              ///< Compact-group-indexed ln(Gamma) for the group in the mixture
      m_Q;                     ///< Compact-group-indexed Q for the group

    std::vector<std::vector<int>> m_group_count;  ///< [component][compact group index] -> group count

    /// Temperature-keyed cache of the fully-computed group tables (dense Psi, per-component pure
    /// reference ln(Gamma), and the mixture group residual ln(Gamma)) at the current composition.
    /// The finite-difference tau-derivatives in ln_gamma_R evaluate at perturbed temperatures
    /// (tau +/- dtau) that repeat across the per-component loop and across density iterations at
    /// fixed (T, composition); keying by T (rather than a single slot) lets those repeats hit the
    /// cache instead of rebuilding the tables.  All entries depend on the composition, so the cache
    /// is cleared whenever the composition changes (set_mole_fractions) -- which, because the flash
    /// Jacobian is analytical (not composition-perturbed), happens only per trial phase.
    struct GroupTables
    {
        std::vector<double> Psi;           ///< dense m_G x m_G
        std::vector<double> pure_lnGamma;  ///< N x m_G, flattened [i*m_G + g]
        std::vector<double> lnGammag;      ///< m_G
    };
    std::map<double, GroupTables> m_T_cache;           ///< T -> tables (at the current composition)
    bool m_tables_valid = false;                       ///< do the working tables hold the current (m_T, composition)?
    static constexpr std::size_t m_T_cache_cap = 256;  ///< safety bound on distinct cached temperatures

    /// A map from (i, j) indices for subgroup, subgroup indices to the interaction parameters for this pair
    std::map<std::pair<int, int>, UNIFACLibrary::InteractionParameters> interaction;

    /// A map from SGI to MGI
    std::map<std::size_t, std::size_t> m_sgi_to_mgi;

    /// The set of unique groups in this mixture
    std::set<std::size_t> unique_groups;

    std::vector<double> mole_fractions;

    std::vector<UNIFACLibrary::Component> components;

    std::vector<ComponentData> pure_data;

   public:
    UNIFACMixture(const UNIFACLibrary::UNIFACParameterLibrary& library, const double T_r) : library(library), T_r(T_r){};

    /** \brief Set all the interaction parameters between groups */
    void set_interaction_parameters();
    void set_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter, const double value);
    /// Get one of the mgi-mgi interaction pairs
    double get_interaction_parameter(const std::size_t mgi1, const std::size_t mgi2, const std::string& parameter);

    /// Set the mole fractions of the components in the mixtures (not the groups)
    void set_mole_fractions(const std::vector<double>& z);

    /// Get the mole fractions of the components in the mixtures (not the groups)
    const std::vector<double>& get_mole_fractions() {
        return mole_fractions;
    }

    /// Set the temperature of the components in the mixtures (not the groups)
    void set_temperature(const double T);

    /// Get the temperature
    double get_temperature() const {
        return m_T;
    }

    double Psi(std::size_t sgi1, std::size_t sgi2) const;

    double theta_pure(std::size_t i, std::size_t sgi) const;

    void activity_coefficients(double tau, const std::vector<double>& z, std::vector<double>& gamma);

    double ln_gamma_R(const double tau, std::size_t i, std::size_t itau);

    std::size_t group_count(std::size_t i, std::size_t sgi) const;

    /// Add a component with the defined groups defined by (count, sgi) pairs
    void add_component(const UNIFACLibrary::Component& comp);

    void set_components(const std::string& identifier_type, const std::vector<std::string>& identifiers);

    const std::vector<UNIFACLibrary::Component>& get_components() {
        return components;
    };

    void set_pure_data();

    /// Modify the surface parameter Q_k of the sub group sgi
    void set_Q_k(const size_t sgi, const double value);

    /// Get the surface parameter Q_k of the sub group sgi
    double get_Q_k(const size_t sgi) const;
};

} /* namespace UNIFAC */

#endif
