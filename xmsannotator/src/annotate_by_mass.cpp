#include <limits>
#include <vector>
#include <tuple>

#include <Rcpp.h>

// [[Rcpp::plugins("cpp11")]]


using namespace Rcpp;


static bool almost_equal(double a, double b, double tolerance) noexcept {
    return std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;
}


template <unsigned Index, typename Pack>
static NumericVector unpack(const std::vector<Pack>& packed) {
    return NumericVector::import_transform(std::begin(packed), std::end(packed), [&] (const Pack& pack) {
        return std::get<Index>(pack);
    });
}


// [[Rcpp::export]]
DataFrame annotate_by_mass(DataFrame peaks, DataFrame adducts, DataFrame compounds, double tolerance) {
    const auto peak_mass = as<NumericVector>(peaks["mz"]);
    const auto adduct_mass = as<NumericVector>(adducts["mass"]);
    const auto adduct_mode = as<NumericVector>(adducts["mode"]);
    const auto adduct_multiplicity = as<NumericVector>(adducts["multiplicity"]);
    const auto compound_mass = as<NumericVector>(compounds["monoisotopic_mass"]);

    const auto n_peaks = peak_mass.size();
    const auto n_adducts = adduct_mass.size();
    const auto n_compounds = compound_mass.size();

    std::vector<std::tuple<R_xlen_t, R_xlen_t, R_xlen_t, double>> annotations;

    for (R_xlen_t adduct = 0; adduct < n_adducts; ++adduct) {
        for (R_xlen_t compound = 0; compound < n_compounds; ++compound) {
            const auto expected_mass = (adduct_multiplicity[adduct] * compound_mass[compound] + adduct_mass[adduct]) / adduct_mode[adduct];

            for (R_xlen_t peak = 0; peak < n_peaks; ++peak) {
                if (almost_equal(peak_mass[peak], expected_mass, tolerance))
                    annotations.emplace_back(peak, adduct, compound, expected_mass);
            }
        }
    }

    const auto peak_name = as<NumericVector>(peaks["peak"]);
    const auto adduct_name = as<CharacterVector>(adducts["adduct"]);
    const auto compound_name = as<NumericVector>(compounds["compound"]);

    return DataFrame::create(
        Named("peak") = peak_name[unpack<0>(annotations)],
        Named("adduct") = adduct_name[unpack<1>(annotations)],
        Named("compound") = compound_name[unpack<2>(annotations)],
        Named("expected_mass") = unpack<3>(annotations)
    );
}
