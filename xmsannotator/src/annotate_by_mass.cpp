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
DataFrame annotate_by_mass(DataFrame peak_table, DataFrame adduct_table, DataFrame compound_table, double tolerance) {
    const auto peak_mz = as<NumericVector>(peak_table["mz"]);
    const auto adduct_mass = as<NumericVector>(adduct_table["mass"]);
    const auto adduct_charge = as<NumericVector>(adduct_table["charge"]);
    const auto adduct_molecules = as<NumericVector>(adduct_table["n_molecules"]);
    const auto compound_mass = as<NumericVector>(compound_table["monoisotopic_mass"]);

    const auto n_peaks = peak_table.nrows();
    const auto n_adducts = adduct_table.nrows();
    const auto n_compounds = compound_table.nrows();

    std::vector<std::tuple<R_xlen_t, R_xlen_t, R_xlen_t, double>> annotations;

    for (R_xlen_t adduct = 0; adduct < n_adducts; ++adduct) {
        for (R_xlen_t compound = 0; compound < n_compounds; ++compound) {
            const auto expected_mass = (adduct_molecules[adduct] * compound_mass[compound] + adduct_mass[adduct]) / adduct_charge[adduct];

            for (R_xlen_t peak = 0; peak < n_peaks; ++peak) {
                if (almost_equal(peak_mz[peak], expected_mass, tolerance))
                    annotations.emplace_back(peak, adduct, compound, expected_mass);
            }
        }
    }

    const auto peak = as<NumericVector>(peak_table["peak"]);
    const auto adduct = as<CharacterVector>(adduct_table["adduct"]);
    const auto compound = as<NumericVector>(compound_table["compound"]);

    return DataFrame::create(
        Named("peak") = peak[unpack<0>(annotations)],
        Named("adduct") = adduct[unpack<1>(annotations)],
        Named("compound") = compound[unpack<2>(annotations)],
        Named("expected_mass") = unpack<3>(annotations)
    );
}
