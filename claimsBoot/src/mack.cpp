#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "constants.h"
#include "helpers.h"

extern "C" {
void mack_boot(int n_dev, double *triangle, int boot_type, int opt, bool cond,
               int n_boot, double *reserve, int seed);

void mack_sim(int n_dev, double *triangle, int sim_type, int boot_type,
              int sim_dist, int n_boot, int n_conf, int m_conf, double *conf,
              double *res, bool show_progress, int seed);
}

//' Simulate Mack CL reserve.
//'
//' @param triangle Cumulative claims triangle
//' @param n_boot Number of bootstrap simulations
//' @param boot_type Type of bootstrap: `"parametric"`, `"residuals"`, or `"pairs"`
//' @param proc_dist Distribution of process error: `"normal"`or `"gamma"`
//' @param cond Specified whether the bootstrap should be conditional or unconditional. Default is `TRUE`.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector
mackBoot(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::String boot_type,
         Rcpp::String opt, bool cond, int seed = 42) {
  int n_dev = triangle.rows();
  triangle = na_to_zero(triangle);

  if (boot_type == "pairs" && !cond) {
    Rcpp::stop("Unconditional method invalid for pairs bootstrap.");
  }

  Rcpp::NumericVector reserve(n_boot);
  mack_boot(n_dev, triangle.begin(), key[boot_type], key[opt],
            cond, n_boot, reserve.begin(), seed);
  return (reserve);
};

//' Simulate Mack CL reserve for different perturbed and excluded points.
//'
//' @param triangle Cumulative claims triangle
//' @param sim_type Simulation type: `"single"` (the default), `"calender"`, or `"origin"`
//' @param n_boot Number of bootstrap simulations
//' @param factor Perturbation factor
//' @param boot_type Type of bootstrap: `"parametric"`, `"residuals"`, or `"pairs"`
//' @param proc_dist Distribution of process error: `"normal"` or `"gamma"`
//' @param cond Specified whether the bootstrap should be conditional or unconditional. Default is `TRUE`.
//'
//' @details The simulation configuration inputs `sim_type`, `factor`, `boot_type`, `proc_dist` and `conditional` can be either strings or character vectors. In the latter case, the simulation is computed for all feasible combinations.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame mackSim(Rcpp::NumericMatrix triangle, Rcpp::String sim_type,
                        int n_boot, Rcpp::NumericVector mean_factors,
                        Rcpp::NumericVector sd_factors,
                        Rcpp::CharacterVector boot_types,
                        Rcpp::String sim_dist = "normal",
                        bool show_progress = true, int seed = 42) {

  int n_dev = triangle.rows();
  triangle = na_to_zero(triangle);
  int n_pts = (pow(n_dev, 2) + n_dev) / 2;

  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("claimsBoot");
  Rcpp::Function mackConfig = pkg["mackConfig"];
  Rcpp::Function mackPost = pkg["mackPost"];

  Rcpp::List res_list(boot_types.length());
  for (Rcpp::CharacterVector::iterator i = boot_types.begin();
       i != boot_types.end(); ++i) {
    Rcpp::String boot_type = *i;
    Rcpp::NumericMatrix conf =
        mackConfig(n_dev, mean_factors, sd_factors, sim_type, boot_type);

    int sim_type_ = key[sim_type];
    int boot_type_ = key[boot_type];
    int sim_dist_ = key[sim_dist];

    // Call Fortran simulation routine.
    int n_conf = conf.nrow();
    int m_conf = conf.ncol();
    Rcpp::NumericMatrix fort_res(n_boot * n_conf, m_conf + 1);
    mack_sim(n_dev, triangle.begin(), sim_type_, boot_type_, sim_dist_, n_boot,
             n_conf, m_conf, conf.begin(), fort_res.begin(), show_progress,
             seed);
    res_list[boot_type_ - 1] = fort_res;
  }
  Rcpp::DataFrame res = mackPost(res_list, boot_types, sim_type);
  return (res);
};