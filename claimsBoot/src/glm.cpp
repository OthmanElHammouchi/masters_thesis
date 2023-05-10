#include "constants.h"
#include "helpers.h"
#include <Rcpp.h>

extern "C" {

void glm_boot(int n_dev, double *triangle, int boot_type, int opt, int n_boot,
              double *reserve, int seed);

void glm_sim(int n_dev, double *triangle, int sim_type, int boot_type,
             int n_conf, int m_conf, double *config, int n_boot,
             double *results, bool show_progress, int seed);

int validate_rng_(int n_samples);
}

//' Simulate Poisson GLM reserve.
//'
//' @param triangle Incremental claims triangle
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector glmBoot(Rcpp::NumericMatrix triangle, int n_boot,
                            Rcpp::String boot_type, Rcpp::String opt,
                            int seed = 42) {
  triangle = na_to_zero(triangle);
  int n_dev = triangle.cols();

  int boot_type_ = key[boot_type];
  int opt_ = key[opt];

  Rcpp::NumericVector reserve(n_boot);
  glm_boot(n_dev, triangle.begin(), boot_type_, n_boot, opt_, reserve.begin(),
           seed);
  return (reserve);
};

//' Simulate Poisson GLM reserve for different perturbed and excluded points.
//'
//' @param triangle Incremental claims triangle
//' @param sim_type Simulation type: `"single"` (the default), `"calender"`, or `"origin"`
//' @param n_boot Number of bootstrap simulations
//' @param factor Perturbation factor
//' @param boot_type Type of bootstrap: `"parametric"`, `"residuals"`, or `"pairs"`
//' @param sim_dist Distribution of process error: `"normal"` or `"gamma"`
//'
//' @details
//' The simulation configuration inputs `sim_type`, `factor`, `boot_type`,
//`proc_dist` and `conditional` can be ' either strings or character vectors. In
// the latter case, the simulation is computed for all feasible combinations.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame glmSim(Rcpp::NumericMatrix triangle, Rcpp::String sim_type,
                       int n_boot, Rcpp::NumericVector factors,
                       Rcpp::CharacterVector boot_types,
                       Rcpp::String sim_dist = "normal",
                       bool show_progress = true, int seed = 42) {
  triangle = na_to_zero(triangle);
  int n_dev = triangle.cols();

  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("claimsBoot");
  Rcpp::Function glmConfig = pkg["glmConfig"];
  Rcpp::Function glmPost = pkg["glmPost"];

  Rcpp::List res_list(boot_types.length());
  for (Rcpp::CharacterVector::iterator i = boot_types.begin();
       i != boot_types.end(); ++i) {
    Rcpp::String boot_type = *i;
    Rcpp::NumericMatrix conf = glmConfig(n_dev, factors, sim_type, boot_type);
    int n_conf = conf.rows();
    int m_conf = conf.cols();

    int sim_type_ = key[sim_type];
    int boot_type_ = key[boot_type];
    int sim_dist_ = key[sim_dist];

    Rcpp::NumericMatrix fort_res(n_conf * n_boot, m_conf + 1);
    glm_sim(n_dev, triangle.begin(), sim_type_, boot_type_, n_conf, m_conf,
            conf.begin(), n_boot, fort_res.begin(), show_progress, seed);
    int idx = i - boot_types.begin();
    res_list[idx] = fort_res;
  }
  Rcpp::DataFrame res = glmPost(res_list, boot_types, sim_type);
  return (res);
};