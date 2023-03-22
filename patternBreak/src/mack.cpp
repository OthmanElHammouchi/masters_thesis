#include <algorithm>
#include <vector>
#include <string>

#include "helpers.h"
#include "constants.h"

extern "C" {
  void mack_boot_(int n_dev, double* triangle, int boot_type, int process_dist, bool conditional, int resids_type, int n_boot, double* reserve);

  void mack_sim_(int n_dev, double* triangle, int sim_type, int n_boot, int n_res, int m_res, double* results, bool show_progress);

  void build_sim_table_(int sim_type, int n_dev, int n_boot, int n_factors, double* factors, int n_boot_types, int* boot_types, int n_proc_dists, int* proc_dists, int n_conds, bool* conds, int n_resids_types, int* resids_types, int* n_res, int* m_res);
}

//' Simulate Mack CL reserve.
//'
//' @param triangle Cumulative claims triangle
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mackBoot(Rcpp::NumericMatrix triangle, int n_boot, Rcpp::String boot_type, Rcpp::String process_dist, bool cond, Rcpp::Nullable<Rcpp::String> resids_type_ = R_NilValue) {
  int n_dev = triangle.rows();
  triangle = na2zero(triangle);

  Rcpp::String resids_type;
  if (resids_type_.isNotNull()) {
    resids_type = resids_type_;
    if (boot_type != "residuals") {
      Rcpp::stop("Argument 'resids_type' only valid for residuals bootstrap.");
      }
    } else {
      resids_type = "none";
    }

    if (boot_type == "pairs" && !cond) {
      Rcpp::stop("Unconditional method invalid for pairs bootstrap.");
    }

  Rcpp::NumericVector reserve(n_boot);
  mack_boot_(n_dev, triangle.begin(), key[boot_type], key[process_dist], cond, key[resids_type], n_boot, reserve.begin());
  return(reserve);
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
//' @details
//' The simulation configuration inputs `sim_type`, `factor`, `boot_type`, `proc_dist` and `conditional` can be
//' either strings or character vectors. In the latter case, the simulation is computed for all feasible combinations.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame mackSim(Rcpp::NumericMatrix triangle, Rcpp::String sim_type, int n_boot, Rcpp::NumericVector factor, Rcpp::CharacterVector boot_type, Rcpp::CharacterVector proc_dists, Rcpp::LogicalVector conds, Rcpp::Nullable<Rcpp::CharacterVector> resids_type = R_NilValue, bool show_progress = true) {

  int n_dev = triangle.rows();
  triangle = na2zero(triangle);

  Rcpp::CharacterVector resids_types__;
  if (resids_type.isNotNull()) {
    if (!contains_str(boot_type, "residuals")) {
      Rcpp::stop("Argument 'resids_type' only valid for residuals bootstrap.");
    }
    resids_types__ = resids_type;
  } else {
    resids_types__ = "none";
  }

  if (boot_type.length() == 1 && boot_type(0) == "pairs") {
    if (conds.length() == 1 && !conds(0)) {
      Rcpp::stop("Unconditional method invalid for pairs bootstrap.");
    }
  }

  // Prepare arguments to be passed to Fortran.
  int n_factors = factor.length();
  int n_boot_types = boot_type.length();
  int n_proc_dists = proc_dists.length();
  int n_conds = conds.length();
  int n_resids_types = resids_types__.length();

  int sim_type_;
  if (sim_type == "single") {
    sim_type_ = SINGLE;
  } else if (sim_type == "calendar") {
    sim_type_ = CALENDAR;
  } else {
    sim_type_ = ORIGIN;
  }

  int boot_types_[n_boot_types];
  for (int i = 0; i < n_boot_types; i++) {
    if (boot_type(i) == "parametric") {
      boot_types_[i] = PARAMETRIC;
    } else if (boot_type(i) == "residuals") {
      boot_types_[i] = RESID;
    } else {
      boot_types_[i] = PAIRS;
    }
  }

  int proc_dists_[n_proc_dists];
  for (int i = 0; i < n_proc_dists; i++) {
    if (proc_dists(i) == "normal") {
      proc_dists_[i] = NORMAL;
    } else if (proc_dists(i) == "gamma") {
      proc_dists_[i] = GAMMA;
    }
  }

  bool conds_[n_conds];
  for (int i = 0; i < n_conds; i++) {
    if (conds(i)) {
      conds_[i] = true;
    } else {
      conds_[i] = false;
    }
  }

  int resids_types_[n_resids_types];
  for (int i = 0; i < n_resids_types; i++) {
    if (resids_types__(i) == "standardised") {
      resids_types_[i] = STANDARDISED;
    } else if (resids_types__(i) == "modified") {
      resids_types_[i] = MODIFIED;
    } else if (resids_types__(i) == "studentised") {
      resids_types_[i] = STUDENTISED;
    } else {
      resids_types_[i] = LOGNORMAL;
    }
  }
  
  // Call Fortran simulation routine.
  int n_res, m_res;
  build_sim_table_(sim_type_, n_dev, n_boot, n_factors, factor.begin(), n_boot_types, boot_types_, n_proc_dists, proc_dists_, n_conds, conds_, n_resids_types, resids_types_, &n_res, &m_res);

  Rcpp::NumericMatrix res_(n_res, m_res);
  mack_sim_(n_dev, triangle.begin(), sim_type_, n_boot, n_res, m_res, res_.begin(), show_progress);

  // Convert result to dataframe with proper column names and descriptive elements.
  Rcpp::Function asDataFrame("as.data.frame");
  Rcpp::DataFrame res = asDataFrame(res_);
  if (sim_type == "single") {
    res.attr("names") = Rcpp::CharacterVector({"outlier.rowidx", "outlier.colidx", "factor", "excl.rowidx", "excl.colidx", "boot.type", "proc.dist", "cond", "resids.type", "reserve"});
  } else if (sim_type == "calendar") {
      res.attr("names") = Rcpp::CharacterVector({"outlier.diagidx", "factor", "excl.diagidx", "boot.type", "proc.dist", "cond", "resids.type", "reserve"});    
  } else if (sim_type == "origin") {
    res.attr("names") = Rcpp::CharacterVector({"outlier.rowidx", "factor", "excl.rowidx", "boot.type", "proc.dist", "cond", "resids.type", "reserve"});
  }

  Rcpp::CharacterVector boot_col = Rcpp::as<Rcpp::CharacterVector>(res["boot.type"]);
  std::replace(boot_col.begin(), boot_col.end(), std::string("1"), std::string("parametric"));
  std::replace(boot_col.begin(), boot_col.end(), std::string("2"), std::string("residuals"));
  std::replace(boot_col.begin(), boot_col.end(), std::string("3"), std::string("pairs"));
  res["boot.type"] = Rcpp::CharacterVector(boot_col.begin(), boot_col.end());

  Rcpp::CharacterVector proc_col = Rcpp::as<Rcpp::CharacterVector>(res["proc.dist"]);
  std::replace(proc_col.begin(), proc_col.end(), std::string("1"), std::string("normal"));
  std::replace(proc_col.begin(), proc_col.end(), std::string("2"), std::string("gamma"));
  res["proc.dist"] = Rcpp::CharacterVector(proc_col.begin(), proc_col.end());

  res["cond"] = Rcpp::as<Rcpp::LogicalVector>(res["cond"]);

  Rcpp::CharacterVector resids_col = Rcpp::as<Rcpp::CharacterVector>(res["resids.type"]);
  std::replace(resids_col.begin(), resids_col.end(), std::string("0"), std::string(Rcpp::String(NA_STRING)));
  std::replace(resids_col.begin(), resids_col.end(), std::string("1"), std::string("standardised"));
  std::replace(resids_col.begin(), resids_col.end(), std::string("2"), std::string("modified"));
  std::replace(resids_col.begin(), resids_col.end(), std::string("3"), std::string("studentised"));
  std::replace(resids_col.begin(), resids_col.end(), std::string("4"), std::string("log-normal"));
  res["resids.type"] = Rcpp::CharacterVector(resids_col.begin(), resids_col.end());

  return(res);
};
