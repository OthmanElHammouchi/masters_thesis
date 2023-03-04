#include <math.h>
#include <Rcpp.h>

#include <RcppThread.h>
#include <dust/random/random.hpp>
#include <trng/special_functions.hpp>
#include <boost/math/special_functions/erf.hpp>

using namespace Rcpp;

extern "C" {

  void mack_boot_f(int n_dev, double* triangle, int resids_type, int boot_type, int dist, int n_boot, double* reserve);

  void glm_boot_f(int n_dev, double* triangle, int n_boot, double* reserve);

  void mack_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  void glm_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  void* init_rng(int n_threads, int seed) {
    using rng_state_type = dust::random::generator<double>;
    void* rng_ptr = (void*) new dust::random::prng<rng_state_type>(n_threads, seed);
    return(rng_ptr);
  }

  double rnorm_par(void* rng_ptr, int i_thread, double mean, double sd) {
    using rng_state_type = dust::random::generator<double>;
    using prng = dust::random::prng<rng_state_type>;
    prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
    auto& state = rng_ptr_->state(0);
    double sample = dust::random::normal<double>(state, mean, sd);
    return(sample);
}

  double rgamma_par(void* rng_ptr, int i_thread, double shape, double scale) {
    using rng_state_type = dust::random::generator<double>;
    using prng = dust::random::prng<rng_state_type>;
    prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
    auto& state = rng_ptr_->state(0);    
    double sample = dust::random::gamma<double>(state, shape, scale);
    return(sample);
  }

  double runif_par(void* rng_ptr, int i_thread) {
    using rng_state_type = dust::random::generator<double>;
    using prng = dust::random::prng<rng_state_type>;
    prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
    auto& state = rng_ptr_->state(0);    
    double sample = dust::random::random_real<double>(state);
    return(sample);
  }

  double pois_cdf(double x, double lambda) {
    double p;
    for (int i; i < x; i++) {
      p += pow(lambda, i) * std::exp(-lambda) / std::tgamma(i + 1);
    };
    return(p);
  }

  long int rpois_large_mean(double u, double lambda) {
    long double w = boost::math::erf_inv(u);
    long double x, delta;
    if (std::abs(w) < 3) {
      long double Q1 = lambda + std::sqrt(lambda) * w + (1/3 + (1/6) * pow(w, 2));
      x = Q1 + (1/std::sqrt(lambda)) * ((-1/36) * w - (1/72) * pow(w, 3));
      delta = (1/40 + (1/80) * pow(w, 2) + (1/160) * pow(w, 4));
    } else {
      long double r = 1 + w / std::sqrt(lambda) + (1/6)*pow((w/std::sqrt(lambda)), 2) - (1/72) * pow(w / std::sqrt(lambda), 3);
      long double c0 = 1/3 - (1/36) * (r - 1);
      x = lambda * r + c0;
      x -= (4.1/805)/(x + 0.025 * lambda);
    }
    long int n = std::floor(x + delta);
    if (x > 10) {
      if ((x - n) > delta) {
        return(n);
      } else if (pois_cdf(n, lambda) < u) {
        return(n);
      }
      else {
        return(n - 1);
      }
    }
  }

  long int rpois_par(void* rng_ptr, int i_thread, double mean) {
    if (mean < 1e7) {
      using rng_state_type = dust::random::generator<double>;
      using prng = dust::random::prng<rng_state_type>;
    prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
    auto& state = rng_ptr_->state(0);      
    int sample = dust::random::poisson<double>(state, mean);
      return(sample);
    } else {
      double u = runif_par(rng_ptr, i_thread);
      long int sample = rpois_large_mean(u, mean);
      return(sample);
    };
  }

  void check_user_input(void) {
    RcppThread::checkUserInterrupt();
  }

  void* pgbar_create(int total, int freq) {
    void* ptr = (void*) new RcppThread::ProgressBar(total, freq);
    return(ptr);
  }

  void pgbar_incr(void* progress_bar) {
    (*static_cast<RcppThread::ProgressBar*>(progress_bar))++;
  }
}

double** array_create(int n, int m) {
    double* data = (double*)malloc(n * m * sizeof(double));
    double** array = (double**)malloc(n*sizeof(double*));
    for (int i=0; i<n; i++) {
      array[i] = &(data[m * i]);
    }
    return(array);
}

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = "mackBoot")]]
NumericVector mack_boot(NumericMatrix triangle, int n_boot, int resids_type, int boot_type, int dist) {


    int n_dev = triangle.rows();
    
    double** triangle_f = array_create(n_dev, n_dev);

    for (int i = 0; i < n_dev; i++) {
      for (int j = 0; j < n_dev; j++) {
        triangle_f[j][i] = triangle(i, j);
      }
    }

    NumericVector reserve(n_boot);

    mack_boot_f(n_dev, &triangle_f[0][0], resids_type, boot_type, dist, n_boot, reserve.begin());

    free(triangle_f);

    return(reserve);
};

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = "glmBoot")]]
NumericVector glm_boot(NumericMatrix triangle, int n_boot) {

    int n_dev = triangle.rows();
    
    double** triangle_f = array_create(n_dev, n_dev);

    for (int i = 0; i < n_dev; i++) {
      for (int j = 0; j < n_dev; j++) {
        triangle_f[j][i] = triangle(i, j);
      }
    }

    NumericVector reserve(n_boot);

    glm_boot_f(n_dev, &triangle_f[0][0], n_boot, reserve.begin());

    free(triangle_f);

  return(reserve);
};

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".mackSim")]]
NumericMatrix mack_sim(NumericMatrix triangle, int n_boot, NumericMatrix config, int type) {

  int n_dev = triangle.rows();

  int n_config = config.rows();
  int m_config = config.cols();

  double** triangle_f = array_create(n_dev, n_dev);

  for (int i = 0; i < n_dev; i++) {
    for (int j = 0; j < n_dev; j++) {
      triangle_f[j][i] = triangle(i, j);
    }
  }

  double** config_f = array_create(m_config, n_config);

  for (int i=0; i < n_config; i++) {
    for (int j=0; j < m_config; j++) {
      config_f[j][i] = config(i, j);
    }
  }

  int n_results = n_boot * n_config;
  int m_results = m_config + 1;

  double** results_f = array_create(m_results, n_results);

  mack_sim_f(n_dev, &triangle_f[0][0], n_config, m_config, &config_f[0][0], type, n_boot, &results_f[0][0]);

  NumericMatrix results(n_results, m_results);

    for (int i=0; i < n_results; i++) {
      for (int j=0; j < m_results; j++) {
        results(i, j) = results_f[j][i];
      }
    }


  free(triangle_f);
  free(config_f);
  free(results_f);

  return(results);
};

//' Leading NA
//' 
//' This function returns a logical vector identifying if 
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export(name = ".glmSim")]]
NumericMatrix glm_sim(NumericMatrix triangle, int n_boot, NumericMatrix config, int type) {

  int n_dev = triangle.rows();

  int n_config = config.rows();
  int m_config = config.cols();

  double** triangle_f = array_create(n_dev, n_dev);

  for (int i = 0; i < n_dev; i++) {
    for (int j = 0; j < n_dev; j++) {
      triangle_f[j][i] = triangle(i, j);
    }
  }

  double** config_f = array_create(m_config, n_config);

  for (int i=0; i < n_config; i++) {
    for (int j=0; j < m_config; j++) {
      config_f[j][i] = config(i, j);
    }
  }

  int n_results = n_boot * n_config;
  int m_results = m_config + 1;

  double** results_f = array_create(m_results, n_results);

  glm_sim_f(n_dev, &triangle_f[0][0], n_config, m_config, &config_f[0][0], type, n_boot, &results_f[0][0]);

  NumericMatrix results(n_results, m_results);

    for (int i=0; i < n_results; i++) {
      for (int j=0; j < m_results; j++) {
        results(i, j) = results_f[j][i];
      }
    }

  free(triangle_f);
  free(config_f);
  free(results_f);

  return(results);
};