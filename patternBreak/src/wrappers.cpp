#include <stdlib.h>
#include <math.h>

#include <Rcpp.h>
#include <RcppThread.h>
#include <dust/random/random.hpp>
#include <boost/math/special_functions/erf.hpp>

extern "C" {
  void* init_rng(int n_threads, int seed);

  double rnorm_par(void* rng_ptr, int i_thread, double mean, double sd);

  double rgamma_par(void* rng_ptr, int i_thread, double shape, double scale);

  double runif_par(void* rng_ptr, int i_thread);

  long int rpois_par(void* rng_ptr, int i_thread, double mean);

  void* pgbar_create(int total, int freq);

  void pgbar_incr(void* progress_bar);

  void rprint_par(char* str);

  void check_user_input(void);
}

void* init_rng(int n_threads, int seed) {
  using rng_state_type = dust::random::generator<double>;
  void* rng_ptr = (void*) new dust::random::prng<rng_state_type>(n_threads, seed);
  return(rng_ptr);
}

double rnorm_par(void* rng_ptr, int i_thread, double mean, double sd) {
  using rng_state_type = dust::random::generator<double>;
  using prng = dust::random::prng<rng_state_type>;
  prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
  auto& state = rng_ptr_->state(i_thread);
  next(state);
  double sample = dust::random::normal<double>(state, mean, sd);
  return(sample);
}

double rgamma_par(void* rng_ptr, int i_thread, double shape, double scale) {
  using rng_state_type = dust::random::generator<double>;
  using prng = dust::random::prng<rng_state_type>;
  prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
  auto& state = rng_ptr_->state(i_thread);
  next(state);
  double sample = dust::random::gamma<double>(state, shape, scale);
  return(sample);
}

double runif_par(void* rng_ptr, int i_thread) {
  using rng_state_type = dust::random::generator<double>;
  using prng = dust::random::prng<rng_state_type>;
  prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
  auto& state = rng_ptr_->state(i_thread);
  next(state);
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
    delta = (1/40 + (1/80) * pow(w, 2) + (1/160) * pow(w, 4)) / lambda;
  } else {
    long double r = 1 + w / std::sqrt(lambda) + (1/6)*pow((w/std::sqrt(lambda)), 2) - (1/72) * pow(w / std::sqrt(lambda), 3);
    long double c0 = 1/3 - (1/36) * (r - 1);
    x = lambda * r + c0;
    x -= (4.1/805) / (x + 0.025 * lambda);
    delta = 0.01 / lambda;
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
  return(n);
}

long int rpois_par(void* rng_ptr, int i_thread, double mean) {
  if (mean < 1e7) {
    using rng_state_type = dust::random::generator<double>;
    using prng = dust::random::prng<rng_state_type>;
  prng* rng_ptr_ = static_cast<prng*>(rng_ptr);
  auto& state = rng_ptr_->state(i_thread);
  next(state);     
  int sample = dust::random::poisson<double>(state, mean);
    return(sample);
  } else {
    double u = runif_par(rng_ptr, i_thread);
    long int sample = rpois_large_mean(u, mean);
    return(sample);
  };
}

// [[Rcpp::export]]
Rcpp::NumericVector test_pois(int n, double lambda) {
  using rng_state_type = dust::random::generator<double>;
  void* rng_ptr = (void*) new dust::random::prng<rng_state_type>(2, 42);
  Rcpp::NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res(i) = rpois_par(rng_ptr, 1, lambda);
  }
  return(res);
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

void rprint_par(char* str) {
  RcppThread::Rcout << str;
}