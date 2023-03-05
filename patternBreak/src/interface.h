#ifndef INTERFACE
#define INTERFACE

extern "C" {

  void mack_boot_f(int n_dev, double* triangle, int resids_type, int boot_type, int dist, int n_boot, double* reserve);

  void glm_boot_f(int n_dev, double* triangle, int n_boot, double* reserve);

  void mack_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  void glm_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  int validate_rng_f(int n_samples);
}

#endif