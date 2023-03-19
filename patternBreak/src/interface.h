#ifndef INTERFACE
#define INTERFACE

extern "C" {

  void mack_boot_(int n_dev, double* triangle, int boot_type, int process_dist, bool conditional, int resids_type, int n_boot, double* reserve);

  void glm_boot_f(int n_dev, double* triangle, int n_boot, double* reserve);

  void mack_sim_(int n_dev, double* triangle, int sim_type, int n_boot, int n_factors, double* factors, int boot_type, int process_dist, bool conditional, int resids_type, bool show_progress, int n_sim, int m_sim, double* results);

  void glm_sim_f(int n_dev, double* triangle, int n_config, int m_config, double* config, int type, int n_boot, double* results);

  int validate_rng_f(int n_samples);
}

#endif