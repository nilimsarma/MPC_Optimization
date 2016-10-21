#ifndef UTILS_HPP
#define UTILS_HPP

#include "ip.hpp"

void print_matrix(double* A, int num_rows, int num_cols);
void print_vector(double* A, int vec_size);
void print_ip_vars(const struct_ip_vars &s_ip_vars);
void print_problem_parameters(void);
void print_mu_value(const double &mu);
void print_nu_value(const double &nu);
void print_alpha_value(const double &alpha);
void print_alpha_max_value(const double &alpha_max);
void print_error_value(const double &error, const double &mu);

#endif //UTILS_HPP
