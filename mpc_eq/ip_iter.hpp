#ifndef IP_ITER_HPP
#define IP_ITER_HPP

#include "ip.hpp"
#include "mpc_discretize.hpp"
#include "ip_primal_dual_dir.hpp"
#include "utils.hpp"

struct_alpha Compute_alpha_max(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir);
double Compute_merit_function(const struct_ip_vars &s_ip_vars);
double Compute_directional_derivative_merit_function(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir);
struct_alpha Backtrack_alpha(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir, const struct_alpha &s_alpha_max);
void Update_IP_vars(struct_ip_vars &s_ip_vars, const struct_alpha &s_alpha, const struct_primal_dual_direction &s_primal_dual_dir);
double Compute_Error(const struct_ip_vars &s_ip_vars, const double mu);
void ip_init(struct_ip_vars &s_ip_vars);
void Compute_mu(struct_ip_vars &s_ip_vars);

#endif //IP_ITER_HPP
