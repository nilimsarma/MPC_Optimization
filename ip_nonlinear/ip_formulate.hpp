#ifndef MPC_DISCRETIZE_HPP
#define MPC_DISCRETIZE_HPP

#include "ip.hpp"

double Compute_objective_value(const struct_ip_vars &s_ip_vars);
void Compute_inequality_constraints(const struct_ip_vars &s_ip_vars, double Inequality_constraints[NUM_INEQUALITY_CONSTRAINTS]);

#endif // MPC_DISCRETIZE_HPP
