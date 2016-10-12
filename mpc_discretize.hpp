#ifndef MPC_DISCRETIZE_HPP
#define MPC_DISCRETIZE_HPP

#include "ip.hpp"
#include "mpc_formulate.hpp"

void Integrator(const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double State_vars_updated[NUM_STATE_VARIABLES]);
double Compute_objective_value(const struct_ip_vars &s_ip_vars);
void Compute_equality_constraints(const struct_ip_vars &s_ip_vars, double Equality_constraints[NUM_EQUALITY_CONSTRAINTS]);
void Compute_inequality_constraints(const struct_ip_vars &s_ip_vars, double Inequality_constraints[NUM_INEQUALITY_CONSTRAINTS]);

#endif // MPC_DISCRETIZE_HPP
