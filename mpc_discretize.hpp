#ifndef MPC_DISCRETIZE_HPP
#define MPC_DISCRETIZE_HPP

#include <stdio.h>
#include <stdlib.h> 


double[NUM_STATE_VARIABLES] Integrator(double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES]);
double Compute_objective_value(struct_ip_vars s_ip_vars);
double[NUM_EQUALITY_CONSTRAINTS] Compute_equality_constraints(struct_ip_vars s_ip_vars);
double[NUM_INEQUALITY_CONSTRAINTS] Compute_inequality_constraints(struct_ip_vars s_ip_vars);

#endif // MPC_DISCRETIZE_HPP
