#ifndef MPC_DISCRETIZE_HPP
#define MPC_DISCRETIZE_HPP

#include "ip.hpp"

#define HESSIAN_LAGRANGIAN_SIZE 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_INEQUALITIES_NUM_ROWS	NUM_INEQUALITY_CONSTRAINTS
#define JACOBIAN_INEQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_EQUALITIES_NUM_ROWS	NUM_EQUALITY_CONSTRAINTS
#define JACOBIAN_EQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

void Compute_objective_value(const struct_ip_vars &s_ip_vars, double &objective_value);
void Compute_inequality_constraints(const struct_ip_vars &s_ip_vars, double Inequality_constraints[NUM_INEQUALITY_CONSTRAINTS]);
void Compute_equality_constraints(const struct_ip_vars &s_ip_vars, double Equality_constraints[NUM_EQUALITY_CONSTRAINTS]);

#endif // MPC_DISCRETIZE_HPP
