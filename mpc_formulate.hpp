#ifndef MPC_FORMULATE_HPP
#define MPC_FORMULATE_HPP

#include "ip.hpp"

double mpc_objective (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES]);
double mpc_objective_end_term(const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES]);
void mpc_initial_value (double State_variables_initial_value[NUM_STATE_VARIABLES]);
void mpc_state_differential (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double dot_State_variables [NUM_STATE_VARIABLES]);
void mpc_path_constraints (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double Constraints [NUM_PATH_CONSTRAINTS]);
double mpc_terminal_constraints(const double State_variables[NUM_STATE_VARIABLES]);

#endif // MPC_FORMULATE_HPP
