#ifndef MPC_FORMULATE_HPP
#define MPC_FORMULATE_HPP

#include "ip.hpp"

double mpc_objective (double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES]);
double mpc_objective_end_term(double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES]);
double[NUM_STATE_VARIABLES] mpc_initial_value (void);
double[NUM_STATE_VARIABLES] mpc_state_differential (double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES]);
double[NUM_PATH_CONSTRAINTS] mpc_path_constraints (double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES]);
double mpc_terminal_constraints(double State_variables[NUM_STATE_VARIABLES]);


#endif // MPC_FORMULATE_HPP
