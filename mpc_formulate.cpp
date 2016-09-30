#ifndef MPC_FORMULATE_CPP
#define MPC_FORMULATE_CPP

#include "mpc_formulate.hpp"

/****************************
	input: x(t), u(t)
	return: f(x(t), u(t))
	
****************************/
double mpc_objective (double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES])
{
	double objective_value = 0;
}

/****************************

	input: x(T)
	return: E(x(T))
	
****************************/

double mpc_objective_end_term(double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES])
{
}


/****************************

	input: void
	output: state_variables_initial_value
	
****************************/


double[NUM_STATE_VARIABLES] mpc_initial_value (void)
{
	double[NUM_STATE_VARIABLES] State_variables_initial_value;
}

/****************************

	Input: x(t), u(t)
	return: dot_x(t)
	
****************************/

double[NUM_STATE_VARIABLES] mpc_state_differential (double State_variables[NUM_STATE_VARIABLES], double Control_variables[NUM_CONTROL_VARIABLES])
{
	int i = 0;
	double dot_State_variables [NUM_STATE_VARIABLES];

	dot_State_variables[i++] = 	
	
}


/****************************

	input: x(t), u(t)
	return: C_I(x(t), u(t)) where C_I(x) >= 0
	
****************************/

double[NUM_PATH_CONSTRAINTS] mpc_path_constraints (double[NUM_STATE_VARIABLES] State_variables, double[NUM_CONTROL_VARIABLES] Control_variables)
{
	double Constraints [NUM_PATH_CONSTRAINTS];	
}


/****************************

	input: x(T)
	return: r(x(T)) where r(x) >= 0
	
****************************/

double mpc_terminal_constraints(double[NUM_STATE_VARIABLES] State_variables)
{
	
}

#endif // MPC_FORMULATE_CPP

