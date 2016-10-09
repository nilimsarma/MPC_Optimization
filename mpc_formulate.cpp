#ifndef MPC_FORMULATE_CPP
#define MPC_FORMULATE_CPP

#include "mpc_formulate.hpp"
 
/****************************
	input: x(t), u(t)
	return: f(x(t), u(t))
	
****************************/
double mpc_objective (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES])
{
	double objective_value = 0;
}

/****************************

	input: x(T)
	return: E(x(T))
	
****************************/

double mpc_objective_end_term(const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES])
{
}


/****************************

	input: void
	output: state_variables_initial_value
	
****************************/


void mpc_initial_value (double[NUM_STATE_VARIABLES] State_variables_initial_value)
{
}

/****************************

	Input: x(t), u(t)
	return: dot_x(t)
	
****************************/

void mpc_state_differential (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double dot_State_variables [NUM_STATE_VARIABLES])
{
	int i = 0;
	dot_State_variables[i++] = 	
	
}


/****************************

	input: x(t), u(t)
	return: C_I(x(t), u(t)) where C_I(x) >= 0
	
****************************/

void mpc_path_constraints (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double Constraints [NUM_PATH_CONSTRAINTS])
{
	
}


/****************************

	input: x(T)
	return: r(x(T)) where r(x) >= 0
	
****************************/

double mpc_terminal_constraints(const double State_variables[NUM_STATE_VARIABLES])
{
	
}

#endif // MPC_FORMULATE_CPP

