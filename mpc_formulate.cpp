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
	int i = 0;
	
	objective_value += (10.0 * State_variables[i] * State_variables[i]); i++;
	objective_value += (10.0 * State_variables[i] * State_variables[i]); i++;
	objective_value += (1.0 * State_variables[i] * State_variables[i]); i++;
	objective_value += (1.0 * State_variables[i] * State_variables[i]); i++;
	objective_value += (1.0e-8 * Control_variables[0] * State_variables[0]);

	return objective_value;
}

/****************************

	input: x(T)
	return: E(x(T))
	
****************************/

double mpc_objective_end_term(const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES])
{
	return 0;
}


/****************************

	input: void
	output: state_variables_initial_value
	
****************************/


void mpc_initial_value (double State_variables_initial_value[NUM_STATE_VARIABLES])
{
	int i;
	for(i = 0; i < NUM_STATE_VARIABLES; i++)
		State_variables_initial_value[i] = 0.5;
}

/****************************

	Input: x(t), u(t)
	return: dot_x(t)
	
****************************/

void mpc_state_differential (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double dot_State_variables [NUM_STATE_VARIABLES])
{
	double mB, mW, kS, kT;
	mB = 350.0;
	mW = 50.0;
	kS = 20,000.0;
	kT = 200,000.0;
	
	int i = 0;
	dot_State_variables[i++] = 	State_variables[2];
	dot_State_variables[i++] = 	State_variables[3];
	dot_State_variables[i++] = 	(-kS*State_variables[0] + kS*State_variables[1] + Control_variables[0])/mB;
	dot_State_variables[i++] = 	(kS*State_variables[0] - (kT + kS)*State_variables[1] + kT*Control_variables[1] - Control_variables[0])/mW;
}


/****************************

	input: x(t), u(t)
	return: C_I(x(t), u(t)) where C_I(x) >= 0
	
****************************/

void mpc_path_constraints (const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double Constraints [NUM_PATH_CONSTRAINTS])
{
	int i = 0;
	
	Constraints[i++] = Control_variables[0] + 200.0;
	Constraints[i++] = -Control_variables[0] + 200.0;
	Constraints[i++] = Control_variables[1];
	Constraints[i++] = -Control_variables[1];
}


/****************************

	input: x(T)
	return: r(x(T)) where r(x) >= 0
	
****************************/

void mpc_terminal_constraints(const double State_variables[NUM_STATE_VARIABLES], double terminal_constraints[NUM_TERMINAL_CONSTRAINTS])
{
	int i;
	for(i = 0; i < NUM_TERMINAL_CONSTRAINTS; i++)	terminal_constraints[i] = 0;
}

#endif // MPC_FORMULATE_CPP

