#ifndef MPC_DISCRETIZE_CPP
#define MPC_DISCRETIZE_CPP

#include "mpc_discretize.hpp"

/****************************************************************

	Gauss Legendre Order 4 Integrator (Explicit) -- may change later to implicit 

	b1 = 1/2
	b2 = 1/2

	k1 = f( t[n], y[n] )
	k2 = f( t[n] + c2*h, y[n] + h*a21*k1 )

	y[n+1] = y[n] + h(b1*k1 + b2*k2)	

*****************************************************************/


/********************

	Input: x(t), u(t)
	output: x(t+1)

*********************/
void Integrator(const double State_variables[NUM_STATE_VARIABLES], const double Control_variables[NUM_CONTROL_VARIABLES], double State_vars_updated[NUM_STATE_VARIABLES])
{
	
	double h  = INTERVAL_SIZE;
	double b1 = 0.5;
	double b2 = 0.5;
	double c2 =  (0.50 + sqrt(3)/6.0);
	double a21 = (0.25 + sqrt(3)/6.0);
	
	double k1[NUM_STATE_VARIABLES], k2[NUM_STATE_VARIABLES];

	double State_vars	[NUM_STATE_VARIABLES];
	double Control_vars [NUM_CONTROL_VARIABLES];

	int i;
	for(i = 0; i < NUM_CONTROL_VARIABLES; i++)
	{
		Control_vars[i] = Control_variables[i];
	}
	
	for(i = 0; i < NUM_STATE_VARIABLES; i++)
	{
		State_vars[i] = State_variables[i];
	}

	mpc_state_differential(State_vars, Control_vars, k1);
	
	for(i = 0; i < NUM_STATE_VARIABLES; i++)
	{
		State_vars[i] = State_variables[i] + h*a21*k1[i];
	}
	
	mpc_state_differential(State_vars, Control_vars, k2);

	for(i = 0; i < NUM_STATE_VARIABLES; i++)
	{
		State_vars_updated[i] = State_vars[i] + h*(b1*k1[i] + b2*k2[i]);
	}
	
}

double Compute_objective_value(const struct_ip_vars &s_ip_vars)
{
	int i,j,k,m;

	double objective_value = 0;
	
	double State_variables[NUM_STATE_VARIABLES];
	double Control_variables[NUM_CONTROL_VARIABLES];

	m = 0;	
	for(i = 0; i < (NUM_DISCRETIZATION-1); i++)
	{
		for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j] = s_ip_vars.Optimization_variables[m++];
		for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];
		objective_value += mpc_objective(State_variables, Control_variables);
	}

	for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j] = s_ip_vars.Optimization_variables[m++];
	for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];	
	objective_value += mpc_objective_end_term(State_variables, Control_variables);

	return objective_value;
}

void Compute_equality_constraints(const struct_ip_vars &s_ip_vars, double Equality_constraints[NUM_EQUALITY_CONSTRAINTS])
{
	int i,j,k,m,p;

	double State_variables[NUM_STATE_VARIABLES];
	double Control_variables[NUM_CONTROL_VARIABLES];
	
	double State_vars_updated[NUM_STATE_VARIABLES];

	m = 0;	p = 0;

	//initial state variables
	
	for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j]	 = s_ip_vars.Optimization_variables[m++];
	for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];

	mpc_initial_value(State_vars_updated);
	
	for(j = 0; j < NUM_STATE_VARIABLES; j++)	Equality_constraints[p++] = State_variables[j] - State_vars_updated[j];

	//updated state variables using integrator
	for(i = 1; i < NUM_DISCRETIZATION; i++)
	{
		Integrator(State_variables, Control_variables, State_vars_updated);

		for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j] 	 = s_ip_vars.Optimization_variables[m++];
		for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];

		for(j = 0; j < NUM_STATE_VARIABLES; j++)	Equality_constraints[p++] = State_variables[j] - State_vars_updated[j];
		
	}
}

void Compute_inequality_constraints(const struct_ip_vars &s_ip_vars, double Inequality_constraints[NUM_INEQUALITY_CONSTRAINTS])
{
	int i,j,k,m,p;

	double path_constraints[NUM_PATH_CONSTRAINTS];
	double terminal_constraints[NUM_TERMINAL_CONSTRAINTS];

	double State_variables[NUM_STATE_VARIABLES];
	double Control_variables[NUM_CONTROL_VARIABLES];

	m = 0; p = 0;
	for(i = 0; i < NUM_DISCRETIZATION; i++)
	{		
		for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j] 	 = s_ip_vars.Optimization_variables[m++];
		for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];

		mpc_path_constraints(State_variables,Control_variables, path_constraints);
		mpc_terminal_constraints(State_variables, terminal_constraints);

		for(j = 0 ; j < NUM_PATH_CONSTRAINTS; j++)	Inequality_constraints[p++] = path_constraints[j];
		for(j = 0 ; j < NUM_TERMINAL_CONSTRAINTS; j++)	Inequality_constraints[p++] = terminal_constraints[j];
	}
}

#endif // MPC_DISCRETIZE_CPP
