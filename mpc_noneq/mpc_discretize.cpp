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

#ifndef TEST
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
	
	double Equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_equality_constraints(s_ip_vars, Equality_constraints);

	p = 0;
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		Inequality_constraints[p++] = Equality_constraints[i];
		Inequality_constraints[p++] = -Equality_constraints[i];
	}

	m = 0;
	for(i = 0; i < NUM_DISCRETIZATION; i++)
	{			
		double State_variables[NUM_STATE_VARIABLES];
		for(j = 0; j < NUM_STATE_VARIABLES; j++)	State_variables[j] = s_ip_vars.Optimization_variables[m++];

		double Control_variables[NUM_CONTROL_VARIABLES];
		for(j = 0; j < NUM_CONTROL_VARIABLES; j++)	Control_variables[j] = s_ip_vars.Optimization_variables[m++];

		double path_constraints[NUM_PATH_CONSTRAINTS];
		mpc_path_constraints(State_variables,Control_variables, path_constraints);
		for(j = 0 ; j < NUM_PATH_CONSTRAINTS; j++)	Inequality_constraints[p++] = path_constraints[j];

#if (NUM_TERMINAL_CONSTRAINTS > 0)	
		double terminal_constraints[NUM_TERMINAL_CONSTRAINTS];
		mpc_terminal_constraints(State_variables, terminal_constraints);
		for(j = 0 ; j < NUM_TERMINAL_CONSTRAINTS; j++)	Inequality_constraints[p++] = terminal_constraints[j];
#endif

	}
}

#else	// ifndef TEST

double Compute_objective_value(const struct_ip_vars &s_ip_vars)
{
	double objective_value = 0;

	double x = s_ip_vars.Optimization_variables[0];
	double y = s_ip_vars.Optimization_variables[1];
	
	objective_value = x*x + 2*y*y - 2*x - 6*y - 2*x*y;

	return objective_value;
}

void Compute_inequality_constraints(const struct_ip_vars &s_ip_vars, double Inequality_constraints[NUM_INEQUALITY_CONSTRAINTS])
{
	double x = s_ip_vars.Optimization_variables[0];
	double y = s_ip_vars.Optimization_variables[1];

	Inequality_constraints[0] = 1 - x/2 - y/2;
	Inequality_constraints[1] = 2 + x - 2*y;
	Inequality_constraints[2] = x;
	Inequality_constraints[3] = y;
//	Inequality_constraints[4] = x + y - 2;
//	Inequality_constraints[5] = 2 - x - y;
}
#endif // ifndef TEST

#endif // MPC_DISCRETIZE_CPP
