#ifndef IP_FORMULATE_CPP
#define IP_FORMULATE_CPP

#include "ip_formulate.hpp"

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
}

#endif // IP_FORMULATE_CPP