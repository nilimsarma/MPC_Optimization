#ifndef IP_FORMULATE_CPP
#define IP_FORMULATE_CPP

#include "ip_formulate.hpp"

void Compute_objective_value(const struct_ip_vars &s_ip_vars, double &objective_value)
{
	double objective_value = 0;

	double x = s_ip_vars.Optimization_variables[0];
	double y = s_ip_vars.Optimization_variables[1];
	
	objective_value = x*x + 2*y*y - 2*x - 6*y - 2*x*y;
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

void Compute_equality_constraints(const struct_ip_vars &s_ip_vars, double Equality_constraints[NUM_EQUALITY_CONSTRAINTS])
{
	double x = s_ip_vars.Optimization_variables[0];
	double y = s_ip_vars.Optimization_variables[1];

	Equality_constraints[0] = x + y - 2;
}

void Compute_Jacobian_Equalities(const struct_ip_vars &s_ip_vars, double Jacobian_Equalities[JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS])
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0[NUM_EQUALITY_CONSTRAINTS], f1[NUM_EQUALITY_CONSTRAINTS];

	s_ip_vars_1 = s_ip_vars;
	for(j = 0; j < JACOBIAN_EQUALITIES_NUM_COLS; j++)
	{
		//f0
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		Compute_equality_constraints(s_ip_vars_1, f0);

		//f1
		s_ip_vars_1.Optimization_variables[j] -= (2*ESP_DIFFERENTIATION);
		Compute_equality_constraints(s_ip_vars_1, f1);

		//restore original
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		
		for(i = 0; i < JACOBIAN_EQUALITIES_NUM_ROWS; i++)
		{
			Jacobian_Equalities[i][j] = (f0[i] - f1[i]) / (2*ESP_DIFFERENTIATION);
		}
	}
}



void Compute_Jacobian_Inequalities(const struct_ip_vars &s_ip_vars, double Jacobian_Inequalities[JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS])
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0[NUM_INEQUALITY_CONSTRAINTS], f1[NUM_INEQUALITY_CONSTRAINTS];

	s_ip_vars_1 = s_ip_vars;
	for(j = 0; j < JACOBIAN_INEQUALITIES_NUM_COLS; j++)
	{
		//f0
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		Compute_inequality_constraints(s_ip_vars_1, f0);

		//f1
		s_ip_vars_1.Optimization_variables[j] -= (2*ESP_DIFFERENTIATION);
		Compute_inequality_constraints(s_ip_vars_1, f1);

		//restore original
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		
		for(i = 0; i < JACOBIAN_INEQUALITIES_NUM_ROWS; i++)
		{
			Jacobian_Inequalities[i][j] = (f0[i] - f1[i]) / (2*ESP_DIFFERENTIATION);
		}
	}
}


void Compute_gradient_objective(const struct_ip_vars &s_ip_vars, double Gradient_Objective[NUM_OPTMIZATION_VARIABLES])
{
	int i,j;
	double f0, f1;
	struct_ip_vars s_ip_vars_1;

	s_ip_vars_1 = s_ip_vars;
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		//f0
		s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION;
		Compute_objective_value(s_ip_vars_1, f0);

		//f1
		s_ip_vars_1.Optimization_variables[i] -= (2*ESP_DIFFERENTIATION);
		Compute_objective_value(s_ip_vars_1, f1);

		//restore original
		s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION;

		Gradient_Objective[i] = (f0 - f1) / (2*ESP_DIFFERENTIATION);
	}
}


void Compute_Hessian_Lagrangian(const struct_ip_vars &s_ip_vars, double Hessian_Lagrangian[HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE])
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0, f1, f2, f3;

	s_ip_vars_1 = s_ip_vars;
	for(i = 0; i < HESSIAN_LAGRANGIAN_SIZE; i++)
	{
		for(j = 0; j < HESSIAN_LAGRANGIAN_SIZE; j++)
		{
		
			//f3
			f3 = Compute_Lagrangian(s_ip_vars_1);

			//f1 (i)
			s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION;
			f1 = Compute_Lagrangian(s_ip_vars_1);

			//f0 (i,j)
			s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
			f0 = Compute_Lagrangian(s_ip_vars_1);

			//f2 (j)
			s_ip_vars_1.Optimization_variables[i] -= ESP_DIFFERENTIATION;
			f2 = Compute_Lagrangian(s_ip_vars_1);

			//restore original 
			s_ip_vars_1.Optimization_variables[j] -= ESP_DIFFERENTIATION;
			
			Hessian_Lagrangian[i][j] = (f0-f1-f2+f3)/(ESP_DIFFERENTIATION * ESP_DIFFERENTIATION);

		}
	}
}

double Compute_Lagrangian(const struct_ip_vars &s_ip_vars)
{
	int i = 0;

	double Lagrangian = 0;
	Compute_objective_value(s_ip_vars, Lagrangian);

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];
	Compute_inequality_constraints(s_ip_vars, inequality_constraints);

	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_equality_constraints(s_ip_vars, equality_constraints);

	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		Lagrangian -= s_ip_vars.Lagrange_multiplier_equality[i]*equality_constraints[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		Lagrangian -= s_ip_vars.Lagrange_multiplier_inequality[i]*(inequality_constraints[i] - s_ip_vars.Slack[i]);
	}

	return Lagrangian;
}

#endif // IP_FORMULATE_CPP
