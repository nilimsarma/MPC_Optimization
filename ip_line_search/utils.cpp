#ifndef UTILS_CPP
#define UTILS_CPP

#include "utils.hpp"

void print_matrix(double* A, int num_rows, int num_cols)
{
	int i,j;

	printf("\n");
	for(i = 0; i < num_rows; i++)
	{
		for(j = 0; j < num_cols; j++)
		{
			printf("%f, ", A[i*num_cols + j]);
		}
		printf(";\n");
	}
}

void print_vector(double* A, int vec_size)
{
	int i,j;

	printf("\n");
	for(i = 0; i < vec_size; i++)
	{
		printf("%f, ", A[i]);
	
	}
	printf("\n");
}

void print_ip_vars(const struct_ip_vars &s_ip_vars)
{
	int i,j,k;

#if (LOG_OPTIMIZATION_VARIABLES == 1)
	printf("\n***Optimization Variables***\n");
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		printf("\n%f", s_ip_vars.Optimization_variables[i]);
	}
	printf("\n");
#endif

#if (LOG_SLACK_VARIABLES == 1)
	printf("\n***Slack Varibles***\n");
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		printf("\n%f", s_ip_vars.W[i]);
	}
	printf("\n");
#endif

#if (LOG_INEQUALITY_CONSTRAINTS == 1)
	printf("\n***Lagrange Multiplier Inequality***\n");
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		printf("\n%f", s_ip_vars.Lagrange_multiplier_inequality[i]);
	}
	printf("\n");	
#endif

}


void print_problem_parameters(void)
{
	//Optimization problem formulation parameters

	printf("\n**************\n");

	printf("\nNUM_OPTMIZATION_VARIABLES = %d", NUM_OPTMIZATION_VARIABLES);
	printf("\nNUM_INEQUALITY_CONSTRAINTS = %d", NUM_INEQUALITY_CONSTRAINTS);

	printf("\n**************\n");
}

void print_mu_value(const double &mu)
{
#if (LOG_MU_VALUE == 1)
	printf("\nmu = %f\n", mu);
#endif	
}

void print_nu_value(const double &nu)
{
#if (LOG_NU_VALUE == 1)
	printf("\nnu = %e\n", nu);
#endif	
}

void print_alpha_value(const double &alpha)
{
#if (LOG_ALPHA_VALUE == 1)
	printf("\nalpha = %f\n",alpha);
#endif
}

void print_alpha_max_value(const double &alpha_max)
{
#if (LOG_ALPHA_MAX_VALUE == 1)
	printf("\nalpha_max = %f\n",alpha_max);
#endif
}

void print_error_value(const double &error, const double &mu)
{
#if (LOG_ERROR_VALUE == 1)
	if(mu == 0)	printf("\n***Error*** = %e\n", error);
	else printf("\nError = %e\n", error);
#endif	
}

#endif //UTILS_CPP
