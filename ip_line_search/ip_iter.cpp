#ifndef IP_ITER_CPP
#define IP_ITER_CPP

#include "ip_iter.hpp"

double Compute_alpha_max(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;	
	double alpha_max;
	double t;
		
	alpha_max = 0.0;		
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		t = ( -s_primal_dual_dir.Vector_Pw[i] / s_ip_vars.W[i] );
		if( t > alpha_max )	alpha_max = t;

		t = ( -s_primal_dual_dir.Vector_Py[i] / s_ip_vars.Lagrange_multiplier_inequality[i] );
		if( t > alpha_max )	alpha_max = t;
	}

	alpha_max = TAU / alpha_max;

	print_alpha_max_value(alpha_max);
	return alpha_max;
}

double Compute_merit_function(const struct_ip_vars &s_ip_vars)
{
	int i,j,k;
	double merit_value;

	merit_value = Compute_objective_value(s_ip_vars);

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value -= s_ip_vars.mu * log(s_ip_vars.W[i]);
	}

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];

	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
	 
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * (inequality_constraints[i] - s_ip_vars.W[i]) * (inequality_constraints[i] - s_ip_vars.W[i]);
	}

	return merit_value;
}

double Compute_directional_derivative_merit_function(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j;

	double dir_deriv_merit_func = 0;
	struct_ip_vars s_ip_vars_1;

	s_ip_vars_1 = s_ip_vars;

	dir_deriv_merit_func = -Compute_merit_function(s_ip_vars);
	
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Px[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars_1.W[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Pw[i];
	}

	dir_deriv_merit_func += Compute_merit_function(s_ip_vars_1);
	dir_deriv_merit_func /= ESP_DIFFERENTIATION;

	return dir_deriv_merit_func;
}

double Backtrack_alpha(struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir, const double &alpha_max)
{
	int i,j,k;

	struct_ip_vars s_ip_vars_1;
	double alpha;
	double merit_function, merit_function_updated, dir_deriv_merit_function;
	
	alpha = alpha_max / ALPHA_BACKTRACK_RATIO;
	
	while( (dir_deriv_merit_function = Compute_directional_derivative_merit_function(s_ip_vars, s_primal_dual_dir)) >= 0)
	{
		s_ip_vars.nu *= 10;	
		print_nu_value(s_ip_vars.nu);
	}
	merit_function = Compute_merit_function(s_ip_vars);

	do
	{
		s_ip_vars_1 = s_ip_vars;
		alpha *= ALPHA_BACKTRACK_RATIO;
		
		for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++) 
			s_ip_vars_1.Optimization_variables[i] += alpha * s_primal_dual_dir.Vector_Px[i];
		
		for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++) 
			s_ip_vars_1.W[i] += alpha * s_primal_dual_dir.Vector_Pw[i];

		merit_function_updated = Compute_merit_function(s_ip_vars_1);

	} while( (merit_function + ETA * alpha * dir_deriv_merit_function - merit_function_updated) < 0);
	
	print_alpha_value(alpha);
	return alpha;
}

void Update_IP_vars(struct_ip_vars &s_ip_vars, const double &alpha, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;
		
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars.Optimization_variables[i] += alpha * s_primal_dual_dir.Vector_Px[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.W[i] += alpha * s_primal_dual_dir.Vector_Pw[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Lagrange_multiplier_inequality[i] += alpha * s_primal_dual_dir.Vector_Py[i];			
	}
	
	print_ip_vars(s_ip_vars);
}

double Compute_Error(const struct_ip_vars &s_ip_vars, double mu)
{
	int i, indx;
	double err[NUM_ERROR_TERMS], error_max;

//err[0]

	double Gradient_Lagrangian[NUM_OPTMIZATION_VARIABLES];
	Compute_Gradient_Lagrangian(s_ip_vars, Gradient_Lagrangian);

	err[0] = 0.0;
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		err[0] += fabs(Gradient_Lagrangian[i]);
	}
	
//err[1]	

	err[1] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[1] += fabs(s_ip_vars.Lagrange_multiplier_inequality[i] * s_ip_vars.W[i] - mu);
	}
	
//err[2]
	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];
	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
	
	err[2] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[2] += fabs(inequality_constraints[i] - s_ip_vars.W[i]);
	}

//total

	error_max = 0;
	for(i = 0; i < NUM_ERROR_TERMS; i++)
	{
		if ( err[i] > error_max)	{error_max = err[i]; indx = i;}
	}

	print_error_value(error_max, mu);

	return error_max;
}

void ip_init(struct_ip_vars &s_ip_vars)
{
	int i,j,k;

	// ++ select intial points ++

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
		s_ip_vars.Optimization_variables[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.W[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Lagrange_multiplier_inequality[i] = (rand() % 1000) / 1000.0;

	print_ip_vars(s_ip_vars);
	// -- select intial points --

	//set parameter values

	s_ip_vars.nu = NU;	
	print_nu_value(s_ip_vars.nu);

	s_ip_vars.mu = 0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.mu += s_ip_vars.W[i] * s_ip_vars.Lagrange_multiplier_inequality[i];
	}
	s_ip_vars.mu *= (SIGMA_MU / NUM_INEQUALITY_CONSTRAINTS);
	print_mu_value(s_ip_vars.mu);
}

int main(int argc, char** argv)
{
	int i,j,k;
	struct_ip_vars s_ip_vars;
	struct_primal_dual_direction s_primal_dual_dir;
	double alpha_max, alpha;

	print_problem_parameters();
	ip_init(s_ip_vars);
	
	//start iterations	
	j = 0;
	while(Compute_Error(s_ip_vars, 0) > ERROR_TOL_TOTAL)
	{
		k = 0;
		while(Compute_Error(s_ip_vars, s_ip_vars.mu) > (s_ip_vars.mu)*10)
		{	
			k++; j++;
			//compute primal dual directions
			Compute_primal_dual_direction (s_ip_vars, s_primal_dual_dir);
			
			//compute alpha_max
			alpha_max  = Compute_alpha_max(s_ip_vars, s_primal_dual_dir);
			
			//compute alpha
			alpha = Backtrack_alpha(s_ip_vars, s_primal_dual_dir, alpha_max);

			//update IP variables
			Update_IP_vars(s_ip_vars, alpha, s_primal_dual_dir);
		}

		//update mu

		if(k == 0)
		{
			s_ip_vars.mu *= SIGMA_MU;
		}
		else
		{
			s_ip_vars.mu = 0;
			for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
			{
				s_ip_vars.mu += s_ip_vars.W[i] * s_ip_vars.Lagrange_multiplier_inequality[i];
			}
			s_ip_vars.mu *= (SIGMA_MU / NUM_INEQUALITY_CONSTRAINTS);
		}		
		print_mu_value(s_ip_vars.mu);
	}

	printf("\n\nIterations = %d\n", j);
}

#endif //IP_ITER_CPP
