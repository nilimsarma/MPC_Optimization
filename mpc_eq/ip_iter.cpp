#ifndef IP_ITER_CPP
#define IP_ITER_CPP

#include "ip_iter.hpp"

struct_alpha Compute_alpha_max(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;	
	struct_alpha s_alpha_max;
	double t;
		
	s_alpha_max.alpha_s = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{		
		t =  -s_primal_dual_dir.Vector_Ps[i] / s_ip_vars.Slack[i];
		if ( t > s_alpha_max.alpha_s )	s_alpha_max.alpha_s = t;
	}
	s_alpha_max.alpha_s = TAU / s_alpha_max.alpha_s;
	if(s_alpha_max.alpha_s > 1.0) s_alpha_max.alpha_s = 1.0;
		
	s_alpha_max.alpha_z = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{				
		t = -s_primal_dual_dir.Vector_Pz[i] / s_ip_vars.Lagrange_multiplier_inequality[i];
		if ( t > s_alpha_max.alpha_z )	s_alpha_max.alpha_z = t;
	}
	s_alpha_max.alpha_z = TAU / s_alpha_max.alpha_z;
	if(s_alpha_max.alpha_z > 1.0) s_alpha_max.alpha_z = 1.0;
	
	return s_alpha_max;
}

double Compute_merit_function(const struct_ip_vars &s_ip_vars)
{
	int i,j,k;
	double merit_value;

	merit_value = Compute_objective_value(s_ip_vars);

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value -= s_ip_vars.mu * log(s_ip_vars.Slack[i]);
	}

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];
	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
		 	
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * (inequality_constraints[i] - s_ip_vars.Slack[i]) * (inequality_constraints[i] - s_ip_vars.Slack[i]);
	}

	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_equality_constraints(s_ip_vars, equality_constraints);
	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * (equality_constraints[i]) * (equality_constraints[i]);
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
		s_ip_vars_1.Slack[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Ps[i];
	}

	dir_deriv_merit_func += Compute_merit_function(s_ip_vars_1);
	dir_deriv_merit_func /= ESP_DIFFERENTIATION;

	return dir_deriv_merit_func;
}

struct_alpha Backtrack_alpha(struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir, const struct_alpha &s_alpha_max)
{
	int i,j,k;

	struct_ip_vars s_ip_vars_1;
	struct_alpha s_alpha;
	double merit_function, merit_function_updated, dir_deriv_merit_function;
	
	s_alpha.alpha_s = (s_alpha_max.alpha_s / ALPHA_BACKTRACK_RATIO);
	s_alpha.alpha_z = s_alpha_max.alpha_z;

	while( (dir_deriv_merit_function = Compute_directional_derivative_merit_function(s_ip_vars, s_primal_dual_dir)) >= 0)
	{
		s_ip_vars.nu *= 10;
		print_nu_value(s_ip_vars.nu);
	}
	merit_function = Compute_merit_function(s_ip_vars);

	do
	{
		s_ip_vars_1 = s_ip_vars;
		s_alpha.alpha_s *= ALPHA_BACKTRACK_RATIO;
		
		for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++) 
			s_ip_vars_1.Optimization_variables[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Px[i];
		
		for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++) 
			s_ip_vars_1.Slack[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];

		merit_function_updated = Compute_merit_function(s_ip_vars_1);

	} while( (merit_function + ETA * s_alpha.alpha_s * dir_deriv_merit_function - merit_function_updated) < 0);

	print_alpha_value(s_alpha);
	return s_alpha;
}

void Update_IP_vars(struct_ip_vars &s_ip_vars, const struct_alpha &s_alpha, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;
		
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars.Optimization_variables[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Px[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Slack[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Lagrange_multiplier_inequality[i] += s_alpha.alpha_z * s_primal_dual_dir.Vector_Pz[i];			
	}

	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Lagrange_multiplier_equality[i] += s_alpha.alpha_z * s_primal_dual_dir.Vector_Py[i];
	}

	print_ip_vars(s_ip_vars);
}

double Compute_Error(const struct_ip_vars &s_ip_vars, const double mu)
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
		err[1] += fabs(s_ip_vars.Lagrange_multiplier_inequality[i] * s_ip_vars.Slack[i] - mu);
	}

//err[2]

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];
	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
	
	err[2] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[2] += fabs(inequality_constraints[i] - s_ip_vars.Slack[i]);
	}

//err[3]

	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_equality_constraints(s_ip_vars, equality_constraints);

	err[3] = 0.0;	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		err[3] += fabs(equality_constraints[i]);
	}

//total

	error_max = 0.0;
	for(i = 0; i < NUM_ERROR_TERMS; i++)
	{
		if (err[i] > error_max)	{error_max = err[i]; indx = i;}
	}

	print_error_value(error_max, mu);

	return error_max;
}

void Compute_mu(struct_ip_vars &s_ip_vars)
{
	int i;
	
	s_ip_vars.mu = 0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS;i++)
		s_ip_vars.mu += s_ip_vars.Slack[i] * s_ip_vars.Lagrange_multiplier_inequality[i];
	
	s_ip_vars.mu = SIGMA_MU * (s_ip_vars.mu / NUM_INEQUALITY_CONSTRAINTS);
}

void ip_init(struct_ip_vars &s_ip_vars)
{
	int i,j,k;

	// ++ select intial points ++

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
		s_ip_vars.Optimization_variables[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Slack[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Lagrange_multiplier_inequality[i] = (rand() % 1000) / 1000.0;

#if (USE_EQUALITY_CONSTRAINTS == 1)	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Lagrange_multiplier_equality[i] = (rand() % 1000) / 1000.0;	
#endif

	print_ip_vars(s_ip_vars);
	
	// -- select intial points --

	//set parameter values

	s_ip_vars.nu = NU;
	print_nu_value(s_ip_vars.nu);
	
	Compute_mu(s_ip_vars);
	print_mu_value(s_ip_vars.mu);
}

int main(int argc, char** argv)
{
	struct_ip_vars s_ip_vars;
	struct_primal_dual_direction s_primal_dual_dir;
	struct_alpha s_alpha_max, s_alpha;

	print_problem_parameters();

	ip_init(s_ip_vars);

	int i,j,k;

	
	//start iterations	
	k = 0;
	while( Compute_Error(s_ip_vars, 0) > ERROR_TOL_TOTAL )
	{
		j = 0;
		while( Compute_Error(s_ip_vars, s_ip_vars.mu) > (s_ip_vars.mu*10) )
		{
			j++; k++;
			
			//compute primal dual directions
			Compute_primal_dual_direction (s_ip_vars, s_primal_dual_dir);
			
			//compute alpha_max
			s_alpha_max  = Compute_alpha_max(s_ip_vars, s_primal_dual_dir);
			
			//compute alpha
			s_alpha = Backtrack_alpha(s_ip_vars, s_primal_dual_dir, s_alpha_max);

			//update IP variables
			Update_IP_vars(s_ip_vars, s_alpha, s_primal_dual_dir);
		}

		//update mu

		if(j == 0)	s_ip_vars.mu *= SIGMA_MU;
		else	Compute_mu(s_ip_vars);		
		print_mu_value(s_ip_vars.mu);
	}

	printf("\nIterations = %d\n", k);
}

#endif //IP_ITER_CPP
