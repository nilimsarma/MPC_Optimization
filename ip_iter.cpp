#ifndef IP_ITER_CPP
#define IP_ITER_CPP

#include "ip_iter.hpp"

struct_alpha Compute_alpha_max(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;	
	struct_alpha s_alpha_max;
	double t;
		
	//alpha_s_max
	  
	s_alpha_max.alpha_s = 1.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		if( s_ip_vars.S[i] < 0 )	printf("\ns_ip_vars.S[%d] = %f", i, s_ip_vars.S[i]);
		
		if(s_primal_dual_dir.Vector_Ps[i] < 0)
		{
			t = ( (-TAU) * s_ip_vars.S[i]) / s_primal_dual_dir.Vector_Ps[i];
			if ( t < s_alpha_max.alpha_s )	s_alpha_max.alpha_s = t;
		}
		//printf("\ns_alpha_max.alpha_s = %f", s_alpha_max.alpha_s);
	}

	//alpha_z_max
	
	s_alpha_max.alpha_z = 1.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		if( s_ip_vars.Lagrange_multiplier_inequality[i] < 0 )	
			printf("\ns_ip_vars.Lagrange_multiplier_inequality[%d] = %f", i, s_ip_vars.Lagrange_multiplier_inequality[i]);
		
		if(s_primal_dual_dir.Vector_Pz[i] < 0)
		{
			t = ( (-TAU) * s_ip_vars.Lagrange_multiplier_inequality[i]) / s_primal_dual_dir.Vector_Pz[i];
			if ( t < s_alpha_max.alpha_z )	s_alpha_max.alpha_z = t;
		}
		//printf("\ns_alpha_max.alpha_z = %f", s_alpha_max.alpha_z);
	}

	return s_alpha_max;
}

double Compute_merit_function(const struct_ip_vars &s_ip_vars)
{
	int i,j,k;
	double merit_value;

	merit_value = Compute_objective_value(s_ip_vars);

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value -= s_ip_vars.mu * log(s_ip_vars.S[i]);
	}

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];
	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];

	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
	Compute_equality_constraints(s_ip_vars, equality_constraints);
	 
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * fabs(equality_constraints[i]);
	}
	
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * fabs(inequality_constraints[i] - s_ip_vars.S[i]);
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
		s_ip_vars_1.S[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Ps[i];
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
	
	s_alpha = s_alpha_max;

	merit_function = Compute_merit_function(s_ip_vars);

	while( (dir_deriv_merit_function = Compute_directional_derivative_merit_function(s_ip_vars, s_primal_dual_dir)) > 0)
	{
		s_ip_vars.nu *= 10;
	}
	
	s_ip_vars_1 = s_ip_vars;

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++) 
		s_ip_vars_1.Optimization_variables[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Px[i];
	
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++) 
		s_ip_vars_1.S[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];

	while(1)
	{
		
		merit_function_updated = Compute_merit_function(s_ip_vars_1);

		if( (merit_function + s_ip_vars_1.eta * s_alpha.alpha_s * dir_deriv_merit_function - merit_function_updated) < 0)
		{
			s_alpha.alpha_s *= ALPHA_BACKTRACK_RATIO;

			for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
				s_ip_vars_1.Optimization_variables[i] -= (1 - ALPHA_BACKTRACK_RATIO) * s_alpha.alpha_s * s_primal_dual_dir.Vector_Px[i];
			
			for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
				s_ip_vars_1.S[i] -= (1 - ALPHA_BACKTRACK_RATIO) * s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];
		}
		else 
		{
			if(s_alpha.alpha_s < 0.0) s_alpha.alpha_s = 0.0;
			if(s_alpha.alpha_z < 0.0) s_alpha.alpha_z = 0.0;
			
			printf("\nalpha_s = %f, alpha_z = %f\n",s_alpha.alpha_s,s_alpha.alpha_z);

			return s_alpha;
		}
	}
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
		s_ip_vars.S[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];
	}

	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Lagrange_multiplier_equality[i] += s_alpha.alpha_z * s_primal_dual_dir.Vector_Py[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.Lagrange_multiplier_inequality[i] += s_alpha.alpha_z * s_primal_dual_dir.Vector_Pz[i];			
	}
}

double Compute_Error(const struct_ip_vars &s_ip_vars, double mu)
{
	int i;
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
		err[1] += fabs(s_ip_vars.Lagrange_multiplier_inequality[i] * s_ip_vars.S[i] - mu);
	}
	
//err[2]

	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_equality_constraints(s_ip_vars, equality_constraints);

	err[2] = 0.0;	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		err[2] += fabs(equality_constraints[i]);
	}

//err[3]

	double inequality_constraints[NUM_EQUALITY_CONSTRAINTS];
	Compute_inequality_constraints(s_ip_vars, inequality_constraints);
	
	err[3] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[3] += fabs(inequality_constraints[i] - s_ip_vars.S[i]);
	}

//total

	error_max = 0;
	for(i = 0; i < NUM_ERROR_TERMS; i++)
	{
		if ( err[i] > error_max)	error_max = err[i];
	}

	if(mu == 0)	printf("\nError = %e", error_max);
	else printf("\n\tError = %e", error_max);

	return error_max;
}

void ip_init(struct_ip_vars &s_ip_vars)
{
	int i,j,k;

	// ++ select intial points ++

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
		s_ip_vars.Optimization_variables[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.S[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Lagrange_multiplier_equality[i] = (rand() % 1000) / 1000.0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		s_ip_vars.Lagrange_multiplier_inequality[i] = (rand() % 1000) / 1000.0;

	// -- select intial points --

	//set parameter values

	s_ip_vars.nu = NU;
	s_ip_vars.eta = ETA;

	
	s_ip_vars.mu = 0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS;i++)
	{
		s_ip_vars.mu += s_ip_vars.S[i] * s_ip_vars.Lagrange_multiplier_inequality[i];
	}
	s_ip_vars.mu /= NUM_INEQUALITY_CONSTRAINTS;
}

void print_problem_parameters(void)
{
	int i,j,k;
	
	//Optimization problem formulation parameters

	printf("\n**************\n");

	printf("\nNUM_OPTMIZATION_VARIABLES = %d", NUM_OPTMIZATION_VARIABLES);
	printf("\nNUM_EQUALITY_CONSTRAINTS = %d", NUM_EQUALITY_CONSTRAINTS);
	printf("\nNUM_INEQUALITY_CONSTRAINTS = %d", NUM_INEQUALITY_CONSTRAINTS);

	printf("\n**************\n");
}

int main(int argc, char** argv)
{
	struct_ip_vars s_ip_vars;
	struct_primal_dual_direction s_primal_dual_dir;
	struct_alpha s_alpha_max, s_alpha;

	ip_init(s_ip_vars);
	print_problem_parameters();

	int i,j ;

	
	//start iterations	
	//i = 0;
	while( Compute_Error(s_ip_vars, 0) > ERROR_TOL_TOTAL )
	{
		//j = 0;
		while( Compute_Error(s_ip_vars, s_ip_vars.mu) > s_ip_vars.mu )
		{
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
		
		//s_ip_vars.mu *= SIGMA_MU;
		
		s_ip_vars.mu = 0;
		
		for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS;i++)
		{
			s_ip_vars.mu += s_ip_vars.S[i] * s_ip_vars.Lagrange_multiplier_inequality[i];
		}
		s_ip_vars.mu /= NUM_INEQUALITY_CONSTRAINTS;
	}
}

#endif //IP_ITER_CPP
