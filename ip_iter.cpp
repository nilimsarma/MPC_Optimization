#ifndef IP_ITER_CPP
#define IP_ITER_CPP

#include "ip_iter.hpp"

struct_alpha Compute_alpha_max(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j,k;	
	struct_alpha s_alpha_max;

	//alpha_s_max
	  
	s_alpha_max.alpha_s = 0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		double t = ( (-TAU) * s_ip_vars.S[i]) / s_primal_dual_dir.Vector_Ps[i];
		if ( t > s_alpha_max.alpha_s )	s_alpha_max.alpha_s = t;
	}

	if(s_alpha_max.alpha_s > 1.0) s_alpha_max.alpha_s = 1.0;

	//alpha_z_max
	
	s_alpha_max.alpha_z = 0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		double t = ( (-TAU) * s_ip_vars.Lagrange_multiplier_inequality[i]) / s_primal_dual_dir.Vector_Pz[i];
		if ( t > s_alpha_max.alpha_z )	s_alpha_max.alpha_z = t;
	}

	if(s_alpha_max.alpha_z > 1.0) s_alpha_max.alpha_z = 1.0;

	return s_alpha_max;
}

double Compute_merit_function(const struct_ip_vars &s_ip_vars)
{
	int i,j,k;
	double merit_value;

	merit_value += Compute_objective_value(s_ip_vars);

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value -= s_ip_vars.mu * log(s_ip_vars.s[i]);
	}

	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS] = Compute_inequality_constraints(s_ip_vars);
	double equality_constraints[NUM_EQUALITY_CONSTRAINTS] = Compute_equality_constraints(s_ip_vars);
	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * fabs(equality_constraints[i]);
	}
	
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		merit_value += s_ip_vars.nu * fabs(inequality_constraints[i] - s_ip_vars.s[i]);
	}

	return merit_value;
}

double Compute_directional_derivative_merit_function(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir)
{
	int i,j;

	double dir_deriv_merit_func = 0;

	dir_deriv_merit_func = -Compute_merit_function(s_ip_vars);
	
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars.Optimization_variables[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Px[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		s_ip_vars.S[i] += ESP_DIFFERENTIATION * s_primal_dual_dir.Vector_Ps[i];
	}

	dir_deriv_merit_func += Compute_merit_function(s_ip_vars);
	dir_deriv_merit_func /= ESP_DIFFERENTIATION;

	return dir_deriv_merit_func;
}

struct_alpha Backtrack_alpha(const struct_ip_vars &s_ip_vars, const struct_primal_dual_direction &s_primal_dual_dir, const struct_alpha &s_alpha_max)
{
	int i,j,k;

	struct_ip_vars s_ip_vars_1;
	struct_alpha s_alpha;
	double merit_function, merit_function_updated, dir_deriv_merit_function;
	
	s_alpha = s_alpha_max;

	merit_function = Compute_merit_function(s_ip_vars);
	dir_deriv_merit_function = Compute_directional_derivative_merit_function(s_ip_vars, s_primal_dual_dir);

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
			s_alpha.alpha_s -= ALPHA_BACKTRACK_VAL;

			for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
				s_ip_vars_1.Optimization_variables[i] -= ALPHA_BACKTRACK_VAL * s_primal_dual_dir.Vector_Px[i];
			
			for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
				s_ip_vars_1.S[i] -= ALPHA_BACKTRACK_VAL * s_primal_dual_dir.Vector_Ps[i];
		}
		else 
		{
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

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars.Lagrange_multiplier_inequality[i] += s_alpha.alpha_z * s_primal_dual_dir.Vector_Pz[i];			
	}

	return s_ip_vars;
}

double Compute_Error(const struct_ip_vars &s_ip_vars, double mu)
{
	int i;
	double err[NUM_ERROR_TERMS], error_max;

	err[0] = fabs(Compute_Gradient_Lagrangian(s_ip_vars));

	err[1] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[1] += fabs(s_ip_vars.Lagrange_multiplier_inequality[i] - mu / s_ip_vars.S[i]);
	}

	double equality_constraints[NUM_EQUALITY_CONSTRAINTS] = Compute_equality_constraints(s_ip_vars);

	err[2] = 0.0;	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		err[2] += fabs(equality_constraints[i]);
	}
	
	double inequality_constraints[NUM_EQUALITY_CONSTRAINTS] = Compute_inequality_constraints(s_ip_vars);
	
	err[3] = 0.0;
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		err[3] += fabs(inequality_constraints[i] - s_ip_vars.S[i]);
	}

	error_max = err[0];
	for(i = 1; i < NUM_ERROR_TERMS; i++)
	{
		if ( err[i] > error_max)	error_max = err[i];
	}

	return error_max;
}

int main(int argc, char** argv)
{
	struct_ip_vars s_ip_vars;
	struct_primal_dual_direction s_primal_dual_dir;
	struct_alpha s_alpha_max, s_alpha;
	
	//select intial points

	while(Compute_Error(s_ip_vars, 0) > ERROR_TOL_TOTAL)
	{
		while(Compute_Error(s_ip_vars, s_primal_dual_dir.mu) > ERROR_TOL_MU)
		{
			//compute primal dual directions
			Compute_primal_dual_direction (s_ip_vars, s_primal_dual_dir);
			
			//compute alpha_max
			s_alpha_max  = Compute_alpha_max(s_ip_vars, s_primal_dual_dir);
			
			//compute alpha
			s_alpha = Backtrack_alpha(s_ip_vars, s_primal_dual_dir, s_alpha_max);

			//update IP variables
			Update_IP_vars(s_ip_vars, s_alpha, s_primal_dual_dir)
		}

		//update mu
		s_primal_dual_dir.mu *= SIGMA_MU;
	}
}

#endif //IP_ITER_CPP
