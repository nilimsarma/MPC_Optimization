#ifndef IP_ITER_CPP
#define IP_ITER_CPP

struct_alpha Compute_alpha_max(struct_ip_vars s_ip_vars, struct_primal_dual_direction s_primal_dual_dir)
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

double Compute_merit_function(struct_ip_vars s_ip_vars)
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

double Compute_directional_derivative_merit_function(struct_ip_vars s_ip_vars, struct_primal_dual_direction s_primal_dual_dir)
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

struct_alpha alpha_backtrack(struct_ip_vars s_ip_vars, struct_primal_dual_direction s_primal_dual_dir, struct_alpha s_alpha_max)
{
	struct_ip_vars s_ip_vars_1;
	struct_alpha s_alpha;
	double merit_function, merit_function_updated, dir_deriv_merit_function;
	int i,j,k;
	
	s_alpha = s_alpha_max;

	while(1)
	{
		s_ip_vars_1 = s_ip_vars;

		merit_function = Compute_merit_function(s_ip_vars_1);
		dir_deriv_merit_function = Compute_directional_derivative_merit_function(s_ip_vars_1, s_primal_dual_dir);
		
		for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
		{
			s_ip_vars_1.Optimization_variables[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Px[i];
		}

		for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
		{
			s_ip_vars_1.S[i] += s_alpha.alpha_s * s_primal_dual_dir.Vector_Ps[i];
		}

		merit_function_updated = Compute_merit_function(s_ip_vars_1);

		if( merit_function + s_ip_vars_1.eta * s_alpha.alpha_s * dir_deriv_merit_function - merit_function_updated < 0)
		{
			s_alpha.alpha_s -= 0.1;
		}
		else 
		{
			return s_alpha;
		}
	}
	
}

struct_ip_vars Update_IP_vars(struct_ip_vars s_ip_vars, struct_alpha s_alpha, struct_primal_dual_direction s_primal_dual_dir)
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
		s_ip_vars.Lagrange_multiplier_equality[i] += s_alpha.alpha_z * s_primal_dual_dir.Py[i];
	}

	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)
	{
		s_ip_vars.Lagrange_multiplier_inequality[i] += s_alpha.alpha_z * s_primal_dual_dir.Pz[i];			
	}

	return s_ip_vars;
}

double Compute_error(struct_ip_vars s_ip_vars, double mu)
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

	while(Compute_error(s_ip_vars, 0) > ERROR_TOL_TOTAL)
	{
		while(Compute_error(s_ip_vars, s_primal_dual_dir.mu) > ERROR_TOL_MU)
		{
			//compute primal dual directions
			s_primal_dual_dir = Compute_primal_dual_direction (s_ip_vars);
			
			//compute alpha_max
			s_alpha_max  = Compute_alpha_max(s_ip_vars, s_primal_dual_dir);
			
			//compute alpha
			s_alpha = alpha_backtrack(s_ip_vars, s_primal_dual_dir, s_alpha_max);
		}

		//update mu
		s_primal_dual_dir.mu *= SIGMA_MU;
	}
}

#endif //IP_ITER_CPP
