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
		double t = ( (-TAU) * s_ip_vars.s[i]) / s_primal_dual_dir.Vector_Ps[i];
		if ( t > s_alpha_max.alpha_s )	s_alpha_max.alpha_s = t;
	}

	if(s_alpha_max.alpha_s > 1) s_alpha_max.alpha_s = 1.0;

	//alpha_z_max
	
	s_alpha_max.alpha_z = 0;

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		double t = ( (-TAU) * s_ip_vars.Lagrange_multiplier_inequality[i]) / s_primal_dual_dir.Vector_Pz[i];
		if ( t > s_alpha_max.alpha_z )	s_alpha_max.alpha_z = t;
	}

	if(s_alpha_max.alpha_z > 1) s_alpha_max.alpha_z = 1.0;

	return s_alpha_max;
}

#endif //IP_ITER_CPP
