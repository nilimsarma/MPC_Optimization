#ifndef IP_PRIMAL_DUAL_DIR_CPP
#define IP_PRIMAL_DUAL_DIR_CPP

//Solve Linear System Ax = b
		
struct_primal_dual_direction Compute_primal_dual_direction (struct_ip_vars s_ip_vars)
{
	int i,j,k,m;
		
	double Hessian_Lagrangian [HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE];
	double Jacobian_Equalities [JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS];
	double Jacobian_Inequalities [JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS];
	double Diag_Matrix_Sigma [DIAG_MATRIX_SIGMA_SIZE];

	Hessian_Lagrangian = Compute_Hessian_Lagrangian(s_ip_vars);
	Jacobian_Equalities	= Compute_Jacobian_Equalities(s_ip_vars);
	Jacobian_Inequalities = Compute_Jacobian_Inequalities(s_ip_vars);
	Diag_Matrix_Sigma = Compute_Diag_Matrix_Sigma(s_ip_vars);

	//Create Sq_Matrix_A
	double Sq_Matrix_A [SQ_MATRIX_A_SIZE][SQ_MATRIX_A_SIZE];

	for(i = 0; i < SQ_MATRIX_A_SIZE; i++)
	for(j = 0; j < SQ_MATRIX_A_SIZE; j++)
	{
		Sq_Matrix_A[i][j] = 0.0;
	}
	
	for(i = 0; i < HESSIAN_LAGRANGIAN_SIZE; i++)
	for(j = 0; j < HESSIAN_LAGRANGIAN_SIZE; j++)
	{
		Sq_Matrix_A[i][j] = Hessian_Lagrangian[i][j];
	}

	for(i = 0, k = ( HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE); i < JACOBIAN_EQUALITIES_NUM_ROWS ; i++, k++)
	for(j = 0; m = 0; j < JACOBIAN_EQUALITIES_NUM_COLS; j++, m++)
	{
		Sq_Matrix_A[k][m] = Jacobian_Equalities[i][j];
	}

	for(i = 0, k = ( HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE + JACOBIAN_EQUALITIES_NUM_ROWS ); i < JACOBIAN_INEQUALITIES_NUM_ROWS ; i++, k++)
	for(j = 0; m = 0; j < JACOBIAN_INEQUALITIES_NUM_COLS; j++, m++)
	{
		Sq_Matrix_A[k][m] = Jacobian_Inequalities[i][j];
	}

	for(i = HESSIAN_LAGRANGIAN_SIZE, j = 0; i < ( HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE ); i++, j++ )
	{
		Sq_Matrix_A[i][i] = Diag_Matrix_Sigma[j];
	}

	for(i = ( HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE + JACOBIAN_EQUALITIES_NUM_ROWS ), j = JACOBIAN_INEQUALITIES_NUM_COLS; i < SQ_MATRIX_A_SIZE; i++, j++)
	{
		Sq_Matrix_A[i][j] = -1.0;
	}

	//copy upper triangular
	for(i = 0; i < SQ_MATRIX_A_SIZE; i++)
	for(j = 0; j < i; j++)
	{
		Sq_Matrix_A[j][i] = Sq_Matrix_A[i][j];
	}

	double Vector_b0 [VECTOR_SIZE_b0];	//Jacobian_Lagrangian
	double Vector_b1 [VECTOR_SIZE_b1];
	double Vector_b2 [VECTOR_SIZE_b2];
	double Vector_b3 [VECTOR_SIZE_b3];

	Vector_b0 = Compute_Gradient_Lagrangian(s_ip_vars);
	Vector_b1 = Compute_vector_b1(s_ip_vars);
	Vector_b2 = Compute_vector_b2(s_ip_vars);
	Vector_b3 = Compute_vector_b3(s_ip_vars);

	//Create Vector b
	double Vector_b	[VECTOR_b_SIZE];
	
	j = 0;
	for(i = 0; i < VECTOR_SIZE_b0; i++, j++)	Vector_b[j] = -Vector_b0[i];
	for(i = 0; i < VECTOR_SIZE_b1; i++, j++)	Vector_b[j] = -Vector_b1[i];
	for(i = 0; i < VECTOR_SIZE_b2; i++, j++)	Vector_b[j] = -Vector_b2[i];
	for(i = 0; i < VECTOR_SIZE_b3; i++, j++)	Vector_b[j] = -Vector_b3[i];		

	double Vector_x [VECTOR_x_SIZE];	

/* ++ Solve Linear System ++ */

	double a[SQ_MATRIX_A_SIZE * SQ_MATRIX_A_SIZE];
	double b[VECTOR_b_SIZE];

	//switch to column major
	for (i = 0; i < SQ_MATRIX_A_SIZE; i++)
		for(j = 0; j < SQ_MATRIX_A_SIZE; j++) 
			a[j * SQ_MATRIX_A_SIZE + i] = Sq_Matrix_A[i][j];	

	for(i = 0; i < VECTOR_b_SIZE; i++)	b[i] = Vector_b[i];
		
	int n = SQ_MATRIX_A_SIZE, nrhs = 1, lda = SQ_MATRIX_A_SIZE, ldb = VECTOR_b_SIZE, info;
	int ipiv[SQ_MATRIX_A_SIZE];
	
	dgesv( n, nrhs, a, lda, ipiv, b, ldb, info );

	if(info == 0) 
	{
		for(i = 0; i < VECTOR_x_SIZE; i++) Vector_x[i] = b[i];
	}
	else 
	{
		printf("\nLapack error!!!\n");
		exit(1);
	}

/* -- Solve Linear System -- */

//	double Vector_Px [VECTOR_SIZE_Px];
//	double Vector_Ps [VECTOR_SIZE_Ps];
//	double Vector_Py [VECTOR_SIZE_Py];
//	double Vector_Pz [VECTOR_SIZE_Pz];

	struct_primal_dual_direction s_primal_dual_dir;

	j = 0;
	for(i = 0; i < VECTOR_SIZE_Px; i++, j++)	s_primal_dual_dir.Vector_Px[i] = Vector_x[j];
	for(i = 0; i < VECTOR_SIZE_Ps; i++, j++)	s_primal_dual_dir.Vector_Ps[i] = Vector_x[j];
	for(i = 0; i < VECTOR_SIZE_Py; i++, j++)	s_primal_dual_dir.Vector_Py[i] = -Vector_x[j];
	for(i = 0; i < VECTOR_SIZE_Pz; i++, j++)	s_primal_dual_dir.Vector_Pz[i] = -Vector_x[j];		

	return s_primal_dual_dir;
	
}


double[HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE] 
Compute_Hessian_Lagrangian(struct_ip_vars s_ip_vars)
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0, f1, f2, f3;

	double Hessian_Lagrangian[HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE];
		
	for(i = 0; i < HESSIAN_LAGRANGIAN_SIZE; i++)
	{
		for(j = 0; j < HESSIAN_LAGRANGIAN_SIZE; j++)
		{
			//f0
			s_ip_vars_1 = s_ip_vars;			
			s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION;
			s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
			f0 = Compute_Lagrangian(s_ip_vars_1);

			//f1
			s_ip_vars_1 = s_ip_vars;			
			s_ip_vars_1.Optimization_variables[i] += ESP_DIFFERENTIATION;
			f1 = Compute_Lagrangian(s_ip_vars_1);

			//f2
			s_ip_vars_1 = s_ip_vars;			
			s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
			f2 = Compute_Lagrangian(s_ip_vars_1);

			//f3
			s_ip_vars_1 = s_ip_vars;			
			f3 = Compute_Lagrangian(s_ip_vars_1);

			Hessian_Lagrangian[i][j] = (f0-f1-f2+f3)/ESP_DIFFERENTIATION;
			
		}
	}

	return Hessian_Lagrangian;
}

double [JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS]
Compute_Jacobian_Equalities(struct_ip_vars s_ip_vars)
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0, f1;

	double Jacobian_Equalities[JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS];
	double f0[NUM_EQUALITY_CONSTRAINTS], f1[NUM_EQUALITY_CONSTRAINTS];
	
	for(j = 0; j < JACOBIAN_EQUALITIES_NUM_COLS; j++)
	{
		//f0
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		f0 = Compute_equality_constraints(s_ip_vars_1);

		//f1
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] -= ESP_DIFFERENTIATION;
		f1 = Compute_equality_constraints(s_ip_vars_1);

		for(i = 0; i < JACOBIAN_EQUALITIES_NUM_ROWS; i++)
		{
			Jacobian_Equalities[i][j] = (f0[i] - f1[i])/(2*ESP_DIFFERENTIATION);
		}
	}

	return Jacobian_Equalities;
	
}

double [JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS]
Compute_Jacobian_Inequalities(struct_ip_vars s_ip_vars)
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0, f1;

	double Jacobian_Inequalities[JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS];
	double f0[NUM_INEQUALITY_CONSTRAINTS], f1[NUM_INEQUALITY_CONSTRAINTS];
	
	for(j = 0; j < JACOBIAN_INEQUALITIES_NUM_COLS; j++)
	{
		//f0
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		f0 = Compute_inequality_constraints(s_ip_vars_1);

		//f1
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] -= ESP_DIFFERENTIATION;
		f1 = Compute_inequality_constraints(s_ip_vars_1);

		for(i = 0; i < JACOBIAN_INEQUALITIES_NUM_ROWS; i++)
		{
			Jacobian_Inequalities[i][j] = (f0[i] - f1[i]) / (2*ESP_DIFFERENTIATION);
		}
	}

	return Jacobian_Inequalities;
	
}


double [DIAG_MATRIX_SIGMA_SIZE]
Compute_Diag_Matrix_Sigma(struct_ip_vars s_ip_vars)
{
	int i;
	double Diag_Matrix_Sigma [DIAG_MATRIX_SIGMA_SIZE];
	
	for(i = 0; i < DIAG_MATRIX_SIGMA_SIZE; i++)
	{
		Diag_Matrix_Sigma[i] = s_ip_vars.Lagrange_multiplier_inequality[i] / s_ip_vars.S[i];
	}
}

double [NUM_OPTMIZATION_VARIABLES] 
Compute_Gradient_Lagrangian(struct_ip_vars s_ip_vars)
{
	int i,j;
	struct_ip_vars s_ip_vars_1;
	double f0, f1;

	double Gradient_Lagrangian[NUM_OPTMIZATION_VARIABLES];
	double f0, f1;
	
	for(j = 0; j < NUM_OPTMIZATION_VARIABLES; j++)
	{
		//f0
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] += ESP_DIFFERENTIATION;
		f0 = Compute_Lagrangian(s_ip_vars_1);

		//f1
		s_ip_vars_1 = s_ip_vars;
		s_ip_vars_1.Optimization_variables[j] -= ESP_DIFFERENTIATION;
		f1 = Compute_Lagrangian(s_ip_vars_1);
		
		Gradient_Lagrangian[j] = (f0 - f1) / (2*ESP_DIFFERENTIATION);
	}

	return Gradient_Lagrangian;
	
}


double [VECTOR_SIZE_b1]
Compute_vector_b1(struct_ip_vars s_ip_vars)
{
	int i;
	double vector_b1[VECTOR_SIZE_b1];

	for(i = 0; i < VECTOR_SIZE_b1; i++)
	{
		vector_b1[i] = s_ip_vars.Lagrange_multiplier_inequality[i] - s_ip_vars.mu / s_ip_vars.S[i];
	}

	return vector_b1;
}

double [VECTOR_SIZE_b2]
Compute_vector_b2(struct_ip_vars s_ip_vars)
{
	double vector_b2[VECTOR_SIZE_b2];
	vector_b2 = Compute_equality_constraints(s_ip_vars);
	return vector_b2;
}

double [VECTOR_SIZE_b3]
Compute_vector_b3(struct_ip_vars s_ip_vars)
{
	int i;
	
	double vector_b3[VECTOR_SIZE_b3];
	vector_b3 = Compute_inequality_constraints(s_ip_vars);

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		vector_b3[i] -= s_ip_vars.S[i];
	}

	return vector_b3;
}

double 
Compute_Lagrangian(struct_ip_vars s_ip_vars)
{
	double Lagrangian = 0;
	double equality_constraints[NUM_EQUALITY_CONSTRAINTS];
	double inequality_constraints[NUM_INEQUALITY_CONSTRAINTS];

	int i = 0;

	Lagrangian = Compute_objective_value(s_ip_vars);
	equality_constraints = Compute_equality_constraints(s_ip_vars);
	inequality_constraints = Compute_inequality_constraints(s_ip_vars);
	
	for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++)
	{
		Lagrangian += s_ip_vars.Lagrange_multiplier_equality[i]*equality_constraints[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		Lagrangian += s_ip_vars.Lagrange_multiplier_inequality[i]*(inequality_constraints[i] - s_ip_vars.S[i]);
	}

	return Lagrangian;
}


#endif // IP_PRIMAL_DUAL_DIR_CPP