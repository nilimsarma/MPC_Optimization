#ifndef IP_ITER_CPP
#define IP_ITER_CPP

void Compute_Lagrange_Multipliers(struct_ip_vars &s_ip_vars)
{
	int i,j,k;
	double c;
	
	double Jacobian_Equalities[JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS];
	Compute_Jacobian_Equalities(s_ip_vars, Jacobian_Equalities);
	
	double Jacobian_Inequalities[JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS];
	Compute_Jacobian_Inequalities(s_ip_vars, Jacobian_Inequalities);

/* ++ Create Matrix A_Cap ++ */

	double Matrix_A_Cap[MAT_A_CAP_NUM_ROWS][MAT_A_CAP_NUM_COLS];

	//init to all zeros
	for(i = 0; i < MAT_A_CAP_NUM_ROWS; i++)
		for(j = 0; j < MAT_A_CAP_NUM_COLS; j++)
			Matrix_A_Cap[i][j] = 0.0;

	//Jacobian Equality
	for(i = 0; i < JACOBIAN_EQUALITIES_NUM_ROWS; i++)
	{	
		for(j = 0; j < JACOBIAN_EQUALITIES_NUM_COLS; j++)
		{
			Matrix_A_Cap[i][j] = Jacobian_Equalities[i][j];
		}
	}

	//Jacobian Inequality
	for(i = 0; i < JACOBIAN_INEQUALITIES_NUM_ROWS; i++)
	{	
		for(j = 0; j < JACOBIAN_INEQUALITIES_NUM_COLS; j++)
		{
			Matrix_A_Cap[JACOBIAN_EQUALITIES_NUM_ROWS + i][j] = Jacobian_Inequalities[i][j];
		}
	}

	//Slack Matrix
	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		Matrix_A_Cap[JACOBIAN_EQUALITIES_NUM_ROWS + i][JACOBIAN_EQUALITIES_NUM_COLS + i] = -s_ip_vars.Slack[i];
	}

/* -- Create Matrix A_Cap -- */

/* ++ Create Vector b  ++ */

	double Gradient_Objective[NUM_OPTMIZATION_VARIABLES];
	Compute_gradient_objective(s_ip_vars, Gradient_Objective);
		
	double Vec_b[MAT_A_CAP_NUM_COLS];
	
	for(i = 0; i < NUM_OPTMIZATION_VARIABLES; i++)	
	{
		Vec_b[i] = Gradient_Objective[i];
	}

	for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++)
	{
		Vec_b[NUM_OPTMIZATION_VARIABLES + i] = -s_ip_vars.mu;
	}
	
/* -- Create Vector b  -- */

/* ++ Solve Linear System ++ */

	//a = A_Cap * trans(A_Cap)
	double a[MAT_A_CAP_NUM_ROWS * MAT_A_CAP_NUM_ROWS];
	
	for(i = 0; i < MAT_A_CAP_NUM_ROWS; i++)
	{
		for(j = 0; j < MAT_A_CAP_NUM_ROWS; j++)
		{
			c = 0;
			for(k = 0; k < MAT_A_CAP_NUM_COLS; k++)
			{
				c += Matrix_A_Cap[i][k] * Matrix_A_Cap[j][k];
			}
			
			//column major
			a[j*MAT_A_CAP_NUM_ROWS + i] = c;
		}
	}
	
	//b = A_Cap * Vec_b
	double b[MAT_A_CAP_NUM_ROWS];
	
	for(i = 0; i < MAT_A_CAP_NUM_ROWS; i++)
	{
		c = 0;
		for(j = 0; j < MAT_A_CAP_NUM_COLS; j++)
		{
			c += Matrix_A_Cap[i][j]*Vec_b[j];
		}
		b[i] = c;
	}

	
	lapack_int n, nrhs, lda, ldb, info;	
	n = MAT_A_CAP_NUM_ROWS;
	nrhs = 1;
	lda = MAT_A_CAP_NUM_ROWS;
	ldb = MAT_A_CAP_NUM_ROWS;
	lapack_int ipiv[MAT_A_CAP_NUM_ROWS];
	
//	lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, lapack_int* ipiv, double* b, lapack_int ldb );
	info = LAPACKE_dgesv( LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );

	if(info == 0) 
	{
		for(i = 0; i < NUM_EQUALITY_CONSTRAINTS; i++) 
		{
			s_ip_vars.Lagrange_multiplier_equality[i] = b[i];
		}
		
		for(i = 0; i < NUM_INEQUALITY_CONSTRAINTS; i++) 
		{
			if(b[i] > 0)	
				s_ip_vars.Lagrange_multiplier_equality[i] = b[i];
			
			else if( (s_ip_vars.mu / s_ip_vars.Slack[i]) < 0.001)	
				s_ip_vars.Lagrange_multiplier_equality[i] = (s_ip_vars.mu / s_ip_vars.Slack[i]);

			else
				s_ip_vars.Lagrange_multiplier_equality[i] = 0.001;
		}
	}
	else 
	{
		printf("\nLapack error!!!\n");
		printf("\ninfo = %d\n", info);
		exit(1);
	}

/* -- Solve Linear System -- */

}

#endif // IP_ITER_CPP

