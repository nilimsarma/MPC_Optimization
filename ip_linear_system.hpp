#ifndef IP_LINEAR_SYSTEM_HPP
#define IP_LINEAR_SYSTEM_HPP

#include <stdio.h>
#include <stdlib.h>
#include "mpc_formulate.hpp"

//Solve Linear System Ax = b

#define HESSIAN_LAGRANGIAN_SIZE 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_EQUALITIES_NUM_ROWS	NUM_EQUALITY_CONSTRAINTS
#define JACOBIAN_EQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_INEQUALITIES_NUM_ROWS	NUM_INEQUALITY_CONSTRAINTS
#define JACOBIAN_INEQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

#define DIAG_MATRIX_SIGMA_SIZE	NUM_INEQUALITY_CONSTRAINTS

#define VECTOR_SIZE_Px	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_Ps 	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Py	NUM_EQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Pz	NUM_INEQUALITY_CONSTRAINTS

#define VECTOR_SIZE_b0	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_b1	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_b2	NUM_EQUALITY_CONSTRAINTS
#define VECTOR_SIZE_b3	NUM_INEQUALITY_CONSTRAINTS

#define SQ_MATRIX_A_SIZE (HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE + JACOBIAN_EQUALITIES_NUM_ROWS + JACOBIAN_INEQUALITIES_NUM_ROWS)
#define VECTOR_x_SIZE	 (VECTOR_SIZE_Px + VECTOR_SIZE_Ps + VECTOR_SIZE_Py + VECTOR_SIZE_Pz)
#define VECTOR_b_SIZE	 (VECTOR_SIZE_b0 + VECTOR_SIZE_b1 + VECTOR_SIZE_b2 + VECTOR_SIZE_b3)

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double mu;
	double S[NUM_INEQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_equality[NUM_EQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS];
} struct_ip_vars;


#endif // IP_LINEAR_SYSTEM_HPP
