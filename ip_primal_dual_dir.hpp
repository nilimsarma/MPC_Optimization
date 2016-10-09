#ifndef IP_PRIMAL_DUAL_DIR_HPP
#define IP_PRIMAL_DUAL_DIR_HPP

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

#define ESP_DIFFERENTIATION

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double mu;
	double nu;
	double eta;
	double S[NUM_INEQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_equality[NUM_EQUALITY_CONSTRAINTS];	//y
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS]; //z
} struct_ip_vars;

typedef struct
{
	double Vector_Px [VECTOR_SIZE_Px];
	double Vector_Ps [VECTOR_SIZE_Ps];
	double Vector_Py [VECTOR_SIZE_Py];
	double Vector_Pz [VECTOR_SIZE_Pz];
		
} struct_primal_dual_direction;

void Linear_System_Setup(struct_ip_vars s_ip_vars);
double [HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE] Compute_Hessian_Lagrangian(struct_ip_vars s_ip_vars);
double [JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS] Compute_Jacobian_Equalities(struct_ip_vars s_ip_vars);
double [JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS] Compute_Jacobian_Inequalities(struct_ip_vars s_ip_vars);
double [DIAG_MATRIX_SIGMA_SIZE] Compute_Diag_Matrix_Sigma(struct_ip_vars s_ip_vars);
double [NUM_OPTMIZATION_VARIABLES] Compute_Gradient_Lagrangian(struct_ip_vars s_ip_vars);
double [VECTOR_SIZE_b1] Compute_vector_b1(struct_ip_vars s_ip_vars);
double [VECTOR_SIZE_b2] Compute_vector_b2(struct_ip_vars s_ip_vars);
double [VECTOR_SIZE_b3] Compute_vector_b3(struct_ip_vars s_ip_vars);
double Compute_Lagrangian(struct_ip_vars s_ip_vars);


#endif // IP_PRIMAL_DUAL_DIR_HPP