#ifndef IP_PRIMAL_DUAL_DIR_HPP
#define IP_PRIMAL_DUAL_DIR_HPP

#include "ip.hpp"
#include "mpc_discretize.hpp"

//Solve Linear System Ax = b 

#define HESSIAN_LAGRANGIAN_SIZE 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_EQUALITIES_NUM_ROWS	NUM_EQUALITY_CONSTRAINTS
#define JACOBIAN_EQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

#define JACOBIAN_INEQUALITIES_NUM_ROWS	NUM_INEQUALITY_CONSTRAINTS
#define JACOBIAN_INEQUALITIES_NUM_COLS 	NUM_OPTMIZATION_VARIABLES

#define DIAG_MATRIX_SIGMA_SIZE	NUM_INEQUALITY_CONSTRAINTS

#define VECTOR_SIZE_b0	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_b1	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_b2	NUM_EQUALITY_CONSTRAINTS
#define VECTOR_SIZE_b3	NUM_INEQUALITY_CONSTRAINTS

#define SQ_MATRIX_A_SIZE (HESSIAN_LAGRANGIAN_SIZE + DIAG_MATRIX_SIGMA_SIZE + JACOBIAN_EQUALITIES_NUM_ROWS + JACOBIAN_INEQUALITIES_NUM_ROWS)
#define VECTOR_x_SIZE	 (VECTOR_SIZE_Px + VECTOR_SIZE_Ps + VECTOR_SIZE_Py + VECTOR_SIZE_Pz)
#define VECTOR_b_SIZE	 (VECTOR_SIZE_b0 + VECTOR_SIZE_b1 + VECTOR_SIZE_b2 + VECTOR_SIZE_b3)


void Compute_primal_dual_direction (const struct_ip_vars &s_ip_vars, struct_primal_dual_direction &s_primal_dual_dir);
void Compute_Hessian_Lagrangian(const struct_ip_vars &s_ip_vars, double Hessian_Lagrangian[HESSIAN_LAGRANGIAN_SIZE][HESSIAN_LAGRANGIAN_SIZE]);
void Compute_Jacobian_Equalities(const struct_ip_vars &s_ip_vars, double Jacobian_Equalities[JACOBIAN_EQUALITIES_NUM_ROWS][JACOBIAN_EQUALITIES_NUM_COLS]);
void Compute_Jacobian_Inequalities(const struct_ip_vars &s_ip_vars, double Jacobian_Inequalities[JACOBIAN_INEQUALITIES_NUM_ROWS][JACOBIAN_INEQUALITIES_NUM_COLS]);
void Compute_Diag_Matrix_Sigma(const struct_ip_vars &s_ip_vars, double Diag_Matrix_Sigma[DIAG_MATRIX_SIGMA_SIZE]);
void Compute_Gradient_Lagrangian(const struct_ip_vars &s_ip_vars, double Gradient_Lagrangian[NUM_OPTMIZATION_VARIABLES]);
void Compute_vector_b1(const struct_ip_vars &s_ip_vars, double vector_b1[VECTOR_SIZE_b1]);
void Compute_vector_b2(const struct_ip_vars &s_ip_vars, double vector_b2[VECTOR_SIZE_b2]);
void Compute_vector_b3(const struct_ip_vars &s_ip_vars, double vector_b3[VECTOR_SIZE_b3]);
double Compute_Lagrangian(const struct_ip_vars &s_ip_vars);

int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

#endif // IP_PRIMAL_DUAL_DIR_HPP
