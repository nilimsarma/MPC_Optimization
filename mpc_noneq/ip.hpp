#ifndef IP_HPP
#define IP_HPP

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

//mpc formulation parameters

#define T_START	0.0
#define T_END	1.0
#define NUM_DISCRETIZATION (5 + 1)
#define INTERVAL_SIZE ( (T_END - T_START) / (NUM_DISCRETIZATION - 1) )

#define NUM_STATE_VARIABLES 	4
#define NUM_CONTROL_VARIABLES 	2
#define NUM_PARAMETERS			4

#define NUM_PATH_CONSTRAINTS		4
#define NUM_TERMINAL_CONSTRAINTS 	0

//Optimization problem formulation parameters

#define TEST

#ifndef TEST

#define NUM_OPTMIZATION_VARIABLES	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES + NUM_CONTROL_VARIABLES) )
#define NUM_EQUALITY_CONSTRAINTS	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES) )
#define NUM_INEQUALITY_CONSTRAINTS  ( (NUM_DISCRETIZATION * (NUM_PATH_CONSTRAINTS + NUM_TERMINAL_CONSTRAINTS) ) + 2*NUM_EQUALITY_CONSTRAINTS)

#else

#define NUM_OPTMIZATION_VARIABLES	2
#define NUM_EQUALITY_CONSTRAINTS	0
#define NUM_INEQUALITY_CONSTRAINTS  (4 + 2 * NUM_EQUALITY_CONSTRAINTS)

#endif

//Interior Point Algrithm parameters

#define TAU	0.95
#define	ERROR_TOL_TOTAL		1.0e-5
#define SIGMA_MU			0.5
#define ESP_DIFFERENTIATION	1.0e-5
#define NU	1.0
#define ETA	0.3
#define ALPHA_BACKTRACK_RATIO 0.5
#define NUM_ERROR_TERMS 3

//primal dual direction parameters

#define VECTOR_SIZE_Px	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_Ps 	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Pz	NUM_INEQUALITY_CONSTRAINTS

//typepdefs

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double Slack[NUM_INEQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS];  //z

	double mu;
	double nu;
	
} struct_ip_vars;

typedef struct
{
	double Vector_Px [VECTOR_SIZE_Px];
	double Vector_Ps [VECTOR_SIZE_Ps];
	double Vector_Pz [VECTOR_SIZE_Pz];

} struct_primal_dual_direction;

typedef struct  
{
	double alpha_s;
	double alpha_z;
	
} struct_alpha;


#define LOG_OPTIMIZATION_VARIABLES 0
#define LOG_INEQUALITY_CONSTRAINTS 0
#define LOG_EQUALITY_CONSTRAINTS 0
#define LOG_SLACK_VARIABLES 0
#define LOG_MU_VALUE 1
#define LOG_NU_VALUE 1
#define LOG_ALPHA_VALUE 0
#define LOG_ERROR_VALUE 1

#endif // IP_HPP
