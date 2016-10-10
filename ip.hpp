#ifndef IP_HPP
#define IP_HPP

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

//Interior Point Algrithm parameters

#define TAU	0.995
#define ERROR_TOL_MU		1.0e-10
#define	ERROR_TOL_TOTAL		1.0e-10
#define SIGMA_MU			0.9
#define ESP_DIFFERENTIATION	1.0e-20	
#define MU	1.0
#define NU	100
#define ETA	0.9
#define ALPHA_BACKTRACK_VAL 0.1

#define NUM_ERROR_TERMS 4

//mpc formulation parameters

#define T_START	0
#define T_END	1
#define NUM_DISCRETIZATION (20 + 1)
#define INTERVAL_SIZE ( (T_END - T_START) / (NUM_DISCRETIZATION - 1) )

#define NUM_STATE_VARIABLES 	5
#define NUM_CONTROL_VARIABLES 	2
#define NUM_PARAMETERS			4

#define NUM_PATH_CONSTRAINTS		4
#define NUM_TERMINAL_CONSTRAINTS 	0

//Optimization problem formulation parameters

#define NUM_OPTMIZATION_VARIABLES	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES + NUM_CONTROL_VARIABLES) )
#define NUM_EQUALITY_CONSTRAINTS	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES) )
#define NUM_INEQUALITY_CONSTRAINTS  ( NUM_DISCRETIZATION * (NUM_PATH_CONSTRAINTS + NUM_TERMINAL_CONSTRAINTS) )

//primal dual direction parameters

#define VECTOR_SIZE_Px	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_Ps 	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Py	NUM_EQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Pz	NUM_INEQUALITY_CONSTRAINTS

//typepdefs

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double S[NUM_INEQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_equality[NUM_EQUALITY_CONSTRAINTS];	//y
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS]; //z
	
	double mu;
	double nu;
	double eta;
	
} struct_ip_vars;

typedef struct
{
	double Vector_Px [VECTOR_SIZE_Px];
	double Vector_Ps [VECTOR_SIZE_Ps];
	double Vector_Py [VECTOR_SIZE_Py];
	double Vector_Pz [VECTOR_SIZE_Pz];
		
} struct_primal_dual_direction;

typedef struct  
{
	double alpha_s;
	double alpha_z;
	
} struct_alpha;


#endif // IP_HPP
