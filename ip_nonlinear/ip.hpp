#ifndef IP_HPP
#define IP_HPP

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "utils.hpp"

//Interior Point Algrithm parameters

//#define TAU	0.995
//#define ERROR_TOL_MU		1.0e-10
#define	ERROR_TOL_TOTAL		1.0e-5
#define SIGMA_MU			0.5
#define ESP_DIFFERENTIATION	1.0e-5
//#define MU	1.0
#define NU	1
#define ETA	0.3
#define ALPHA_BACKTRACK_RATIO 0.5

#define NUM_ERROR_TERMS 3

//Optimization problem formulation parameters

#define NUM_OPTMIZATION_VARIABLES	2
#define NUM_INEQUALITY_CONSTRAINTS  6

//primal dual direction parameters

#define VECTOR_SIZE_Px	NUM_OPTMIZATION_VARIABLES
#define VECTOR_SIZE_Pw 	NUM_INEQUALITY_CONSTRAINTS
#define VECTOR_SIZE_Py	NUM_INEQUALITY_CONSTRAINTS

//typepdefs

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double W[NUM_INEQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS]; //y
	
	double mu;
	double nu;
	double eta;
	
} struct_ip_vars;

typedef struct
{
	double Vector_Px [VECTOR_SIZE_Px];
	double Vector_Pw [VECTOR_SIZE_Pw];
	double Vector_Py [VECTOR_SIZE_Py];
		
} struct_primal_dual_direction;


#endif // IP_HPP
