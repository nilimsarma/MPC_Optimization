#ifndef IP_HPP
#define IP_HPP

#include <stdio.h>
#include <stdlib.h>
#include "math.h"


//Interior Point Algrithm parameters

#define	ERROR_TOL_TOTAL		1.0e-5
#define SIGMA_MU			0.5
#define ESP_DIFFERENTIATION	1.0e-5

#define NUM_ERROR_TERMS 3

//Optimization problem formulation parameters

#define NUM_OPTMIZATION_VARIABLES	2
#define NUM_INEQUALITY_CONSTRAINTS  4
#define NUM_EQUALITY_CONSTRAINTS	1

//typepdefs

typedef struct
{
	double Optimization_variables[NUM_OPTMIZATION_VARIABLES];
	double Lagrange_multiplier_equality[NUM_EQUALITY_CONSTRAINTS];
	double Lagrange_multiplier_inequality[NUM_INEQUALITY_CONSTRAINTS];
	double Slack[NUM_INEQUALITY_CONSTRAINTS];
	
} struct_ip_vars;


#endif	// IP_HPP
