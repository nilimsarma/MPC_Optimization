#ifndef MPC_FORMULATE_HPP
#define MPC_FORMULATE_HPP

#include <stdio.h>
#include <stdlib.h>

#define T_START	0
#define T_END	6
#define NUM_DISCRETIZATION (40 + 1)
#define INTERVAL_SIZE ( (T_END - T_START) / (NUM_DISCRETIZATION - 1) )

#define NUM_STATE_VARIABLES 	11
#define NUM_CONTROL_VARIABLES 	6
#define NUM_PARAMETERS			87

#define NUM_PATH_CONSTRAINTS	
#define NUM_TERMINAL_CONSTRAINTS 	1

/********************************/

#define NUM_OPTMIZATION_VARIABLES	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES + NUM_CONTROL_VARIABLES) )
#define NUM_EQUALITY_CONSTRAINTS	( NUM_DISCRETIZATION * (NUM_STATE_VARIABLES) )
#define NUM_INEQUALITY_CONSTRAINTS  ( NUM_DISCRETIZATION * (NUM_PATH_CONSTRAINTS + NUM_TERMINAL_CONSTRAINTS) )

double* Evaluate_state_differential (double* State_variables, double* Control_variables);
double* Evaluate_constraints (double* State_variables, double* Control_variables);
double Evaluate_objective (double* State_variables, double* Control_variables);

#endif // MPC_FORMULATE_HPP