#ifndef UTILS_CPP
#define UTILS_CPP

#include "utils.hpp"

void print_matrix(double* A, int num_rows, int num_cols)
{
	int i,j;

	printf("\n");
	for(i = 0; i < num_rows; i++)
	{
		for(j = 0; j < num_cols; j++)
		{
			printf("%f, ", A[i*num_cols + j]);
		}
		printf(";\n");
	}
}

void print_vector(double* A, int vec_size)
{
	int i,j;

	printf("\n");
	for(i = 0; i < vec_size; i++)
	{
		printf("%f, ", A[i]);
	
	}
	printf("\n");
}

#endif //UTILS_CPP
