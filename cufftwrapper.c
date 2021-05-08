#include "wrap_common.h"
#include <stdio.h>

void cud2z_( int *nx, int *ny, int *nz, cufftDoubleReal *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, (*nx+2), ((*nx)/2+1), nb, CUFFT_D2Z );
	if (cufftExecD2Z(*plan, a, (cufftDoubleComplex*)a) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecD2Z Forward failed");
		return;
	}
}

void cuz2d_( int *nx, int *ny, int *nz, cufftDoubleComplex *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, ((*nx)/2+1), (*nx), nb, CUFFT_Z2D );
	if (cufftExecZ2D(*plan, a, (cufftDoubleReal*)a) != CUFFT_SUCCESS){
		fprintf(stderr, "CUFFT error: ExecZ2D Forward failed");
		return;
	}
}
