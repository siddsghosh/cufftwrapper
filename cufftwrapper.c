#include "wrap_common.h"
#include <stdio.h>

void cud2z_( int *nx, int *ny, int *nz, cufftDoubleReal *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, (*nx+2), ((*nx)/2+1), nb, CUFFT_D2Z );
	gpuErrchk( cufftExecD2Z(*plan, a, (cufftDoubleComplex*)a) );
}

void cuz2d_( int *nx, int *ny, int *nz, cufftDoubleComplex *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, ((*nx)/2+1), (*nx), nb, CUFFT_Z2D );
	gpuErrchk( cufftExecZ2D(*plan, a, (cufftDoubleReal*)a) );
}

void cuz2zf_( int *nx, int *ny, int *nz, cufftDoubleComplex *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, *nx, *nx, nb, CUFFT_Z2Z );
	gpuErrchk( cufftExecZ2Z(*plan, a, a, CUFFT_FORWARD ) );
}

void cuz2zb_( int *nx, int *ny, int *nz, cufftDoubleComplex *a )
{
	cufftHandle *plan; 
	int nb;

	nb = (*ny)*(*nz);
	plan = getPlan( *nx, *nx, *nx, nb, CUFFT_Z2Z );
	gpuErrchk( cufftExecZ2Z(*plan, a, a, CUFFT_FORWARD ) );
}
