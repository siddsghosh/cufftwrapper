#include <stdlib.h>
#include "wrap_common.h"

#define MAXPLANS 8 // Number of plans needed to be stored for LES

struct plantblt{
   int n, idist, odist, nb;
   cufftType type;
   cufftHandle plan;
};
int init_state = 0;
struct plantblt* plantbl;

void *getaccstrm(void);

void initialize()
{
	int i;

	plantbl = (struct plantblt*)malloc(MAXPLANS*sizeof(struct plantblt));
	for(i=0; i<MAXPLANS; i++)plantbl[i].n = 0;
	init_state = 1;
}

int isEqualPlan( int i, int n, int idist, int odist, int nb, cufftType tfft )
{
	if ( plantbl[i].n     != 0     && plantbl[i].n     == n     && 
             plantbl[i].idist == idist && plantbl[i].odist == odist && 
             plantbl[i].nb    == nb    && plantbl[i].type  == tfft ){
		return 1;
	}else{
		return 0;
	}
}

cufftHandle * setPlan( int i, int n, int idist, int odist, int nb, cufftType tfft )
{
	cufftHandle tplan;
	int ns[1], iem[] = {1}, oem[] = {1}, idx;

	ns[0] = n;
	gpuErrchk( cufftPlanMany(&plantbl[i].plan,  // plan 
                                                1,  // rank
                                               ns,  // size
                                              iem,  // pointer for later input layout
                                                1,  // Input stride of consecutive data in a sample
                                            idist,  // Input distance between 2 consecutive samples
                                              oem,  // pointer for later output layout
                                                1,  // Output stride of consecutive data in the transform
                                            odist,  // Output distance between 2 consecutive transforms
                                             tfft,  // FFT type 
                                               nb));// Number of 1-D FFTs
	plantbl[i].n     = n;
	plantbl[i].idist = idist;
	plantbl[i].odist = odist;
	plantbl[i].nb    = nb;
	plantbl[i].type  = tfft;
	cufftSetStream(plantbl[i].plan, (cudaStream_t)getaccstrm());
	return &plantbl[i].plan;
}

cufftHandle * getPlan( int n, int idist, int odist, int nb, cufftType tfft )
{
	int i;

	if (init_state == 0)initialize();
	for(i=0; i<MAXPLANS; i++){
		if ( isEqualPlan( i, n, idist, odist, nb, tfft ) == 1){
			return &plantbl[i].plan;
		}else{
			return setPlan( i, n, idist, odist, nb, tfft );
		}
	}
	printf("Please increase MAXPLANS and rerun...!");
	exit(-1);
}

