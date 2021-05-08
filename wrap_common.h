#include <stdlib.h>
#include <stdio.h>
#include <cufft.h>

#define gpuErrchk(ans) { checkCuda((ans), __FILE__, __LINE__); }
static void checkCuda(const cudaError_t code, const char *file, int line )
{  
   if (code != cudaSuccess)
   {  
      fprintf(stderr,"CUDA Error name: %s string: %s Src: %s Line: %d\n", 
         cudaGetErrorName(code),cudaGetErrorString(code), file, line);
	fflush(stderr);
	exit(code);
   }
}

cufftHandle *getPlan( int , int , int , int , cufftType );
