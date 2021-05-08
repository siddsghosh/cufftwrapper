#include "openacc.h"

void *getaccstrm()
{
	return( acc_get_cuda_stream(acc_async_sync) );
}
