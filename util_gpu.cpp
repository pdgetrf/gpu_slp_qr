#include "util_gpu.h"
#include "stdio.h"

extern "C"
void printout_devices( )
{
	int ndevices, idevice;
	cudaGetDeviceCount( &ndevices );

	for( idevice = 0; idevice < 1; idevice++ ) 
	//for( idevice = 0; idevice < ndevices; idevice++ ) 
	{
		cudaDeviceProp prop;
		cudaGetDeviceProperties( &prop, idevice );
		printf( "device %d: %s, %.1f MHz clock, %.1f MB memory, capability %d.%d\n",
				idevice,
				prop.name,
				prop.clockRate / 1000.,
				prop.totalGlobalMem / (1024.*1024.),
				prop.major,
				prop.minor );
	}
}

