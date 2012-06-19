#include "cuda.h"
#include "cuda_runtime_api.h"
#include "cublas.h"
#include "math.h"
#include "ctype.h"

#define TESTING_CUDA_INIT()                     \
	  CUdevice  dev;                            \
  CUcontext context;                            \
  if( CUDA_SUCCESS != cuInit( 0 ) ) {                   \
	      fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);      \
	    }                                 \
  if( CUDA_SUCCESS != cuDeviceGet( &dev, 0 ) ) {            \
	      fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);     \
	    }                                 \
  if( CUDA_SUCCESS != cuCtxCreate( &context, 0, dev ) ) {       \
	      fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1); \
	    }                                 \
  if( CUBLAS_STATUS_SUCCESS != cublasInit( ) ) {                        \
	      fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);     \
	    }                                 \
  printout_devices( );

#define TESTING_CUDA_FINALIZE()         \
	  cuCtxDetach( context );           \
  cublasShutdown();

