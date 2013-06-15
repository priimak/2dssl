#include <stdio.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_specfunc.h>

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf("%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

/* Below are constants specific to this problem */

// mass of the electron
#define mE 1.0

// electron charge
#define e 1.0

// hbar - Plank constant
#define h 1.0

// Boltzman constant
#define k 1.0

#define PI 3.141592653589793115998
#define TWO_PI (2*PI)

// Temperature of the system. Was 5.0
#ifdef T_EXT

#define T T_EXT

#else

#define T 0.1

#endif

// we cannot use long double on our GPU
#define ffloat double

