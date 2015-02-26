// #include <thrust/version.h>

#ifndef BRUTEGPU
#define BRUTEGPU

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>

#include <cuda_runtime.h>

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

void bruteGPU(std::vector<std::vector<double> > &columns, double min, double max, int nbins, std::vector<double> &xi, long chunksize);

#endif