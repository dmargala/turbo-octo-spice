// Created 25-Jan-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include <iostream>


int main(int argc, char **argv) {

    int deviceCount;

    cudaGetDeviceCount(&deviceCount);

    if (deviceCount == 0) {
        std::cout << "No CUDA GPU has been detected" << std::endl;
        return -1;
    } 
    else {
        std::cout << "Number CUDA GPU devices detected: " << deviceCount << std::endl;
    }

    for (int dev = 0; dev < deviceCount; dev++) {
        cudaDeviceProp deviceProp;

        cudaGetDeviceProperties(&deviceProp, dev);

        std::cout << "Device " << dev << " name: " << deviceProp.name << std::endl;
        std::cout << " Computational Capabilities: " << deviceProp.major << "." << deviceProp.minor << std::endl;
        std::cout << " Maximum global memory size: " << deviceProp.totalGlobalMem << std::endl;
        std::cout << " Maximum constant memory size: " << deviceProp.totalConstMem << std::endl;
        std::cout << " Maximum shared memory size per block: " << deviceProp.sharedMemPerBlock << std::endl;
        std::cout << " Maximum block dimensions: " << deviceProp.maxThreadsDim[0] << " x "
            << deviceProp.maxThreadsDim[1] << " x " << deviceProp.maxThreadsDim[2] << std::endl;
        std::cout << " Maximum grid dimensions: " << deviceProp.maxGridSize[0] << " x "
            << deviceProp.maxGridSize[1] << " x " << deviceProp.maxGridSize[2] << std::endl;
        std::cout << " Warp size: " << deviceProp.warpSize << std::endl;
    }

    return 0;
}
