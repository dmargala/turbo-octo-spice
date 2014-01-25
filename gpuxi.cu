// Created 08-Aug-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Correlation function estimator using GPU.

// Compile using:
// nvcc -m64 -arch=sm_20 -lboost_program_options -llikely -lcosmo gpuxi.cu -o gpuxi

// Example usage:
// time -p ./gpuxi -i /Users/daniel/Cosmo/LyAlpha/cosmo/build/delta.dat --verbose

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include <thrust/version.h>

#include "/Users/daniel/source/gpu/cuda_by_example/common/book.h"

namespace po = boost::program_options;
namespace lk = likely;

struct DataStruct {
    float x, y, z, d, w;
};

__global__ void histo_kernel(DataStruct *b1, DataStruct *b2, long size, float *dsum, float *wsum, 
float min, float max, int nbins, int maxbins, bool diag) {
    // Temp histogram is dynamically allocated
    extern __shared__ float shared[];
    float *tempd = (float*) &shared[0];
    float *tempw = (float*) &shared[maxbins];
    tempd[threadIdx.x] = 0; 
    tempw[threadIdx.x] = 0;
    __syncthreads();

    unsigned long i = threadIdx.x + blockIdx.x * blockDim.x; 
    unsigned long offset = blockDim.x * gridDim.x;

    float spacing = (max-min)/nbins;

    while (i < size) {
        float xi = b1[i].x;
        float yi = b1[i].y;
        float zi = b1[i].z;
        float di = b1[i].d;
        float wi = b1[i].w;
        for(unsigned long j = 0; j < size; ++j) {
            float dx = xi - b2[j].x;
            float dy = yi - b2[j].y;
            float dz = zi - b2[j].z;

            float separation = std::sqrt(dx*dx+dy*dy+dz*dz);
            float wgt = wi*b2[j].w;
            int index;
            if(separation < min){
                index = 0;
            }
            else if(separation >= max) {
                index = nbins+1;
            }
            else {
                index = floor((separation-min)/spacing)+1;
            }
            if(diag && j <= i) wgt = 0;
            atomicAdd(&tempd[index], wgt*di*b2[j].d);
            atomicAdd(&tempw[index], wgt);
        }
        i += offset;
    }

    __syncthreads();

    atomicAdd(&(dsum[threadIdx.x]), tempd[threadIdx.x]);
    atomicAdd(&(wsum[threadIdx.x]), tempw[threadIdx.x]);
}

void bruteGPU(std::vector<std::vector<double> > const &columns, double min, double max, 
int nbins, std::vector<double> &xi, long chunksize) {

    long nrows = columns[0].size();

    assert(nrows % chunksize == 0);

    int nchunks = nrows / chunksize;// / chunksize;

    std::cout << "nchunks: " << nchunks << std::endl;
    std::cout << "chunksize: " << chunksize << std::endl;

    DataStruct *data = (DataStruct*) malloc(nrows * sizeof(DataStruct));

    for(int i = 0; i < nrows; ++i) {
        data[i].x = columns[0][i];
        data[i].y = columns[1][i];
        data[i].z = columns[2][i];
        data[i].d = columns[3][i];
        data[i].w = columns[4][i];
    }

    std::cout << "sizeof data: " << nrows*sizeof(DataStruct)/1024./1024. << " MB" << std::endl;

    // Look up device properties
    cudaDeviceProp prop;
    HANDLE_ERROR( cudaGetDeviceProperties( &prop, 0 ) );
    int blocks = 2*prop.multiProcessorCount; 
    int threadsPerBlock = 256;

    std::cout << "num blocks: " << blocks << std::endl;
    std::cout << "threadsPerBlock: " << threadsPerBlock << std::endl;

    int nhistbins = threadsPerBlock;

    assert(nhistbins >= nbins+2);

    float dsum[nhistbins];
    float wsum[nhistbins];

    std::vector<double> tempxi(nbins,0);
    std::vector<double> counts(nbins,0);

    // allocate memory on the GPU for the file's data
    DataStruct *dev_b1, *dev_b2;
    float *dev_dsum;
    float *dev_wsum;

    HANDLE_ERROR( cudaMalloc( (void**)&dev_b1, chunksize * sizeof(DataStruct) ) ); 
    HANDLE_ERROR( cudaMalloc( (void**)&dev_b2, chunksize * sizeof(DataStruct) ) ); 
    HANDLE_ERROR( cudaMalloc( (void**)&dev_dsum, nhistbins * sizeof( float ) ) );
    HANDLE_ERROR( cudaMalloc( (void**)&dev_wsum, nhistbins * sizeof( float ) ) );

    long totalcounts = 0;

    cudaEvent_t start, stop;
    HANDLE_ERROR( cudaEventCreate( &start ) );
    HANDLE_ERROR( cudaEventCreate( &stop ) );

    long sharedMemoryPerBlock = 2*threadsPerBlock*sizeof(float);

    std::cout << "Shared memory per block: " << sharedMemoryPerBlock/1024. << " KB" << std::endl;

    double totalElapsedTime = 0;
    for(int ichunk = 0; ichunk < nchunks; ++ichunk) {
        
        HANDLE_ERROR( cudaEventRecord( start, 0 ) );

        for(int jchunk = 0; jchunk <= ichunk; ++jchunk) {

            //std::cout << "Starting chunk (" << ichunk << "," << jchunk << ")..." << std::endl;

            HANDLE_ERROR( cudaMemcpy( dev_b1, &data[ichunk*chunksize], 
                chunksize * sizeof(DataStruct), cudaMemcpyHostToDevice ) );
            HANDLE_ERROR( cudaMemcpy( dev_b2, &data[jchunk*chunksize], 
                chunksize * sizeof(DataStruct), cudaMemcpyHostToDevice ) );
            HANDLE_ERROR( cudaMemset( dev_dsum, 0, nhistbins * sizeof( float ) ) );
            HANDLE_ERROR( cudaMemset( dev_wsum, 0, nhistbins * sizeof( float ) ) );
        
            histo_kernel<<<blocks, threadsPerBlock, sharedMemoryPerBlock>>>(dev_b1, dev_b2, 
                chunksize, dev_dsum, dev_wsum, min, max, nbins, nhistbins, ichunk == jchunk);

            HANDLE_ERROR( cudaMemcpy( dsum, dev_dsum, nhistbins * sizeof( float ), cudaMemcpyDeviceToHost ) );
            HANDLE_ERROR( cudaMemcpy( wsum, dev_wsum, nhistbins * sizeof( float ), cudaMemcpyDeviceToHost ) );

            long chunkcounts = 0;
            // Save results from chunk
            //std::cout << wsum[0] << " " << wsum[nbins+1] << std::endl;
            for(int i = 0; i < nhistbins; ++i) {
                chunkcounts += wsum[i];
                if (i <= nbins && i > 0) {
                    //std::cout << i-1 << " " << dsum[i] << std::endl;
                    tempxi[i-1] += dsum[i];
                    counts[i-1] += wsum[i];
                }
            }
            totalcounts += chunkcounts;

            //std::cout << "Chunk (" << ichunk << "," << jchunk << ") counts: " << chunkcounts << std::endl;

            cudaDeviceSynchronize();

        }
        // get stop time, and display the timing results
        HANDLE_ERROR( cudaEventRecord( stop, 0 ) ); 
        HANDLE_ERROR( cudaEventSynchronize( stop ) ); 
        float elapsedTime;
        HANDLE_ERROR( cudaEventElapsedTime( &elapsedTime, start, stop ) );
        totalElapsedTime += elapsedTime;
        printf( "Time to generate (%d):  %3.1f ms\n", ichunk, elapsedTime );
    }

    std::cout << "Total elapsed time: " << totalElapsedTime << " ms" << std::endl;

    long usedcounts = 0;
    for(int i = 0; i < nbins; ++i) {
        usedcounts += counts[i];
        if(counts[i] > 0) tempxi[i] /= counts[i];
    }

    tempxi.swap(xi);

    std::cout << "used " << usedcounts << " of " << totalcounts << " pairs." << std::endl;

    // Free host and device memory
    HANDLE_ERROR( cudaEventDestroy( start ) );
    HANDLE_ERROR( cudaEventDestroy( stop ) );
    cudaFree( dev_dsum ); 
    cudaFree( dev_wsum );
    cudaFree( dev_b1 ); 
    cudaFree( dev_b2 ); 
    free(data);

}

int main(int argc, char **argv) {

    // Configure command-line option processing
    std::string infile,outfile,axis1,axis2;
    long chunksize;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value("xi.dat"),
            "Filename to write correlation function to")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
            "Axis-1 binning")
        ("axis2", po::value<std::string>(&axis2)->default_value("[0:200]*50"),
            "Axis-2 binning")
        ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")
        ("chunksize", po::value<long>(&chunksize)->default_value(4096),
            "Number of chunks to split the dataset into.")
        ;

    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose")),rmu(vm.count("rmu"));

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > columns(5);
    try {
        std::ifstream in(infile.c_str());
        lk::readVectors(in,columns);
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << infile << ": " << e.what() << std::endl;
        return -3;
    }
    if(verbose) {
        std::cout << "Read " << columns[0].size() << " rows from " << infile
            << std::endl;
    }

    // Generate the correlation function grid and run the estimator
    std::vector<double> xi;
    try {
        lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2);
        double x1min(bins1->getBinLowEdge(0)), x1max(bins1->getBinHighEdge(bins1->getNBins()-1));
        double x2min(bins2->getBinLowEdge(0)), x2max(bins2->getBinHighEdge(bins2->getNBins()-1));
        lk::BinnedGrid grid(bins1,bins2);
        int x1nbins = bins1->getNBins();

        bruteGPU(columns,x1min,x1max,x1nbins,xi,chunksize);
    }
    catch(std::exception const &e) {
        std::cerr << "Error while running the estimator: " << e.what() << std::endl;
    }

    // Save the estimator results
    try {
        std::ofstream out(outfile.c_str());
        for(int index = 0; index < xi.size(); ++index) {
            out << index << ' ' << xi[index] << std::endl;
        }
        out.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while saving results: " << e.what() << std::endl;
    }

    return 0;
}
