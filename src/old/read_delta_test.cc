//clang++ -lhdf5 -lhdf5_cpp read_delta.cc

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "H5Cpp.h"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

const double lyA = 1216.;

template<typename T> void readDataSet(H5::Group& grp, const std::string& name, std::vector<T>& data){
    // Open the data set
    H5::DataSet dataset = grp.openDataSet(name);
    // Need to figure out dataset dimensions
    H5::DataSpace dataspace = dataset.getSpace();
    // Get the number of dimensions in the dataspace.
    // int rank = dataspace.getSimpleExtentNdims();
    // std::cout << "\tndims: " << rank << std::endl;
    hsize_t dims_out[1];
    dataspace.getSimpleExtentDims(dims_out, NULL);
    int size = dims_out[0];
    // create a vector the same size as the dataset
    data.resize(size);
    // read dataset into vector
    dataset.read(data.data(), H5::PredType::NATIVE_FLOAT);
    // close dataset
    dataset.close();
}

template<typename T> void readAttribute(H5::Group& grp, const std::string& name, T& value){
    H5::Attribute attr;
    attr = grp.openAttribute(name);
    attr.read(attr.getDataType(), &value);
}

template<typename T> double mean(const std::vector<T>& data){
    unsigned long size(data.size());
    double sum(0);
    for(const auto& val : data) {
        sum += val;
    }
    return sum/size;
}

int main(int argc, char **argv) {

    std::string target_id, filename;

    //std::string filename("/Users/daniel/source/qusp/absorption/lc-all-absorber-delta-field.hdf5");
    std::string delta_field_grp("delta_field");

    po::options_description cli("Read HDF5 file");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose,v", "Prints more information.")
        ("target", po::value<std::string>(&target_id)->default_value("3615-55445-8"),
            "target id string")
        ("input", po::value<std::string>(&filename)->default_value(""),
            "input filename")
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
    bool verbose(vm.count("verbose"));

    // Try block to detect exceptions raised by any of the calls inside it
    try {
    
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        H5::Exception::dontPrint();

        // open file and list top level objects
        H5::H5File file(filename, H5F_ACC_RDONLY);
        if (verbose) std::cout << "file '" << filename << "' is open" << std::endl;
        if (verbose) std::cout << "\tnumber of objs: " << file.getNumObjs() << std::endl;

        // open delta field group and show number of objects
        H5::Group group = file.openGroup(delta_field_grp);
        if (verbose) std::cout << "group '" << delta_field_grp << "' is open" << std::endl;
        if (verbose) std::cout << "\tnumber of objs: " << group.getNumObjs() << std::endl;

        // open group corresponding to specified target 
        H5::Group target_group;
        try{
            target_group = group.openGroup(target_id);
        }
        catch(H5::GroupIException not_found_error) {
            std::cout << "group for specified target not found." << std::endl;
            return -1;
        }
        if (verbose) std::cout << "group '" << target_id << "' is open" << std::endl;
        if (verbose) std::cout << "\tnumber of objs: " << target_group.getNumObjs() << std::endl;
        if (verbose) std::cout << "\tnumber of attrs: " << target_group.getNumAttrs() << std::endl;
        if (verbose) std::cout << std::endl;

        // Read attributes
        double ra, dec, z;
        readAttribute(target_group, "ra", ra);
        readAttribute(target_group, "dec", dec);
        readAttribute(target_group, "z", z);
        std::cout << "ra: " << ra << ", dec: " << dec << ", z: " << z << std::endl;

        // Read datasets
        std::vector<float> redshifts, delta, ivar, loglam;
        readDataSet(target_group, "absorber_z", redshifts);
        readDataSet(target_group, "absorber_delta", delta);
        readDataSet(target_group, "absorber_ivar", ivar);
        std::cout << "mean redshift: " << mean(redshifts) << std::endl;
        std::cout << "mean delta: " << mean(delta) << std::endl;
        std::cout << "mean ivar: " << mean(ivar) << std::endl;

        // convert redshift to loglam
        for(const auto& redshift : redshifts) {
            loglam.push_back(std::log10(lyA*(1+redshift)));
        }
        std::cout << "mean loglam: " << mean(loglam) << std::endl;

        // clean up
        target_group.close();
        group.close();
        file.close();

    } // end of try block
    // catch failure caused by the H5File operations
    catch(H5::FileIException error) {
        error.printError();
        return -1;
    }
    catch(H5::GroupIException error) {
        error.printError();
        return -1;
    }
    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error) {
        error.printError();
        return -1;
    }
    // catch failure caused by the DataSpace operations
    catch(H5::DataSpaceIException error) {
        error.printError();
        return -1;
    }
    // catch failure caused by the Attribute operations
    catch(H5::AttributeIException error){
        error.printError();
        return -1;
    }
    return 0;
}
