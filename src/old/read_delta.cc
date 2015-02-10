//clang++ -lhdf5 -lhdf5_cpp read_delta.cc

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "H5Cpp.h"
#include "boost/program_options.hpp"

namespace po = boost::program_options;

const double lyA = 1216.;
unsigned long testcounter = 0;

// Operator function
extern "C" herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);

herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
    hid_t group;
    /*
     * Open the group using its name.
     */
    group = H5Gopen2(loc_id, name, H5P_DEFAULT);
    /*
     * Display group name.
     */
    //cout << "Name : " << name << endl;
    if((testcounter % 10000) == 0) std::cout << testcounter << " : " << name << std::endl;

    testcounter++;

    H5Gclose(group);
    return 0;
}


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

class LOSDeltaField {
public:
    LOSDeltaField(const std::string& targetId, double ra, double dec, double z,
    std::vector<float>& loglam, std::vector<float>& delta, std::vector<float>& ivar) 
    : _targetId(targetId), _ra(ra), _dec(dec), _z(z) {
        _loglam.swap(loglam);
        _delta.swap(delta);
        _ivar.swap(ivar);
    }
    double meanDelta() {
        return mean(_delta);
    }
    double meanLogLam() {
        return mean(_loglam);
    }
    double meanIVar() {
        return mean(_ivar);
    }
    double ra() { return _ra; }
    double dec() { return _dec; }
    double z() { return _z; }
private:
    double _ra, _dec, _z;
    std::string _targetId;
    std::vector<float> _loglam, _delta, _ivar;
};

class HDF5DeltaField {
public:
    HDF5DeltaField(const std::string& filename) {
        _file = H5::H5File(filename, H5F_ACC_RDONLY);
    }
    ~HDF5DeltaField() {
        _file.close();
    }

    int numTargets() {
        H5::Group grp = _file.openGroup("delta_field");
        int n(grp.getNumObjs());
        grp.close();
        return n;
    }

    void getTargetNames() {
        H5::Group grp = _file.openGroup("delta_field");
        for(int i(0); i < 5; ++i){
            std::string targetId(grp.getObjnameByIdx(i));
            std::cout << targetId << std::endl;
        }
        grp.close();
    }

    std::vector<LOSDeltaField> loadTargets(bool verbose=false) {
        int n = numTargets();
        H5::Group grp = _file.openGroup("delta_field");
        std::string targetId;
        std::vector<LOSDeltaField> losVector;
        // for(int i(0); i < n; ++i){
        //     targetId = grp.getObjnameByIdx(i);
        //     if(verbose & ((i % 10000) == 0)) std::cout << i << " : " << targetId << std::endl;
        //     //loadTarget(targetId);
        // }

        herr_t idx = H5Literate(grp.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);

        grp.close();
        return losVector;
    }

    LOSDeltaField loadTarget(const std::string& targetId) {
        H5::Group targetGroup = _file.openGroup("delta_field/" + targetId);

        // Read attributes
        double ra, dec, z;
        readAttribute(targetGroup, "ra", ra);
        readAttribute(targetGroup, "dec", dec);
        readAttribute(targetGroup, "z", z);

        // Read datasets
        std::vector<float> redshifts, delta, ivar, loglam;
        //readDataSet(targetGroup, "absorber_z", redshifts);
        //readDataSet(targetGroup, "absorber_delta", delta);
        //readDataSet(targetGroup, "absorber_ivar", ivar);
        // Convert redshifts back to loglam
        //for(const auto& redshift : redshifts) {
            //loglam.push_back(std::log10(lyA*(1+redshift)));
        //}

        // clean up
        targetGroup.close();

        return LOSDeltaField(targetId, ra, dec, z, loglam, delta, ivar);
    }

private:
    H5::H5File _file;
};

int main(int argc, char **argv) {

    std::string target_id, filename;
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

    // open delta field file
    HDF5DeltaField deltaField(filename);

    std::cout << "num targets: " << deltaField.numTargets() << std::endl;

    deltaField.getTargetNames();

    std::vector<LOSDeltaField>&& losVector = deltaField.loadTargets(true);

    // load target's line of sight delta field
    LOSDeltaField&& losDF = deltaField.loadTarget(target_id);

    // print info
    std::cout << "ra: " << losDF.ra() << ", dec: " << losDF.dec() << ", z: " << losDF.z() << std::endl;
    std::cout << "mean loglam: " << losDF.meanLogLam() << std::endl;
    std::cout << "mean delta: " << losDF.meanDelta() << std::endl;
    std::cout << "mean ivar: " << losDF.meanIVar() << std::endl;

    return 0;
}
