// Created 16-May-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "XiEstimator.h"
#include "HDF5Delta.h"
#include "HealpixBins.h"
#include "AbsTwoPointGrid.h"
#include "RuntimeError.h"
#include "ThreadPool.h"

#include "boost/bind.hpp"

#include <fstream>

namespace local = turbooctospice;

local::XiEstimator::XiEstimator(int order, std::string infile, cosmo::AbsHomogeneousUniversePtr cosmology,
    local::AbsTwoPointGridPtr grid, local::XiEstimator::BinningCoordinateType type, bool skip_ngc, bool skip_sgc):
    healbins_(local::HealpixBinsI(order)), grid_(grid), coordinate_type_(type), show_progress_(0) {

    // load forest sight lines
    local::HDF5Delta file(infile);
    sightlines_ = file.loadForests(!skip_ngc, !skip_sgc);

    // add sight lines to healpix bins
    unsigned long totalpixels(0);
    double min_loglam(sightlines_[0].pixels[0].loglam), max_loglam(sightlines_[0].pixels[0].loglam);
    for(auto &sightline : sightlines_) {
        int npixels = sightline.pixels.size();
        totalpixels += npixels;
        healbins_.addItem(sightline.ra, sightline.dec, sightline.forest_id);

        // find minimum loglam
        auto loglam_first(sightline.pixels[0].loglam);
        if(loglam_first < min_loglam) {
            min_loglam = loglam_first;
        }
        // find maximum loglam
        auto loglam_last(sightline.pixels[npixels-1].loglam);
        if(loglam_last > max_loglam) {
            max_loglam = loglam_last;
        }
    }

    double zmin(std::pow(10, min_loglam-local::logLyA)-1);
    double zmax(std::pow(10, max_loglam-local::logLyA)-1);

    int numHealBinsOccupied(healbins_.getNBins());
    std::cout << "Read " << totalpixels << " from " << sightlines_.size() << " lines of sight (LOS)" << std::endl;
    std::cout << "Average number of pixels per LOS: " << static_cast<double>(totalpixels)/sightlines_.size() << std::endl;
    std::cout << "Number of Healpix bins occupied: " << numHealBinsOccupied
        << " (" << static_cast<double>(numHealBinsOccupied)/(12*std::pow(4, order)) << ")" << std::endl;

    // the minimum redshift sets the angular scale we will need to consider
    double scale(cosmology->getTransverseComovingScale(zmin));
    max_ang_ = grid_->maxAngularScale(scale);
    std::cout << "Transverse comoving scale at z = " << zmin <<  " (Mpc/h): " << scale << std::endl;
    std::cout << "Max angular scale at z = " << zmin <<  " (rad): " << max_ang_  << std::endl;
    std::cout << "Transverse comoving scale at z = " << zmax <<  " (Mpc/h): "
        << cosmology->getTransverseComovingScale(zmax) << std::endl;
    std::cout << "Max angular scale at z = " << zmax <<  " (rad): "
        << grid_->maxAngularScale(cosmology->getTransverseComovingScale(zmax))  << std::endl;

    // axis binning limits
    r_min = grid_->getAxisMin(0); r_max = grid_->getAxisMax(0); r_spacing = grid_->getAxisBinWidth(0);
    rsq_min = r_min*r_min; rsq_max = r_max*r_max;
    num_mu_bins = grid_->getAxisNBins(1);

    rperp_min = grid_->getAxisMin(1); rperp_max = grid->getAxisMax(1); rperp_spacing = grid_->getAxisBinWidth(1);
    num_rperp_bins = grid_->getAxisNBins(1);
    sep_min = grid_->getAxisMin(1); sep_spacing = grid_->getAxisNBins(1);
    num_sep_bins = grid_->getAxisNBins(1);

    z_min = grid_->getAxisMin(2); z_max = grid_->getAxisMax(2); z_spacing = grid_->getAxisBinWidth(2);
    num_z_bins = grid_->getAxisNBins(2);

    // print stats
    // std::cout << "Number of Healpix bins searched: " << numHealpixBinsSearched << std::endl;
    // line of sight pair statistics
    // unsigned long numLOSPairsTotal = (numLOS*(numLOS-1))/2;
    // double fracLOSPairsConsidered = (double) numLOSPairs / numLOSPairsTotal;
    // double fracLOSPairsUsed = (double) numLOSPairsUsed / numLOSPairs;
    //
    // std::cout << "Number of distinct los pairs " << numLOSPairsTotal << std::endl;
    // std::cout << "considered " << numLOSPairs << " of distinct los pairs. (" << fracLOSPairsConsidered << ")" << std::endl;
    // std::cout << "used " << numLOSPairsUsed << " of los pairs considered. (" << fracLOSPairsUsed << ")" << std::endl;
    //
    // // pixel pair statistics
    // unsigned long numPixelPairsTotal = (numPixels*(numPixels-1))/2;
    // double fracPixelPairsConsidered = (double) numPixelPairs / numPixelPairsTotal;
    // double fracPixelPairsUsed = (double) numPixelPairsUsed / numPixelPairs;
    //
    // std::cout << "Number of distinct pixel pairs " << numPixelPairsTotal << std::endl;
    // std::cout << "considered " << numPixelPairs << " of distinct pixel pairs. (" << fracPixelPairsConsidered << ")" << std::endl;
    // std::cout << "used " << numPixelPairsUsed << " of pixel pairs considered. (" << fracPixelPairsUsed << ")" << std::endl;
};

void local::XiEstimator::run(int nthreads) {
    local::ThreadPool pool(nthreads);
    std::list<local::ThreadPool::Task> tasks;

    auto occupied_bins = healbins_.getOccupiedBins();
    for(auto id : occupied_bins) {
        tasks.push_back(boost::bind(&XiEstimator::healxi_task, this, id));
    }
    show_progress_.restart(sightlines_.size());
    pool.run(tasks);
    // task_finalize();
    tasks.clear();

    // Finalize xi
    int num_xi_bins = grid_->getNBinsTotal();
    xi_.resize(num_xi_bins);
    for(int i = 0; i < num_xi_bins; ++i) {
        tasks.push_back(boost::bind(&XiEstimator::xi_finalize_task, this, i));
    }
    pool.run(tasks);
    tasks.clear();

    // Estimate covariance
    cov_.resize(num_xi_bins);
    for(int a = 0; a < num_xi_bins; ++a) {
        cov_[a].resize(num_xi_bins);
        for(int b = 0; b < num_xi_bins; ++b) {
            tasks.push_back(boost::bind(&XiEstimator::cov_task, this, a, b));
        }
    }

}

bool local::XiEstimator::xi_finalize_task(int id) {
    for(const auto& healxi_entry : healxis_) {
        xi_[id].didj += healxis_[healxi_entry.first][id].didj;
        xi_[id].wgt += healxis_[healxi_entry.first][id].wgt;
    }
    if(xi_[id].wgt > 0) {
        xi_[id].didj /= xi_[id].wgt;
    }
    return true;
}

bool local::XiEstimator::cov_task(int a, int b) {
    for(const auto& healxi_entry : healxis_) {
        cov_[a][b] += healxis_[healxi_entry.first][a].didj*healxis_[healxi_entry.first][b].didj
            - healxis_[healxi_entry.first][a].wgt*healxis_[healxi_entry.first][b].wgt*xi_[a].didj*xi_[b].didj;
    }
    double wgt = xi_[a].wgt*xi_[b].wgt;
    if(wgt > 0) {
        cov_[a][b] /= wgt;
    }
    return true;
}

void local::XiEstimator::task_finalize(){
    int num_xi_bins = grid_->getNBinsTotal();
    // copy data to output vectors
    std::vector<local::XiBin> xisum(num_xi_bins, {});
    for(int index = 0; index < num_xi_bins; ++index) {
        for(const auto& pair : healxis_) {
            xisum[index].didj += healxis_[pair.first][index].didj;
            xisum[index].wgt += healxis_[pair.first][index].wgt;
        }
        if(xisum[index].wgt > 0) {
            xisum[index].didj /= xisum[index].wgt;
        }
    }
    xi_ = xisum;

    // Estimate covariance matrix
    std::vector<std::vector<double> > cov(num_xi_bins, std::vector<double>(num_xi_bins, 0));
    for(int a = 0; a < num_xi_bins; ++a){
        for(int b = 0; b < num_xi_bins; ++b) {
            for(const auto& pair : healxis_) {
                cov[a][b] += healxis_[pair.first][a].didj*healxis_[pair.first][b].didj
                    - healxis_[pair.first][a].wgt*healxis_[pair.first][b].wgt*xisum[a].didj*xisum[b].didj;
            }
            double wgt = xisum[a].wgt*xisum[b].wgt;
            if(wgt > 0) {
                cov[a][b] /= wgt;
            }
        }
    }
    cov_ = cov;
}

void local::XiEstimator::increment_progress() {
    std::lock_guard<std::mutex> lock(show_progress_mutex_);
    ++show_progress_;
    // show_progress_mutex_ is automatically released when lock
    // goes out of scope
}

bool local::XiEstimator::healxi_task(int id) {
    increment_progress();
    // create internal accumulation vectors
    int num_xi_bins = grid_->getNBinsTotal();
    std::vector<double> dsum(num_xi_bins,0), wsum(num_xi_bins,0);
    // to avoid calling trig functions inside loop
    double cos_max_ang(std::cos(max_ang_));
    // allocate temporary vectors before loop
    std::vector<int> neighbors;
    // task_id's xi container
    std::vector<local::XiBin> xi(num_xi_bins, {});
        // Iterate over all sight lines in this healpixel
    for(int primary_los_index : healbins_.getBinContents(id)) {
        auto primary_los = sightlines_[primary_los_index]; // should probably avoid copy
        // Find healpixels within max angular separation of interest
        auto neighbors = healbins_.getBinIndicesWithinRadius(primary_los.ra, primary_los.dec, max_ang_);
        // Iterate over all neighboring healpixels
        for(int neighbor : neighbors) {
            // Check that if there are sightlines in the neighboring healpixel
            if(!healbins_.checkBinExists(neighbor)) continue;
            // Iterate over candidate sightline pairs
            for(int pair_los_index : healbins_.getBinContents(neighbor)) {
                // only count pairs once
                if(pair_los_index <= primary_los_index) continue;
                auto other_los = sightlines_[pair_los_index];
                // check angular separation
                double cos_separation(primary_los.angularSeparation(other_los));
                if(cos_separation <= cos_max_ang) continue;
                double separation = std::acos(cos_separation);
                // accumulate statistics for pixel pairs
                for(auto& primary_pixel : primary_los.pixels) {
                    float los_pixel_dist_sq = primary_pixel.distance*primary_pixel.distance;
                    float los_pixel_projection_times_two = 2*primary_pixel.distance*cos_separation;
                    for(auto& other_pixel : other_los.pixels) {
                        // check pairs are within our binning grid
                        int pair_bin_index;
                        // scope for determining xi bin
                        {
                            // check parallel separation
                            float r_sq = los_pixel_dist_sq + (other_pixel.distance - los_pixel_projection_times_two)*other_pixel.distance;
                            float r = std::sqrt(r_sq);
                            // check transverse separation
                            // why is r == 0 sometimes? typically these cases also have
                            // very close separations.
                            float mu = (r == 0 ? 0 : std::fabs(primary_pixel.distance-other_pixel.distance)/r);
                            // check average pair distance
                            float z = 0.5*(primary_pixel.loglam + other_pixel.loglam) - local::logLyA;
                            if(coordinate_type_ == PolarCoordinates) {
                                if(r_sq >= rsq_max || r_sq < rsq_min) continue;
                                pair_bin_index = int((r - r_min)/r_spacing);
                                int mubin = int(mu * num_mu_bins);
                                if(mubin < 0) mubin = 0;
                                if(mubin >= num_mu_bins) mubin = num_mu_bins-1;
                                pair_bin_index = mubin + pair_bin_index*num_mu_bins;
                            }
                            else if(coordinate_type_ == CartesianCoordinates) {
                                // Cartesian Coordinates
                                float rparl = r*mu;
                                if(rparl >= r_max || rparl < r_min) continue;
                                float rperp = r*std::sqrt(1-mu*mu);
                                if(rperp >= rperp_max || rperp < rperp_min) continue;
                                pair_bin_index = int((rparl - r_min)/r_spacing);
                                int perp_bin_index = int((rperp - rperp_min)/rperp_spacing);
                                pair_bin_index = perp_bin_index + pair_bin_index*num_rperp_bins;
                            }
                            else {
                                // Observing coordinates
                                pair_bin_index = std::fabs(primary_pixel.loglam - other_pixel.loglam);
                                int sep_bin_index = int((separation - sep_min)/sep_spacing);
                                pair_bin_index = sep_bin_index + pair_bin_index*num_sep_bins;
                            }
                            int zbin = int((z-z_min)/z_spacing);
                            pair_bin_index = zbin + pair_bin_index*num_z_bins;
                        }
                        // this should never happen, useful for debugging though
                        if(pair_bin_index < 0 || pair_bin_index >= num_xi_bins) {
                            printf("%d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", pair_bin_index,
                                primary_pixel.distance, primary_los.ra, primary_los.dec,
                                other_pixel.distance, other_los.ra, other_los.dec,
                                cos_separation, separation);
                            throw local::RuntimeError("invalid bin index");
                        }
                        // accumulate pixel pair
                        float weight = primary_pixel.weight*other_pixel.weight;
                        float product = weight*primary_pixel.value*other_pixel.value;
                        xi[pair_bin_index].didj += product;
                        xi[pair_bin_index].wgt += weight;
                        // accumulate di, dj ???
                    }
                }
            }
        }
    }
    healxis_[id] = xi; // potentially not thread safe
    return true;
};

void local::XiEstimator::save_results(std::string outfile) {
    std::vector<double> xi_bin_centers(3);

    std::string estimator_filename(outfile + ".dat");
    std::cout << "Saving correlation function to: " << estimator_filename << std::endl;
    std::ofstream estimator_file(estimator_filename.c_str());
    for(int index = 0; index < xi_.size(); ++index) {
        grid_->getBinCenters(index, xi_bin_centers);
        estimator_file << index << ' ' << xi_bin_centers[0] << ' ' << xi_bin_centers[1] << ' ' << xi_bin_centers[2] << ' '
            << xi_[index].didj << ' ' << xi_[index].di << ' ' << xi_[index].dj << ' ' << xi_[index].wgt << std::endl;
    }
    estimator_file.close();

    std::string covariance_filename(outfile + ".cov");
    std::cout << "Saving covariance matrix to: " << covariance_filename << std::endl;
    std::ofstream covariance_file(covariance_filename.c_str());
    for(int a = 0; a < cov_.size(); ++a) {
        for(int b = 0; b < cov_[0].size(); ++b) {
            covariance_file << b + cov_[0].size()*a << ' ' << a << ' ' << b << ' ' << cov_[a][b] << std::endl;
        }
    }
    covariance_file.close();
};
