// Created 16-May-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "XiEstimator.h"
#include "AbsTwoPointGrid.h"
#include "RuntimeError.h"
#include "ThreadPool.h"
#include "SkyBins.h"

#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"

#include <fstream>

namespace local = turbooctospice;

local::XiEstimator::XiEstimator(double scale,
    local::AbsTwoPointGridPtr grid, local::XiEstimator::BinningCoordinateType type,
    std::vector<Forest> sightlines, SkyBinsIPtr skybins):
    grid_(grid), coordinate_type_(type),
    sightlines_(sightlines), skybins_(skybins),
    num_pixels_(0), num_sightline_pairs_(0), num_sightline_pairs_used_(0),
    num_pixel_pairs_(0), num_pixel_pairs_used_(0) {

    max_ang_ = grid_->maxAngularScale(scale);
    cos_max_ang_ = std::cos(max_ang_);

    // axis binning limits
    num_xi_bins_ = grid_->getNBinsTotal();
};

void local::XiEstimator::run(int nthreads) {
    local::ThreadPool pool(nthreads);
    std::list<local::ThreadPool::Task> tasks;

    // Estimate xi, one healpixel at a time
    std::cout << "Estimating xi..." << std::endl;
    auto occupied_bins = skybins_->getOccupiedBins();
    for(auto id : occupied_bins) {
        tasks.push_back(boost::bind(&XiEstimator::healxi_task, this, id));
    }
    show_progress_.reset(new boost::progress_display(tasks.size()));
    pool.run(tasks);
    tasks.clear();
    std::cout << std::endl;

    // Finalize xi
    std::cout << "Finalizing xi estimate..." << std::endl;
    xi_.resize(num_xi_bins_);
    for(int i = 0; i < num_xi_bins_; ++i) {
        tasks.push_back(boost::bind(&XiEstimator::xi_finalize_task, this, i));
    }
    show_progress_.reset(new boost::progress_display(tasks.size()));
    pool.run(tasks);
    tasks.clear();
    std::cout << "Xi estimation complete!" << std::endl;
    std::cout << std::endl;

    print_stats();

    // Estimate covariance
    std::cout << "Estimating covariance..." << std::endl;
    try {
        likely::CovarianceAccumulator cov_accum(num_xi_bins_);
        for(const auto& healxi_entry : healxis_) {
            std::vector<double> xi(num_xi_bins_);
            double weight(0);//Sky.getBinContents(healxi_entry.first).size());
            bool valid(true);
            for(int i = 0; i < num_xi_bins_; ++i) {
                if (healxi_entry.second[i].wgt > 0) {
                    xi[i] = healxi_entry.second[i].didj/healxi_entry.second[i].wgt;
                    weight += healxi_entry.second[i].wgt;
                }
            }
            cov_accum.accumulate(xi, weight);
        }
        cov_matrix_ = cov_accum.getCovariance();
        bool is_pos_def(cov_matrix_->isPositiveDefinite());
        std::cout << "is pos def: " << boost::lexical_cast<std::string>(is_pos_def) << std::endl;
        if(is_pos_def) {
            double detC = cov_matrix_->getLogDeterminant();
            std::cout << "log(|C|): " << boost::lexical_cast<std::string>(detC) << std::endl;
        }
    }
    catch(likely::RuntimeError const &e) {
        throw RuntimeError(e.what());
    }
    std::cout << "Covariance estimation complete!" << std::endl;
    std::cout << std::endl;

}

void local::XiEstimator::increment_progress() {
    boost::unique_lock<boost::mutex> scoped_lock(show_progress_mutex_);
    ++(*show_progress_);
    // show_progress_mutex_ is automatically released when lock
    // goes out of scope
}

void local::XiEstimator::accumulate_stats(unsigned long const &num_sightline_pairs,
    unsigned long const &num_sightline_pairs_used, unsigned long const &num_pixel_pairs,
    unsigned long const &num_pixel_pairs_used, unsigned long const &num_pixels) {
    boost::unique_lock<boost::mutex> scoped_lock(pair_stats_mutex_);
    num_sightline_pairs_ += num_sightline_pairs;
    num_sightline_pairs_used_ += num_sightline_pairs_used;
    num_pixel_pairs_ += num_pixel_pairs;
    num_pixel_pairs_used_ += num_pixel_pairs_used;
    num_pixels_ += num_pixels;
    // pair_stats_mutex_ is automatically released when lock
    // goes out of scope
}

bool local::XiEstimator::healxi_task(int id) {
    increment_progress();
    // create internal accumulation vectors
    std::vector<double> dsum(num_xi_bins_,0), wsum(num_xi_bins_,0);
    unsigned long num_sightline_pairs(0), num_sightline_pairs_used(0),
        num_pixel_pairs(0), num_pixel_pairs_used(0), num_pixels(0);
    // task_id's xi container
    std::vector<local::XiBin> xi(num_xi_bins_);

    const float
        r_min(grid_->getAxisMin(0)), r_max(grid_->getAxisMax(0)), one_over_dr(1.0/grid_->getAxisBinWidth(0)),
        rperp_min(grid_->getAxisMin(1)), rperp_max(grid_->getAxisMax(1)), one_over_drperp(1.0/grid_->getAxisBinWidth(1)),
        sep_min(grid_->getAxisMin(1)), sep_spacing(grid_->getAxisNBins(1)),
        z_min(grid_->getAxisMin(2)), z_max(grid_->getAxisMax(2)), z_spacing(grid_->getAxisBinWidth(2));
    const unsigned
        num_r_bins(grid_->getAxisNBins(0)),
        num_mu_bins(grid_->getAxisNBins(1)), num_rperp_bins(grid_->getAxisNBins(1)),
        num_sep_bins(grid_->getAxisNBins(1)),
        num_z_bins(grid_->getAxisNBins(2));

    const float rsq_min(r_min*r_min), rsq_max(r_max*r_max);

    // Iterate over all sight lines in this healpixel
    for(int primary_los_index : skybins_->getBinContents(id)) {
        auto primary_los = sightlines_[primary_los_index]; // should probably avoid copy
        num_pixels += primary_los.pixels.size();
        // Find healpixels within max angular separation of interest
        auto neighbors = skybins_->getBinIndicesWithinRadius(primary_los.ra, primary_los.dec, max_ang_);
        // Iterate over all neighboring healpixels
        for(int neighbor : neighbors) {
            // Check that if there are sightlines in the neighboring healpixel
            if(!skybins_->checkBinExists(neighbor)) continue;
            // Iterate over candidate sightline pairs
            for(int pair_los_index : skybins_->getBinContents(neighbor)) {
                // only count pairs once
                if(pair_los_index <= primary_los_index) continue;
                const auto other_los = sightlines_[pair_los_index];
                // check angular separation
                const double cos_separation(primary_los.angularSeparation(other_los));
                ++num_sightline_pairs;
                if(cos_separation <= cos_max_ang_) continue;
                ++num_sightline_pairs_used;
                const double separation = std::acos(cos_separation);
                // accumulate statistics for pixel pairs
                for(const auto& primary_pixel : primary_los.pixels) {
                    const float primary_pixel_dist_sq = primary_pixel.distance*primary_pixel.distance;
                    const float primary_pixel_projection_times_two = 2*primary_pixel.distance*cos_separation;
                    for(const auto& other_pixel : other_los.pixels) {
                        ++num_pixel_pairs;
                        // check pairs are within our binning grid
                        unsigned pair_bin_index;
                        // pair_bin_index = grid_.getBinIndex(primary_pixel, other_pixel, cos_separation, separation);
                        // scope for determining xi bin
                        {
                            // check parallel separation
                            const float r_sq = primary_pixel_dist_sq + (other_pixel.distance - primary_pixel_projection_times_two)*other_pixel.distance;
                            const float r = std::sqrt(r_sq);
                            // check transverse separation
                            // why is r == 0 sometimes? these cases also seem to have very close separations.
                            const float mu = (r == 0 ? 0 : std::fabs(primary_pixel.distance-other_pixel.distance)/r);
                            // check average pair distance
                            const float z = 0.5*(primary_pixel.loglam + other_pixel.loglam) - local::logLyA;
                            // const float z = std::pow(10, 0.5*(primary_pixel.loglam + other_pixel.loglam) - local::logLyA) - 1.0;
                            if(coordinate_type_ == PolarCoordinates) {
                                // if(r_sq >= rsq_max || r_sq < rsq_min) continue;
                                pair_bin_index = static_cast<unsigned>((r - r_min)*one_over_dr);
                                if(pair_bin_index >= num_r_bins) continue;
                                const unsigned mubin = static_cast<unsigned>(mu * num_mu_bins);
                                if(mubin >= num_mu_bins) continue;
                                pair_bin_index = mubin + pair_bin_index*num_mu_bins;
                            }
                            else if(coordinate_type_ == CartesianCoordinates) {
                                // Cartesian Coordinates
                                pair_bin_index = static_cast<unsigned>((r*mu - r_min)*one_over_dr);
                                if(pair_bin_index >= num_r_bins) continue;
                                const unsigned perp_bin_index = static_cast<unsigned>((r*std::sqrt(1-mu*mu) - rperp_min)*one_over_drperp);
                                if(perp_bin_index >= num_rperp_bins) continue;
                                pair_bin_index = perp_bin_index + pair_bin_index*num_rperp_bins;
                            }
                            else {
                                // Observing coordinates
                                // not finished
                                pair_bin_index = static_cast<unsigned>(std::fabs(primary_pixel.loglam - other_pixel.loglam)*one_over_dr);
                                const unsigned sep_bin_index = static_cast<unsigned>((separation - sep_min)/sep_spacing);
                                pair_bin_index = sep_bin_index + pair_bin_index*num_sep_bins;
                            }
                            const unsigned zbin = static_cast<unsigned>((z-z_min)/z_spacing);
                            pair_bin_index = zbin + pair_bin_index*num_z_bins;
                        }
                        // this should never happen, useful for debugging though
                        if(pair_bin_index >= num_xi_bins_) {
                            std::cout << boost::lexical_cast<std::string>(pair_bin_index) << " "
                                << boost::lexical_cast<std::string>(primary_pixel.distance) << " "
                                << boost::lexical_cast<std::string>(primary_pixel.loglam) << " "
                                << boost::lexical_cast<std::string>(primary_los.ra) << " "
                                << boost::lexical_cast<std::string>(primary_los.dec) << " "
                                << boost::lexical_cast<std::string>(other_pixel.distance) << " "
                                << boost::lexical_cast<std::string>(other_pixel.loglam) << " "
                                << boost::lexical_cast<std::string>(other_los.ra) << " "
                                << boost::lexical_cast<std::string>(other_los.dec) << " "
                                << boost::lexical_cast<std::string>(cos_separation) << " "
                                << boost::lexical_cast<std::string>(separation) << " "
                                << std::endl;
                            throw local::RuntimeError("invalid bin index");
                        }
                        // accumulate pixel pair
                        ++num_pixel_pairs_used;
                        xi[pair_bin_index].accumulate_pair(primary_pixel, other_pixel);
                    }
                }
            }
        }
    }
    healxis_[id] = xi; // thread safe as long as id is unique? do we need a mutex/lock here?
    accumulate_stats(num_sightline_pairs, num_sightline_pairs_used,
        num_pixel_pairs, num_pixel_pairs_used, num_pixels);
    return true;
};

bool local::XiEstimator::xi_finalize_task(int id) {
    increment_progress();
    for(const auto& healxi_entry : healxis_) {
        xi_[id] += healxi_entry.second[id];
    }
    xi_[id].finalize();
    return true;
}

void local::XiEstimator::save_results(std::string outfile) {
    std::string estimator_filename(outfile + ".data");
    std::cout << "Saving correlation function to: " << estimator_filename << std::endl;
    std::ofstream estimator_file(estimator_filename.c_str());
    std::vector<double> xi_bin_centers(3);
    for(int i = 0; i < num_xi_bins_; ++i) {
        grid_->getBinCenters(i, xi_bin_centers);
        estimator_file << i << ' ' << xi_bin_centers[0] << ' ' << xi_bin_centers[1] << ' ' << xi_bin_centers[2] << ' '
            << boost::lexical_cast<std::string>(xi_[i].didj) << ' '
            << boost::lexical_cast<std::string>(xi_[i].di) << ' '
            << boost::lexical_cast<std::string>(xi_[i].dj) << ' '
            << boost::lexical_cast<std::string>(xi_[i].wgt) << ' '
            << boost::lexical_cast<std::string>(xi_[i].num_pairs)
            << std::endl;
    }
    estimator_file.close();

    std::string covariance_filename(outfile + ".cov");
    std::cout << "Saving covariance matrix to: " << covariance_filename << std::endl;
    std::ofstream covariance_file(covariance_filename.c_str());
    for(int col = 0; col < num_xi_bins_; ++col) {
        for(int row = 0; row <= col; ++row) {
            double value = cov_matrix_->getCovariance(row,col);
            // print matrix elements with full precision
            covariance_file << row << ' ' << col << ' '
                << boost::lexical_cast<std::string>(value) << std::endl;
        }
    }
    covariance_file.close();

    std::string weight_filename(outfile + ".wgt");
    std::cout << "Saving weights matrix to: " << weight_filename << std::endl;
    std::ofstream weight_file(weight_filename.c_str());
    for(int col = 0; col < num_xi_bins_; ++col) {
            double value = 1.0/xi_[col].wgt;
            // print matrix elements with full precision
            weight_file << col << ' ' << col << ' '
                << boost::lexical_cast<std::string>(value) << std::endl;
    }
    weight_file.close();
};

void local::XiEstimator::save_subsamples(std::string outfile_base) {
    for(const auto& healxi_entry : healxis_) {
        std::string estimator_filename(outfile_base + "-" + boost::lexical_cast<std::string>(healxi_entry.first) + ".data");
        std::ofstream estimator_file(estimator_filename.c_str());
        for(int i = 0; i < num_xi_bins_; ++i) {
            estimator_file << i << ' ' << boost::lexical_cast<std::string>(healxi_entry.second[i].didj) << std::endl;
        }
        estimator_file.close();
    }
}

void local::XiEstimator::print_stats() {
    // line of sight pair statistics
    unsigned num_sightlines = sightlines_.size();
    unsigned long num_sightline_pairs_total = (num_sightlines*(num_sightlines-1))/2;
    double frac_sightline_pairs_considered = static_cast<double>(num_sightline_pairs_)/num_sightline_pairs_total;
    double frac_sightline_pairs_used = static_cast<double>(num_sightline_pairs_used_)/num_sightline_pairs_;

    std::cout << "Number of distinct los pairs " << num_sightline_pairs_total << std::endl;
    std::cout << "considered " << num_sightline_pairs_ << " of distinct los pairs. (" << frac_sightline_pairs_considered << ")" << std::endl;
    std::cout << "used " << num_sightline_pairs_used_ << " of los pairs considered. (" << frac_sightline_pairs_used << ")" << std::endl;

    // pixel pair statistics
    unsigned long num_pixel_pairs_total = (num_pixels_*(num_pixels_-1))/2;
    double frac_pixel_pairs_considered = static_cast<double>(num_pixel_pairs_)/num_pixel_pairs_total;
    double frac_pixel_pairs_used = static_cast<double>(num_pixel_pairs_used_)/num_pixel_pairs_;

    std::cout << "Number of distinct pixel pairs " << num_pixel_pairs_total << std::endl;
    std::cout << "considered " << num_pixel_pairs_ << " of distinct pixel pairs. (" << frac_pixel_pairs_considered << ")" << std::endl;
    std::cout << "used " << num_pixel_pairs_used_ << " of pixel pairs considered. (" << frac_pixel_pairs_used << ")" << std::endl;
}
