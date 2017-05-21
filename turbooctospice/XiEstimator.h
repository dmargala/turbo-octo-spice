#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <map>

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/thread/mutex.hpp>

#include "likely/likely.h"

#include "types.h"
#include "constants.h"

namespace turbooctospice {
    /// The XiEstimator
    class XiEstimator {
    public:
        enum BinningCoordinateType {
            PolarCoordinates, CartesianCoordinates, ObservingCoordinates
        };
        /// Create a XiEstimator object
        /// @param scale The transverse comoving scale at the minimum redshift.
        /// @param grid The binning grid to use.
        /// @param type The binning coordinate type.
        /// @param sightlines The forest sightlines to use for calculation.
        /// @param skybins The skybinning bins.
        XiEstimator(double scale, AbsTwoPointGridPtr grid,
            BinningCoordinateType type, std::vector<Forest> sightlines, SkyBinsIPtr skybins);
        /// Run estimator.
        /// @param nthreads The number of threads to use for calcualtion.
        void run(int nthreads, bool do_cov=false);
        /// Save results to the specified outfile.
        /// @param outfile The filename to save results to.
        void save_results(std::string outfile) const;
        /// Save individual sky bin subsamples to files.
        /// @param outfile_base The filename base to use when saving individual filenames.
        void save_subsamples(std::string outfile_base) const;
        /// Print summary of pixel and pair statistics.
        void print_stats() const;
    private:
        /// Increment progress bar.
        void increment_progress();
        /// Process individual sky bin.
        /// @param skybin_index The index of the skybin.
        bool skybin_xi_task(int skybin_index);
        /// Accumulate correlation between two sightlines.
        /// @param primary_los The primary line of sight.
        /// @param other_los The other line of sight.
        /// @param cos_separation The cosine of the angular seperation between the primary and other line of sight.
        /// @param xi Correlation function bins.
        unsigned long accumulate_pixel_pairs(const Forest &primary_los, const Forest &other_los,
            const double &cos_separation, std::vector<XiBin> &xi) const;
        /// Finalize xi bin by accumulating results accross skybins.
        /// @param xi_bin_index The index of the correlation function bin.
        bool xi_finalize_task(int xi_bin_index);
        /// Accumulate stats from individual sky bins.
        /// @param num_sightline_pairs The total number of sightline pairs.
        /// @param num_sightline_pairs_used The number of sightline pairs that contribute to correlation function.
        /// @param num_pixel_pairs The total number of pixel pairs.
        /// @param num_pixel_pairs_used The number of pixel pairs that contribute to correlation function.
        /// @param num_pixels The total number of pixels.
        void accumulate_stats(unsigned long const &num_sightline_pairs,
            unsigned long const &num_sightline_pairs_used, unsigned long const &num_pixel_pairs,
            unsigned long const &num_pixel_pairs_used, unsigned long const &num_pixels);
        // private member variables
        // counters
        unsigned long num_sightline_pairs_, num_sightline_pairs_used_,
            num_pixels_, num_pixel_pairs_, num_pixel_pairs_used_;
        // inputs
        double max_ang_, cos_max_ang_;
        const BinningCoordinateType coordinate_type_;
        const AbsTwoPointGridPtr grid_;
        const SkyBinsIPtr skybins_;
        const std::vector<Forest> sightlines_;
        const unsigned num_xi_bins_;
        // results
        std::map<int, std::vector<XiBin> > skybin_xis_;
        std::vector<XiBin> xi_;
        likely::CovarianceMatrixCPtr cov_matrix_;
        bool cov_good_;
        // helpers
        std::unique_ptr<boost::progress_display> show_progress_;
        boost::mutex show_progress_mutex_; // protects show_progress_
        boost::mutex pair_stats_mutex_; // protects counters
    };
} // turbooctospice

#endif // TURBOOCTOSPICE_XI_ESTIMATOR
