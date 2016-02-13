#ifndef TURBOOCTOSPICE_THREAD_POOL
#define TURBOOCTOSPICE_THREAD_POOL

#include <list>

#include <boost/smart_ptr.hpp>
#include <boost/function.hpp>

namespace turbooctospice {
	class ThreadPool {
	public:
		/// Creates a new thread pool with the specified number of threads.
		/// @param nthreads Number of threads to use.
		ThreadPool(int nthreads);
		virtual ~ThreadPool();
		/// A task takes no arguments and returns a status bool.
		typedef boost::function<bool ()> Task;
		/// Runs the list of tasks provided and returns the logical AND of all status bools
		/// when they are all complete.
		// @param tasks List of tasks to process.
		bool run(std::list<Task> tasks);
	private:
		class Implementation;
		boost::scoped_ptr<Implementation> _pimpl;
	}; // ThreadPool
} // turbooctospice

#endif // TURBOOCTOSPICE_THREAD_POOL
