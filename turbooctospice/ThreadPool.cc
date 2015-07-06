// Created 16-May-2014 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "ThreadPool.h"
#include "RuntimeError.h"

#include <boost/bind.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/future.hpp>

#include <list>
#include <iostream>

namespace local = turbooctospice;

class local::ThreadPool::Implementation {
public:
	Implementation(int nthreads,bool verbose=false) : _verbose(verbose) {
		if(_verbose) std::cout << "Starting " << nthreads << " threads." << std::endl;
		_work.reset(new boost::asio::io_service::work(_service));
		for(int i = 0; i < nthreads; ++i) {
			_pool.create_thread(boost::bind(&boost::asio::io_service::run, &_service));
		}
	}
	~Implementation() {
		if(_verbose) std::cout << "Stopping the thread pool service." << std::endl;
		_work.reset();
		if(_verbose) std::cout << "Waiting for threads to complete." << std::endl;
		_pool.join_all();
		if(_verbose) std::cout << "Thread pool is now closed." << std::endl;
	}
	bool run(std::list<Task> tasks) {
		// Submit each task to our thread pool, keeping pointers to their future results.
		// http://stackoverflow.com/questions/13157502/how-do-you-post-a-boost-packaged-task-to-an-io-service-in-c03
		typedef boost::shared_future<bool> future;
		std::vector<future> pending;
		for(std::list<Task>::iterator iter = tasks.begin(); iter != tasks.end(); ++iter) {
			typedef boost::packaged_task<bool> task_t;
			boost::shared_ptr<task_t> task = boost::make_shared<task_t>(*iter);
			boost::shared_future<bool> result(task->get_future());
			pending.push_back(result);
			_service.post(boost::bind(&task_t::operator(),task));
		}
		// Wait for all of the tasks to complete.
		boost::wait_for_all(pending.begin(),pending.end());
		// Return the logical AND of the individual task results.
		bool result(true);
		for(std::vector<future>::iterator iter = pending.begin(); iter != pending.end(); ++iter) {
			result &= iter->get();
		}
		return result;
	}
private:
	bool _verbose;
	boost::asio::io_service _service;
	boost::thread_group _pool;
	boost::scoped_ptr<boost::asio::io_service::work> _work;
};

local::ThreadPool::ThreadPool(int nthreads) : _pimpl(new Implementation(nthreads,true)) { }

local::ThreadPool::~ThreadPool() { }

bool local::ThreadPool::run(std::list<Task> tasks) {
	return _pimpl->run(tasks);
}
