// Demo of ThreadPool. The iostream output is a bit garbled because we are not properly
// synchronizing access to stdout, but does demonstrate that the pool is operating correctly.

#include "turbooctospice.h"

#include "boost/thread.hpp"
#include "boost/bind.hpp"

#include <iostream>

namespace tos = turbooctospice;


bool task(int id) {
	std::cout << "starting task " << id << std::endl;
	boost::posix_time::seconds delay(2);
	boost::this_thread::sleep(delay);
	std::cout << "ending task " << id << std::endl;
	return id % 2;
}

int main(int argc, char **argv) {
	tos::ThreadPool pool(3);
	std::list<tos::ThreadPool::Task> tasks;
	for(int id = 0; id < 5; ++id) tasks.push_back(boost::bind(task,id));
	for(int batch = 0; batch < 3; ++batch) {
		std::cout << "Running batch " << batch << std::endl;
		pool.run(tasks);
	}
}
