#include <iostream>


#include "boost/bind.hpp"
#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include "boost/coroutine/coroutine.hpp"

#include <vector>

template <typename T>
	class Pair : public std::pair<T,T> {
	public:
		typedef std::pair<T,T> Base;
		typedef boost::coroutines::coroutine<Pair> PairGenerator;
		Pair() : Base() {};
		Pair(const Pair &pair) : Base(pair) {};
		Pair(T first, T second) : Base(first,second) {};
	private:
}; // Pair

typedef Pair<int> PairType;

typedef boost::coroutines::coroutine<PairType> coro_t;

class SearchPolicy {
public:
	SearchPolicy(std::vector<int> numbers) {
		_numbers = numbers;
	}
	void findInts(coro_t::push_type &yield) const {
		BOOST_FOREACH(int i, _numbers) {
			BOOST_FOREACH(int j, _numbers) {
				yield(Pair<int>(i,j));
			}
		}
	}
private:
	std::vector<int> _numbers;
};

int main(int argc, char **argv) {

	std::vector<int> myInts;

	myInts.push_back(1);
	myInts.push_back(2);
	myInts.push_back(3);

	boost::shared_ptr<const SearchPolicy> sp(new SearchPolicy(myInts));

	coro_t::pull_type generator(boost::bind(&SearchPolicy::findInts, sp, _1));

	BOOST_FOREACH(PairType& i, generator){
		std::cout << i.first << " " << i.second << std::endl;
	}

	std::cout << "hello!" << std::endl;
	return 0;
}