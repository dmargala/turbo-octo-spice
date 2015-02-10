#include <iostream>


#include "boost/bind.hpp"
#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include "boost/coroutine/coroutine.hpp"
#include "boost/generator_iterator.hpp"

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

template <typename T> class BoostGenerator {
public:
	typedef Pair<T> result_type;
	BoostGenerator(std::vector<T> numbers, int _i=0) : 
		_numbers(numbers), i(_i), j(0) { }
	Pair<T> operator()() { 
		Pair<T> return_value(_numbers[i], _numbers[j]);
		next();
		return return_value; 
	}
	bool valid() {
		if (i == _numbers.size()) return false;
    	return true;
	}
	void next() {
		++j;
		if (j == _numbers.size()) {
			++i;
			j = 0;
		}
	}
	BoostGenerator begin() {
		return BoostGenerator(*this, _numbers.size());
	}
	BoostGenerator end() {
		return BoostGenerator(*this, _numbers.size());
	}
private:
	int i, j;
	std::vector<T> _numbers;
};

template <typename T> class SimpleGenerator {
	public:
	SimpleGenerator(std::vector<T> numbers) : 
		_numbers(numbers), i(0), j(0) {
	}
    bool valid() {
    	if (i == _numbers.size()) return false;
    	return true;
    }
    void next() {
    	++j;
    	if (j == _numbers.size()) {
    		++i;
    		j = 0;
    	}
    }
    Pair<T> get() {
    	return Pair<T>(_numbers[i],_numbers[j]);
    }
    std::vector<T> _numbers;
    int i, j;
};

int main(int argc, char **argv) {

	std::vector<int> myInts;

	for(int i = 0; i < 10000; ++i){
		myInts.push_back(i);
	}

	// myInts.push_back(1);
	// myInts.push_back(2);
	// myInts.push_back(3);

	std::cout << "Testing coroutine generator..." << std::endl;

	boost::shared_ptr<const SearchPolicy> sp(new SearchPolicy(myInts));
	coro_t::pull_type generator(boost::bind(&SearchPolicy::findInts, sp, _1));
	BOOST_FOREACH(PairType& i, generator){
		//std::cout << i.first << " " << i.second << std::endl;
	}

	// std::cout << "Testing boost generator..." << std::endl;

	// BoostGenerator<int> boostGen(myInts);
	// boost::generator_iterator_generator<BoostGenerator<int> >::type it = boost::make_generator_iterator(boostGen);
	// BOOST_FOREACH(PairType& i, it){
	// //for(int i = 0; i < 10; ++i, ++it){
	// 	std::cout << i.first << " " << i.second << std::endl;
	// }

	std::cout << "Testing simple generator..." << std::endl;

	SimpleGenerator<int> simpleGen(myInts);

	while(simpleGen.valid()) {
		PairType i(simpleGen.get());
		//std::cout << i.first << " " << i.second << std::endl;
		simpleGen.next();
	}

	return 0;
}