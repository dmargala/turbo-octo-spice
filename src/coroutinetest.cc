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

class my_generator {
public:
	typedef int result_type;
	my_generator() : state(0) { }
	int operator()() { return ++state; }
private:
	int state;
};

template <typename T> class Generator {
	public:
	Generator(std::vector<T> numbers) : _numbers(numbers) {
		i = j = 0;
		std::cout << "Generator constructor" << std::endl;
	}
	~Generator(){
		std::cout << "Generator destructor" << std::endl;
	}
    bool valid() {
    	if (i == _numbers.size()) {
    		return false;
    	}
    	else {
    		return true;
    	}
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

	myInts.push_back(1);
	myInts.push_back(2);
	myInts.push_back(3);

	boost::shared_ptr<const SearchPolicy> sp(new SearchPolicy(myInts));

	coro_t::pull_type generator(boost::bind(&SearchPolicy::findInts, sp, _1));

	BOOST_FOREACH(PairType& i, generator){
		std::cout << i.first << " " << i.second << std::endl;
	}

	my_generator gen;
	boost::generator_iterator_generator<my_generator>::type it = boost::make_generator_iterator(gen);
	for(int i = 0; i < 10; ++i, ++it){
		std::cout << *it << std::endl;
	}

	Generator<int> mygen(myInts);
	while(mygen.valid()) {
		PairType i(mygen.get());
		std::cout << i.first << " " << i.second << std::endl;
		mygen.next();
	}

	std::cout << "hello!" << std::endl;
	return 0;
}