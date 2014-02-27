#ifndef TOS_XI_ESTIMATOR
#define TOS_XI_ESTIMATOR

#include <iostream>
#include <string>

namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : private PairSearchPolicy, private BinPolicy {
        using PairSearchPolicy::print;
        using BinPolicy::message;
    public:
        void run() const {
        	print(message());
        }
    };
     
    class PairSearchPolicyBrute {
    protected:
        template<typename MessageType> void print(MessageType const &message) const {
        std::cout << message << std::endl;
        }
    };
     
    class BinPolicyDummy {
    protected:
        std::string message() const;
    };

    class BinPolicyWeighted {
    protected:
        std::string message() const;
    };

}

#endif