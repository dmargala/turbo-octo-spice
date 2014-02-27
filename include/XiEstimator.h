#ifndef TOS_XI_ESTIMATOR
#define TOS_XI_ESTIMATOR

#include <iostream>
#include <string>
#include <vector>

namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : private PairSearchPolicy, private BinPolicy {
        using PairSearchPolicy::findPairs;
        using BinPolicy::binPair;
    public:
        void run(std::vector<float> &pixels) const {
        	findPairs(pixels, binPair);
        }
    };
     
    class PairSearchPolicyBrute {
    protected:
        template<typename BinMethod> void findPairs(std::vector<float> &pixels, BinMethod binPair) const {
            for(int i = 0; i < pixels.size()-1; ++i) {
                float a = pixels[i];
                for(int j = i+1; j < pixels.size(); ++j) {
                    float b = pixels[j];
                    std::cout << a << "," << b << " -> " << binPair(a, b) << std::endl;
                }
            }
        }
    };
     
    class BinPolicyDummy {
    protected:
        static float binPair(float a, float b);
    };

    class BinPolicyWeighted {
    protected:
        static float binPair(float a, float b);
    };

}

#endif