#ifndef TOS_MOCK_SPECTRUM
#define TOS_MOCK_SPECTRUM

#include <vector>
#include <string>

namespace turbooctospice {

	struct QuasarPixel {
    	float frac, lam, wgt, dist;
	};

	class MockSpectrum {

	public:
	    MockSpectrum(std::string target, bool verbose=false);

	    float getZ();
	    float getRA();
	    float getDec();
	    float getCoeff0();
	    float getCoeff1();

	    std::vector<QuasarPixel> getTrimmedSpectrum(
	    	int ncombine=1, float forestlo=1040, float foresthi=1200, float speclo=3650);
	private:
		void loadTarget(bool verbose);
	    float _z, _ra, _dec, _coeff0, _coeff1;
	    std::string _target;
	    std::vector<float> _frac;
	};

	inline float MockSpectrum::getZ() { return _z; };
	inline float MockSpectrum::getRA() { return _ra; };
	inline float MockSpectrum::getDec() { return _dec; };
	inline float MockSpectrum::getCoeff0() { return _coeff0; };
	inline float MockSpectrum::getCoeff1() { return _coeff1; };

	std::string getMockFilename(std::string target);
}

#endif
