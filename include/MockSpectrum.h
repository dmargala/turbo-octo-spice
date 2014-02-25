#ifndef TOS_MOCK_SPECTRUM
#define TOS_MOCK_SPECTRUM

#include <vector>
#include <string>

namespace turbooctospice {

	struct QuasarPixel {
    	float flux, lam, wgt, dist;
	};

	class MockSpectrum {

	public:
	    MockSpectrum(float z, float ra, float dec);

	    float getZ();
	    float getRA();
	    float getDec();

	    std::vector<QuasarPixel> pixels;
	private:
	    float _z, _ra, _dec;
	};

	inline float MockSpectrum::getZ() { return _z; };
	inline float MockSpectrum::getRA() { return _ra; };
	inline float MockSpectrum::getDec() { return _dec; };

	std::string getMockFilename(std::string target);
	void readMockTargets(std::vector<std::string> &targetlist, 
		std::vector<MockSpectrum> &quasars, bool const verbose);

}

#endif
