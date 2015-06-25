#ifndef TURBOOCTOSPICE_RUNTIME_ERROR
#define TURBOOCTOSPICE_RUNTIME_ERROR

#include <stdexcept>
#include <string>

namespace turbooctospice {
	class RuntimeError : public std::runtime_error {
	public:
		explicit RuntimeError(std::string const &reason);
		virtual ~RuntimeError() throw ();
	private:
	}; // RuntimeError
	
	inline RuntimeError::RuntimeError(std::string const &reason)
	: std::runtime_error(reason) { }

    inline RuntimeError::~RuntimeError() throw () { }

} // turbooctospice

#endif // TURBOOCTOSPICE_RUNTIME_ERROR