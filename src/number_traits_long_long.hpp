#include <cmath>
#include <complex>
#include <cassert>

namespace GeNuSys {

	template<>
	const long long NumberTraits<long long>::zero = 0;

	template<>
	const long long NumberTraits<long long>::one = 1;

	template<>
	long long NumberTraits<long long>::abs(const long long& value) {
		return std::abs(value);
	}

	template<>
	long long NumberTraits<long long>::absSqr(const long long& value) {
		return value * value;
	}

	template<>
	double NumberTraits<long long>::sqrt(const long long& value) {
		return std::sqrt((double) value);
	}
	
	template<>
	double NumberTraits<long long>::div(const long long& a, const long long& b) {
		return (double) a / b;
	}

}
