#include <cmath>
#include <complex>
#include <cassert>

namespace GeNuSys {

	template<>
	const int NumberTraits<int>::zero = 0;

	template<>
	const int NumberTraits<int>::one = 1;

	template<>
	int NumberTraits<int>::abs(const int& value) {
		return std::abs(value);
	}

	template<>
	int NumberTraits<int>::absSqr(const int& value) {
		return value * value;
	}

	template<>
	double NumberTraits<int>::sqrt(const int& value) {
		return std::sqrt((double) value);
	}

	template<>
	double NumberTraits<int>::div(const int& a, const int& b) {
		return (double) a / b;
	}

}
