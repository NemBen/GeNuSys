#include <math.h>
#include <complex>
#include <cassert>

namespace GeNuSys {

	template<>
	const double NumberTraits<double>::zero = 0.0;

	template<>
	const double NumberTraits<double>::one = 1.0;

	template<>
	const double NumberTraits<double>::epsilon = 0.000001;

	template<>
	double NumberTraits<double>::abs(const double& value) {
		return std::abs(value);
	}

	template<>
	double NumberTraits<double>::absSqr(const double& value) {
		return value * value;
	}

	template<>
	double NumberTraits<double>::sqrt(const double& value) {
		return std::sqrt(value);
	}

	template<>
	double NumberTraits<double>::div(const double& a, const double& b) {
		return a / b;
	}

}
