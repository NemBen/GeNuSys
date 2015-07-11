#include <complex>

namespace GeNuSys {

	template<>
	struct TypeTraitsBase<double> {

		typedef double RealType;

		typedef double RationalType;

		typedef std::complex<double> ComplexType;

		typedef double AbsType;

		typedef double AbsSqrType;

	};

	template<>
	template<>
	double TypeTraits<double>::asType<double>(const double& value) {
		return value;
	}

	template<>
	template<>
	std::complex<double> TypeTraits<double>::asType<std::complex<double> >(const double& value) {
		return std::complex<double>(value, 0);
	}

	template<>
	template<>
	int TypeTraits<double>::asTypeUnsafe<int>(const double& value) {
		return (int)floor(value + 0.5);
	}

	template<>
	template<>
	long long TypeTraits<double>::asTypeUnsafe<long long>(const double& value) {
		return (long long)floor(value + 0.5);
	}

}
