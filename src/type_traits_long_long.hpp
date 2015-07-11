#include <complex>

namespace GeNuSys {
	
	template<>
	struct TypeTraitsBase<long long> {

		typedef double RealType;

		typedef double RationalType;

		typedef std::complex<long long> ComplexType;

		typedef long long AbsType;

		typedef long long AbsSqrType;

	};

	template<>
	template<>
	long long TypeTraits<long long>::asType<long long>(const long long& value) {
		return value;
	}

	template<>
	template<>
	double TypeTraits<long long>::asType<double>(const long long& value) {
		return (double) value;
	}

	template<>
	template<>
	std::complex<long long> TypeTraits<long long>::asType<std::complex<long long> >(const long long& value) {
		return std::complex<long long>(value, 0);
	}

	template<>
	template<>
	std::complex<double> TypeTraits<long long>::asType<std::complex<double> >(const long long& value) {
		return std::complex<double>((double) value, 0);
	}

	template<>
	template<>
	int TypeTraits<long long>::asTypeUnsafe<int>(const long long& value) {
		return (int) value;
	}

}
