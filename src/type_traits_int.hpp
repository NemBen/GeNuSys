#include <complex>

namespace GeNuSys {
	
	template<>
	struct TypeTraitsBase<int> {

		typedef double RealType;

		typedef double RationalType;

		typedef std::complex<int> ComplexType;

		typedef int AbsType;

		typedef int AbsSqrType;

	};

	template<>
	template<>
	int TypeTraits<int>::asType<int>(const int& value) {
		return value;
	}

	template<>
	template<>
	long long TypeTraits<int>::asType<long long>(const int& value) {
		return value;
	}

	template<>
	template<>
	double TypeTraits<int>::asType<double>(const int& value) {
		return (double) value;
	}

	template<>
	template<>
	std::complex<int> TypeTraits<int>::asType<std::complex<int> >(const int& value) {
		return std::complex<int>(value, 0);
	}

	template<>
	template<>
	std::complex<long long> TypeTraits<int>::asType<std::complex<long long> >(const int& value) {
		return std::complex<long long>(value, 0);
	}

	template<>
	template<>
	std::complex<double> TypeTraits<int>::asType<std::complex<double> >(const int& value) {
		return std::complex<double>((double) value, 0);
	}

}
