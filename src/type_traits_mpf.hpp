#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {
	
	template<>
	struct TypeTraitsBase<mpf_class> {

		typedef mpf_class RealType;

		typedef mpf_class RationalType;

		typedef std::complex<mpf_class> ComplexType;

		typedef mpf_class AbsType;

		typedef mpf_class AbsSqrType;

	};

	template<>
	template<>
	mpf_class TypeTraits<mpf_class>::asType<mpf_class>(const mpf_class& value) {
		return value;
	}

	template<>
	template<>
	std::complex<mpf_class> TypeTraits<mpf_class>::asType<std::complex<mpf_class> >(const mpf_class& value) {
		return std::complex<mpf_class>(value, mpf_class(0));
	}

	template<>
	template<>
	mpf_class TypeTraits<int>::asType<mpf_class>(const int& value) {
		return mpf_class(value);
	}
	
	template<>
	template<>
	std::complex<mpf_class> TypeTraits<int>::asType<std::complex<mpf_class> >(const int& value) {
		return std::complex<mpf_class>(mpf_class(value), mpf_class(0));
	}

}
