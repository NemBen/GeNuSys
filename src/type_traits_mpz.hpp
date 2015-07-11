#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {
	
	template<>
	struct TypeTraitsBase<mpz_class> {

		typedef mpf_class RealType;

		typedef mpq_class RationalType;

		typedef std::complex<mpz_class> ComplexType;

		typedef mpz_class AbsType;

		typedef mpz_class AbsSqrType;

	};

	template<>
	template<>
	mpz_class TypeTraits<mpz_class>::asType<mpz_class>(const mpz_class& value) {
		return value;
	}

	template<>
	template<>
	mpf_class TypeTraits<mpz_class>::asType<mpf_class>(const mpz_class& value) {
		return mpf_class(value);
	}

	template<>
	template<>
	mpq_class TypeTraits<mpz_class>::asType<mpq_class>(const mpz_class& value) {
		return mpq_class(value);
	}

	template<>
	template<>
	std::complex<mpz_class> TypeTraits<mpz_class>::asType<std::complex<mpz_class> >(const mpz_class& value) {
		return std::complex<mpz_class>(value, mpz_class(0));
	}

	template<>
	template<>
	std::complex<mpf_class> TypeTraits<mpz_class>::asType<std::complex<mpf_class> >(const mpz_class& value) {
		return std::complex<mpf_class>(value, mpz_class(0));
	}

	template<>
	template<>
	std::complex<mpq_class> TypeTraits<mpz_class>::asType<std::complex<mpq_class> >(const mpz_class& value) {
		return std::complex<mpq_class>(value, mpz_class(0));
	}

	template<>
	template<>
	mpz_class TypeTraits<int>::asType<mpz_class>(const int& value) {
		return mpz_class(value);
	}

	template<>
	template<>
	std::complex<mpz_class> TypeTraits<int>::asType<std::complex<mpz_class> >(const int& value) {
		return std::complex<mpz_class>(mpz_class(value), mpz_class(0));
	}

}
