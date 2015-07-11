#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {

	template<>
	struct TypeTraitsBase<mpq_class> {

		typedef mpf_class RealType;

		typedef mpq_class RationalType;

		typedef std::complex<mpq_class> ComplexType;

		typedef mpq_class AbsType;

		typedef mpq_class AbsSqrType;

	};

	template<>
	template<>
	mpq_class TypeTraits<mpq_class>::asType<mpq_class>(const mpq_class& value) {
		return value;
	}

	template<>
	template<>
	mpf_class TypeTraits<mpq_class>::asType<mpf_class>(const mpq_class& value) {
		return mpf_class(value);
	}

	template<>
	template<>
	std::complex<mpq_class> TypeTraits<mpq_class>::asType<std::complex<mpq_class> >(const mpq_class& value) {
		return std::complex<mpq_class>(value, mpq_class(0));
	}

	template<>
	template<>
	std::complex<mpf_class> TypeTraits<mpq_class>::asType<std::complex<mpf_class> >(const mpq_class& value) {
		return std::complex<mpf_class>(value, mpq_class(0));
	}

	template<>
	template<>
	mpz_class TypeTraits<mpq_class>::asTypeUnsafe<mpz_class>(const mpq_class& value) {
		mpz_class result;
		mpz_fdiv_q(result.get_mpz_t(), value.get_num_mpz_t(), value.get_den_mpz_t());

		return result;
	}

	template<>
	template<>
	mpq_class TypeTraits<int>::asType<mpq_class>(const int& value) {
		return mpq_class(value);
	}

	template<>
	template<>
	std::complex<mpq_class> TypeTraits<int>::asType<std::complex<mpq_class> >(const int& value) {
		return std::complex<mpq_class>(value, mpq_class(0));
	}

}
