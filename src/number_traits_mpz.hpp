#include <cassert>
#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {

	template<>
	const mpz_class NumberTraits<mpz_class>::zero = mpz_class(0);

	template<>
	const mpz_class NumberTraits<mpz_class>::one = mpz_class(1);

	template<>
	mpz_class NumberTraits<mpz_class>::abs(const mpz_class& value) {
		mpz_class result;
		mpz_abs(result.get_mpz_t(), value.get_mpz_t());

		return result;
	}

	template<>
	mpz_class NumberTraits<mpz_class>::absSqr(const mpz_class& value) {
		return value * value;
	}

	template<>
	mpz_class NumberTraits<mpz_class>::pow(const mpz_class& value, unsigned int n) {
		mpz_class result;
		mpz_pow_ui(result.get_mpz_t(), value.get_mpz_t(), n);

		return result;
	}

	template<>
	mpf_class NumberTraits<mpz_class>::sqrt(const mpz_class& value) {
		mpf_class result;
		mpf_sqrt(result.get_mpf_t(), mpf_class(value).get_mpf_t());

		return result;
	}

	template<>
	mpq_class NumberTraits<mpz_class>::div(const mpz_class& a, const mpz_class& b) {
		mpq_class result;
		mpq_set_num(result.get_mpq_t(), a.get_mpz_t());
		mpq_set_den(result.get_mpq_t(), b.get_mpz_t());
		mpq_canonicalize(result.get_mpq_t());

		return result;
	}

}
