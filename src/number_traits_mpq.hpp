#include <cassert>
#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {

	template<>
	const mpq_class NumberTraits<mpq_class>::zero = mpq_class(0);

	template<>
	const mpq_class NumberTraits<mpq_class>::one = mpq_class(1);

	template<>
	mpq_class NumberTraits<mpq_class>::abs(const mpq_class& value) {
		mpq_class result;
		mpq_abs(result.get_mpq_t(), value.get_mpq_t());

		return result;
	}

	template<>
	mpq_class NumberTraits<mpq_class>::absSqr(const mpq_class& value) {
		return value * value;
	}

	template<>
	mpq_class NumberTraits<mpq_class>::pow(const mpq_class& value, unsigned int n) {
		mpq_class result;
		mpz_pow_ui(mpq_numref(result.get_mpq_t()), mpq_numref(value.get_mpq_t()), n);
		mpz_pow_ui(mpq_denref(result.get_mpq_t()), mpq_denref(value.get_mpq_t()), n);

		return result;
	}

	template<>
	mpf_class NumberTraits<mpq_class>::sqrt(const mpq_class& value) {
		mpf_class sqrt_num, sqrt_den;
		mpf_sqrt(sqrt_num.get_mpf_t(), mpf_class(mpz_class(mpq_numref(value.get_mpq_t()))).get_mpf_t());
		mpf_sqrt(sqrt_den.get_mpf_t(), mpf_class(mpz_class(mpq_denref(value.get_mpq_t()))).get_mpf_t());

		return sqrt_num / sqrt_den;
	}

	template<>
	mpq_class NumberTraits<mpq_class>::div(const mpq_class& a, const mpq_class& b) {
		mpq_class result;
		mpq_div(result.get_mpq_t(), a.get_mpq_t(), b.get_mpq_t());

		return result;
	}

}
