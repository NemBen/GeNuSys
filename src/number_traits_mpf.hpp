#include <complex>

#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {

	template<>
	const mpf_class NumberTraits<mpf_class>::zero = mpf_class(0);

	template<>
	const mpf_class NumberTraits<mpf_class>::one = mpf_class(1);
	
	template<>
	const mpf_class NumberTraits<mpf_class>::epsilon = mpf_class(NumberTraits<double>::epsilon);

	template<>
	mpf_class NumberTraits<mpf_class>::abs(const mpf_class& value) {
		mpf_class result;
		mpf_abs(result.get_mpf_t(), value.get_mpf_t());

		return result;
	}

	template<>
	mpf_class NumberTraits<mpf_class>::absSqr(const mpf_class& value) {
		return value * value;
	}

	template<>
	mpf_class NumberTraits<mpf_class>::pow(const mpf_class& value, unsigned int n) {
		mpf_class result;
		mpf_pow_ui(result.get_mpf_t(), value.get_mpf_t(), n);

		return result;
	}

	template<>
	mpf_class NumberTraits<mpf_class>::sqrt(const mpf_class& value) {
		mpf_class result;
		mpf_sqrt(result.get_mpf_t(), value.get_mpf_t());

		return result;
	}

	template<>
	mpf_class NumberTraits<mpf_class>::div(const mpf_class& a, const mpf_class& b) {
		mpf_class result;
		mpf_div(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());

		return result;
	}

}
