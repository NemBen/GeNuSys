#include <gmp.h>
#include <gmpxx.h>

namespace GeNuSys {

	template<>
	bool IntegerTraits<mpz_class>::divisible(const mpz_class& a, const mpz_class& b) {
		return mpz_divisible_p(a.get_mpz_t(), b.get_mpz_t()) != 0;
	}

	template<>
	mpz_class IntegerTraits<mpz_class>::div(const mpz_class& a, const mpz_class& b) {
		mpz_class result;
		mpz_fdiv_q(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

		return result;
	}

	template<>
	mpz_class IntegerTraits<mpz_class>::mod(const mpz_class& a, const mpz_class& b) {
		mpz_class result;
		mpz_mod(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

		return result;
	}

	template<>
	mpz_class IntegerTraits<mpz_class>::mods(const mpz_class& a, const mpz_class& b) {
		mpz_class result;
		mpz_mod(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

		return (result > b / 2) ? result - b : result;
	}

	template<>
	ExtendedGCD<mpz_class> IntegerTraits<mpz_class>::egcd(const mpz_class& a, const mpz_class& b) {
		mpz_class gcd, aC, bC;
		mpz_gcdext(gcd.get_mpz_t(), aC.get_mpz_t(), bC.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

		return ExtendedGCD<mpz_class>(a, aC, b, bC, gcd);
	}

}
