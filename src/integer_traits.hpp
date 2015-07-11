#include <cassert>

namespace GeNuSys {

	template<typename Type>
	ExtendedGCD<Type>::ExtendedGCD(
		const Type& a,
		const Type& cA,
		const Type& b,
		const Type& cB,
		const Type& gcd):a(a),cA(cA),b(b),cB(cB),gcd(gcd) {
	}

	template<typename Type>
	Type IntegerTraits<Type>::gcd(Type a, Type b) {
		while (b != 0) {
			Type tmp = b;
			b = IntegerTraits<Type>::mod(a,b);
			a = tmp;
		}

		return a;
	}

	template<typename Type>
	ExtendedGCD<Type> IntegerTraits<Type>::egcd(const Type& a, const Type& b) {
		Type s = 0, old_s = 1;
		Type t = 1, old_t = 0;
		Type r = a, old_r = b;
		Type tmp, q;
		while (r != 0) {
			q = old_r / r;

			tmp = r;
			r = old_r - q * r;
			old_r = tmp;

			tmp = s;
			s = old_s - q * s;
			old_s = tmp;

			tmp = t;
			t = old_t - q * t;
			old_t = tmp;
		}

		assert(a * old_t + b * old_s == old_r);

		return ExtendedGCD<Type>(a, old_t, b, old_s, old_r);
	}

}
