#include <cmath>
#include <cassert>

namespace GeNuSys {

	template<typename Type>
	Type NumberTraits<Type>::sgn(const Type& value) {
		Type sign = NumberTraits<Type>::zero;
		if (value > NumberTraits<Type>::zero) {
			sign = NumberTraits<Type>::one;
		} else if (value < NumberTraits<Type>::zero) {
			sign = -NumberTraits<Type>::one;
		}

		return sign;
	}

	template<typename Type>
	Type NumberTraits<Type>::conj(const Type& value) {
		return value;
	}

	template<typename Type>
	Type NumberTraits<Type>::pow(const Type& value, unsigned int n) {
		Type result = NumberTraits<Type>::one;
		for (Type acc = value; n > 0; n /= 2, acc *= acc) {
			if (n % 2 == 1) {
				result *= acc;
			}
		}

		return result;
	}

	template<typename Type>
	typename TypeTraits<Type>::RealType NumberTraits<Type>::root(const Type& value, unsigned int n) {
		typedef typename TypeTraits<Type>::RealType RealType;

		RealType r = TypeTraits<Type>::template asType<RealType>(value);
		RealType g = NumberTraits<Type>::one;
		RealType f;

		do {
			f = NumberTraits<RealType>::pow(g, n) - r;
			g -= f / (NumberTraits<RealType>::pow(g, n-1) * n);
		} while (!NumberTraits<RealType>::isEpsilon(f));

		return g;
	}

	template<typename Type>
	bool NumberTraits<Type>::isEpsilon(const Type& value) {
		return abs(value) <= NumberTraits<typename TypeTraits<Type>::AbsType>::epsilon;
	}

}
