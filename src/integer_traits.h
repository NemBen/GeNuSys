#ifndef GENUSYS_INTEGER_TRAITS_H_
#define GENUSYS_INTEGER_TRAITS_H_

namespace GeNuSys {

	template<typename Type>
	struct ExtendedGCD {

		const Type a, coefA;

		const Type b, coefB;

		const Type gcd;

		ExtendedGCD(const Type& a, const Type& coefA, const Type& b, const Type& coefB, const Type& gcd);

	};

	template<typename Type>
	struct IntegerTraits {
		
		static bool divisible(const Type& a, const Type& b);

		static Type div(const Type& a, const Type& b);

		static Type mod(const Type& a, const Type& b);

		static Type mods(const Type& a, const Type& b);

		static Type gcd(Type a, Type b);

		static ExtendedGCD<Type> egcd(const Type& a, const Type& b);

	};

}

// Include implementation
#include "integer_traits.hpp"
#include "integer_traits_int.hpp"
#include "integer_traits_long_long.hpp"
#ifdef __unix__
	#include "integer_traits_mpz.hpp"
#endif // __unix__

#endif // GENUSYS_INTEGER_TRAITS_H_
