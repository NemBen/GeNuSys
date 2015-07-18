#ifndef GENUSYS_ELEMENT_TRAITS_H_
#define GENUSYS_ELEMENT_TRAITS_H_

#include "type_traits.h"
#include "integer_traits.h"

namespace GeNuSys {

	template<typename Type>
	struct NumberTraits {

		static const Type zero;

		static const Type one;

		static const Type epsilon;

		static Type sgn(const Type& value);

		static typename TypeTraits<Type>::AbsType abs(const Type& value);

		static typename TypeTraits<Type>::AbsSqrType absSqr(const Type& value);

		static Type conj(const Type& value);
		
		static Type pow(const Type& value, unsigned int n);

		static typename TypeTraits<Type>::RealType sqrt(const Type& value);

		static typename TypeTraits<Type>::RealType root(const Type& value, unsigned int n);

		static typename TypeTraits<Type>::RationalType div(const Type& a, const Type& b);

		static bool isEpsilon(const Type& value);

	};



}

// Include implementation
#include "number_traits.hpp"
#include "number_traits_int.hpp"
#include "number_traits_long_long.hpp"
#include "number_traits_double.hpp"
#include "number_traits_complex.hpp"
#ifdef __unix__
	#include "number_traits_mpz.hpp"
	#include "number_traits_mpq.hpp"
	#include "number_traits_mpf.hpp"
#endif // __unix__

#endif // GENUSYS_ELEMENT_TRAITS_H_
