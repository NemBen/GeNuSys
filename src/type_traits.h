#ifndef GENUSYS_TYPE_TRAITS_H_
#define GENUSYS_TYPE_TRAITS_H_

namespace GeNuSys {

	template<typename Type>
	struct TypeTraitsBase {

	};

	template<typename Type>
	struct TypeTraits : TypeTraitsBase<Type> {

		template<typename ConvertedType>
		static ConvertedType asType(const Type& value);

		template<typename ConvertedType>
		static ConvertedType asTypeUnsafe(const Type& value) {
			return TypeTraits<Type>::template asType<ConvertedType>(value);
		}

	};

}

// Include implementation
#include "type_traits_int.hpp"
#include "type_traits_long_long.hpp"
#include "type_traits_double.hpp"
#include "type_traits_complex.hpp"
#ifdef __unix__
	#include "type_traits_mpz.hpp"
	#include "type_traits_mpq.hpp"
	#include "type_traits_mpf.hpp"
#endif // __unix__


#endif // GENUSYS_TYPE_TRAITS_H_
