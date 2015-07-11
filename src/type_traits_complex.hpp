#include <complex>

namespace GeNuSys {

	template<typename Type>
	struct TypeTraitsBase<std::complex<Type> > {

		typedef std::complex<typename TypeTraitsBase<Type>::RealType> RealType;

		typedef std::complex<typename TypeTraitsBase<Type>::RationalType> RationalType;

		typedef std::complex<Type> ComplexType;

		typedef typename TypeTraitsBase<Type>::RealType AbsType;

		typedef Type AbsSqrType;

	};

	template<typename Type>
	struct TypeTraits<std::complex<Type> > : TypeTraitsBase<std::complex<Type> > {

		template<typename ConvertedType>
		static std::complex<ConvertedType> asType(const std::complex<Type>& value) {
			return std::complex<ConvertedType>(
				TypeTraits<Type>::template asType<ConvertedType>(value.real()),
				TypeTraits<Type>::template asType<ConvertedType>(value.imag()));
		}

		template<typename ConvertedType>
		static std::complex<ConvertedType> asTypeUnsafe(const std::complex<Type>& value) {
			return std::complex<ConvertedType>(
				TypeTraits<Type>::template asTypeUnsafe<ConvertedType>(value.real()),
				TypeTraits<Type>::template asTypeUnsafe<ConvertedType>(value.imag()));
		}

	};

}
