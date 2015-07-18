#include <complex>

namespace GeNuSys {
	
	template<typename Type>
	struct NumberTraits<std::complex<Type> > {

		static const std::complex<Type> zero;

		static const std::complex<Type> one;

		static std::complex<typename TypeTraits<Type>::RealType> sgn(const std::complex<Type>& value);

		static typename TypeTraits<std::complex<Type> >::AbsType abs(const std::complex<Type>& value);

		static typename TypeTraits<std::complex<Type> >::AbsSqrType absSqr(const std::complex<Type>& value);

		static std::complex<Type> conj(const std::complex<Type>& value);
		
		static std::complex<Type> pow(const std::complex<Type>& value, unsigned int n);

		static typename TypeTraits<std::complex<Type> >::RealType sqrt(const std::complex<Type>& value);

		static typename TypeTraits<std::complex<Type> >::RationalType div(const std::complex<Type>& a, const std::complex<Type>& b);

	};

	template<typename Type>
	const std::complex<Type> NumberTraits<std::complex<Type> >::zero = std::complex<Type>(NumberTraits<Type>::zero, NumberTraits<Type>::zero);

	template<typename Type>
	const std::complex<Type> NumberTraits<std::complex<Type> >::one = std::complex<Type>(NumberTraits<Type>::one, NumberTraits<Type>::zero);

	template<typename Type>
	std::complex<typename TypeTraits<Type>::RealType> NumberTraits<std::complex<Type> >::sgn(const std::complex<Type>& value) {
		typedef typename TypeTraits<Type>::RealType RealType;

		RealType l = abs(value);
		if (l <= NumberTraits<RealType>::epsilon) {
			return zero;
		}

		return std::complex<RealType>(NumberTraits<Type>::div(value.real(), l), NumberTraits<Type>::div(value.imag(), l));
	}

	template<typename Type>
	typename TypeTraits<std::complex<Type> >::AbsType NumberTraits<std::complex<Type> >::abs(const std::complex<Type>& value) {
		return NumberTraits<Type>::sqrt(value.real() * value.real() + value.imag() * value.imag());
	}

	template<typename Type>
	typename TypeTraits<std::complex<Type> >::AbsSqrType NumberTraits<std::complex<Type> >::absSqr(const std::complex<Type>& value) {
		return value.real() * value.real() + value.imag() * value.imag();
	}

	template<typename Type>
	std::complex<Type> NumberTraits<std::complex<Type> >::conj(const std::complex<Type>& value) {
		return std::conj(value);
	}

	template<typename Type>
	typename TypeTraits<std::complex<Type> >::RealType NumberTraits<std::complex<Type> >::sqrt(const std::complex<Type>& value) {
		typedef typename TypeTraits<Type>::RealType RealType;

		if (NumberTraits<Type>::isEpsilon(value.imag())) {
			if (value.real() < 0) {
				return std::complex<RealType>(NumberTraits<RealType>::zero, NumberTraits<Type>::sqrt(-value.real()));
			} else {
				return std::complex<RealType>(NumberTraits<Type>::sqrt(value.real()), NumberTraits<RealType>::zero);
			}
		}

		RealType l = abs(value);
		RealType re = NumberTraits<RealType>::sqrt((value.real() + l) / 2);
		RealType im = NumberTraits<RealType>::sgn(value.imag()) * NumberTraits<RealType>::sqrt((-value.real() + l) / 2);

		return std::complex<RealType>(re, im);
	}

	template<typename Type>
	typename TypeTraits<std::complex<Type> >::RationalType NumberTraits<std::complex<Type> >::div(const std::complex<Type>& a, const std::complex<Type>& b) {
		std::complex<Type> mul = a * std::conj(b);
		Type lenSqr = absSqr(b);

		return  std::complex<typename TypeTraits<Type>::RationalType>(NumberTraits<Type>::div(mul.real(), lenSqr), NumberTraits<Type>::div(mul.imag(), lenSqr));
	}

	template <class Type>
	bool complex_comparator(const std::complex<Type>& a, const std::complex<Type>& b) {
		return real(a) == real(b) ? imag(a) < imag(b) : real(a) < real(b);
	}

}