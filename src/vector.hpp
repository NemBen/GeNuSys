#include <cassert>

#include "linalg_operations.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		Vector<Type>::Vector():length(0),elem(length) {
		}

		template<typename Type>
		Vector<Type>::Vector(unsigned int length, int):length(length),elem(length) {
		}

		template<typename Type>
		Vector<Type>::Vector(unsigned int length):length(length),elem(length, NumberTraits<Type>::zero) {
		}

		template<typename Type>
		Vector<Type>::Vector(const Vector<Type>& vct):length(vct.length),elem(vct.elem) {
		}

		template<typename Type>
		template<typename SourceType>
		Vector<Type>::Vector(const Vector<SourceType>& vct):length(vct.length),elem(vct.length) {
			for (unsigned int i = 0; i < length; ++i) {
				elem[i] = TypeTraits<SourceType>::template asType<Type>(vct.elem[i]);
			}
		}
		
		template<typename Type>
		template<typename SourceType>
		Vector<Type>::Vector(const SparseVector<SourceType>& vct):length(vct.length),elem(vct.length, NumberTraits<SourceType>::zero) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				elem[vct.elem[i].idx] = TypeTraits<SourceType>::template asType<Type>(vct.elem[i].value);
			}
		}

		template<typename Type>
		Vector<Type>::~Vector() {
		}

		template<typename Type>
		Vector<Type>& Vector<Type>::operator =(const Vector<Type>& vct) {
			if (this != &vct) {
				length = vct.length;
				elem = vct.elem;
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		Vector<Type>& Vector<Type>::operator =(const Vector<SourceType>& vct) {
			length = vct.length;
			elem = std::vector<Type>(vct.length);
			for (unsigned int i = 0; i < length; ++i) {
				elem[i] = TypeTraits<SourceType>::template asType<Type>(vct.elem[i]);
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		Vector<Type>& Vector<Type>::operator =(const SparseVector<SourceType>& vct) {
			length = vct.length;
			elem = std::vector<Type>(vct.length);
			for (unsigned int i = 0; i < vct.size(); ++i) {
				elem[vct.elem[i].idx] = TypeTraits<SourceType>::template asType<Type>(vct.elem[i].value);
			}

			return *this;
		}

		template<typename Type>
		unsigned int Vector<Type>::getLength() const {
			return length;
		}

		template<typename Type>
		void Vector<Type>::set(unsigned int idx, const Type& value) {
			assert(idx < length);

			elem[idx] = value;
		}

		template<typename Type>
		Type& Vector<Type>::operator [](unsigned int idx) {
			assert(idx < length);

			return elem[idx];
		}

		template<typename Type>
		const Type& Vector<Type>::operator [](unsigned int idx) const {
			assert(idx < length);

			return elem[idx];
		}

		template<typename Type>
		Vector<Type> Vector<Type>::operator *(const Type& value) const {
			Vector<Type> result(length, 00);
			Operations::vct_mul(*this, value, result);

			return result;
		}

		template<typename Type>
		Vector<typename TypeTraits<Type>::RationalType> Vector<Type>::operator /(const Type& value) const {
			Vector<typename TypeTraits<Type>::RationalType> result(length, 00);
			Operations::vct_div(*this, value, result);

			return result;
		}

		template<typename Type>
		Type Vector<Type>::operator *(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_mul(*this, vct);
		}
		
		template<typename Type>
		Type Vector<Type>::operator *(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_mul(*this, vct);
		}

		template<typename Type>
		Vector<Type> Vector<Type>::operator +(const Vector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_add(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type> Vector<Type>::operator +(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_add(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type>& Vector<Type>::operator +=(const Vector<Type>& vct) {
			assert(length == vct.length);

			Operations::vct_add(*this, vct);

			return *this;
		}

		template<typename Type>
		Vector<Type>& Vector<Type>::operator +=(const SparseVector<Type>& vct) {
			assert(length == vct.length);

			Operations::vct_add(*this, vct);

			return *this;
		}

		template<typename Type>
		Vector<Type> Vector<Type>::operator -(const Vector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_sub(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type> Vector<Type>::operator -(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_sub(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type>& Vector<Type>::operator -=(const Vector<Type>& vct) {
			assert(length == vct.length);

			Operations::vct_sub(*this, vct);

			return *this;
		}

		template<typename Type>
		Vector<Type>& Vector<Type>::operator -=(const SparseVector<Type>& vct) {
			assert(length == vct.length);

			Operations::vct_sub(*this, vct);

			return *this;
		}

		template<typename Type>
		bool Vector<Type>::operator ==(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator ==(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator !=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return !Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator !=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return !Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator <(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_lt(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator <(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_lt(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator <=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_leq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator <=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_leq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator >(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_gt(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator >(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_gt(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator >=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_geq(*this, vct);
		}

		template<typename Type>
		bool Vector<Type>::operator >=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_geq(*this, vct);
		}

	}
}
