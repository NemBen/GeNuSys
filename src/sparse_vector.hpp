#include <cassert>

#include "linalg_operations.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		SparseVector<Type>::Entry::Entry() {
		}

		template<typename Type>
		SparseVector<Type>::Entry::Entry(unsigned int idx, const Type& value):idx(idx),value(value) {
		}

		template<typename Type>
		SparseVector<Type>::SparseVector():length(0),elem() {
		}

		template<typename Type>
		SparseVector<Type>::SparseVector(unsigned int length, unsigned int size):length(length),elem(size) {
		}

		template<typename Type>
		SparseVector<Type>::SparseVector(unsigned int length):length(length),elem() {
		}

		template<typename Type>
		SparseVector<Type>::SparseVector(const SparseVector<Type>& vct):length(vct.length),elem(vct.elem) {
		}

		template<typename Type>
		template<typename SourceType>
		SparseVector<Type>::SparseVector(const SparseVector<SourceType>& vct):length(vct.length),elem(vct.elem.size()) {
			for (unsigned int i = 0; i < elem.size(); ++i) {
				elem[i] = Entry(i, TypeTraits<SourceType>::template asType<Type>(vct.elem[i].value));
			}
		}

		template<typename Type>
		template<typename SourceType>
		SparseVector<Type>::SparseVector(const Vector<SourceType>& vct):length(vct.length) {
			unsigned int nnz = 0;
			for (unsigned int i = 0; i < length; ++i) {
				if (vct.elem[i] != NumberTraits<SourceType>::zero) {
					++nnz;
				}
			}

			elem = std::vector<Entry>(nnz);
			unsigned int size = 0;
			for (unsigned int i = 0; i < length; ++i) {
				if (vct.elem[i] != NumberTraits<SourceType>::zero) {
					elem[size++] = Entry(i, TypeTraits<SourceType>::template asType<Type>(vct.elem[i]));
				}
			}
		}

		template<typename Type>
		SparseVector<Type>::~SparseVector() {
		}

		template<typename Type>
		SparseVector<Type>& SparseVector<Type>::operator =(const SparseVector<Type>& vct) {
			if (this != &vct) {
				length = vct.length;
				elem = vct.elem;
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		SparseVector<Type>& SparseVector<Type>::operator =(const SparseVector<SourceType>& vct) {
			length = vct.length;
			elem = std::vector<Entry>(vct.elem.size());
			for (unsigned int i = 0; i < elem.size(); ++i) {
				elem[i] = Entry(i, TypeTraits<SourceType>::template asType<Type>(vct.elem[i].value));
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		SparseVector<Type>& SparseVector<Type>::operator =(const Vector<SourceType>& vct) {
			length = vct.length;
			
			unsigned int nnz = 0;
			for (unsigned int i = 0; i < length; ++i) {
				if (vct.elem[i] != NumberTraits<SourceType>::zero) {
					++nnz;
				}
			}

			elem = std::vector<Entry>(nnz);
			unsigned int size = 0;
			for (unsigned int i = 0; i < length; ++i) {
				if (vct.elem[i] != NumberTraits<SourceType>::zero) {
					elem[size++] = Entry(i, TypeTraits<SourceType>::template asType<Type>(vct.elem[i]));
				}
			}

			return *this;
		}

		template<typename Type>
		unsigned int SparseVector<Type>::size() const {
			return elem.size();
		}

		template<typename Type>
		unsigned int SparseVector<Type>::search(unsigned int idx) const {
			unsigned int min = 0;
			unsigned int max = size();
			unsigned int mid;
			while (min < max) {
				mid = (min + max) / 2;
				if (elem[mid].idx < idx) {
					min = mid + 1;
				} else {
					max = mid;
				}
			}
			return min;
		}

		template<typename Type>
		void SparseVector<Type>::push(unsigned int idx, const Type& value) {
			elem.push_back(Entry(idx, value));
		}

		template<typename Type>
		unsigned int SparseVector<Type>::getLength() const {
			return length;
		}

		template<typename Type>
		void SparseVector<Type>::set(unsigned int idx, const Type& value) {
			assert(idx < length);

			unsigned int entry_idx = search(idx);
			if (value != NumberTraits<Type>::zero) {
				if (entry_idx < size() && elem[entry_idx].idx == idx) {
					elem[entry_idx].value = value;
				} else {
					elem.insert(elem.begin() + entry_idx, Entry(idx, value));
				}
			} else {
				if (entry_idx < size() && elem[entry_idx].idx == idx) {
					elem.erase(elem.begin() + entry_idx);
				}
			}
		}

		template<typename Type>
		const Type& SparseVector<Type>::operator [](unsigned int idx) const {
			assert(idx < length);

			unsigned int entry_idx = search(idx);
			if (entry_idx < size() && elem[entry_idx].idx == idx) {
				return elem[entry_idx].value;
			} else {
				return NumberTraits<Type>::zero;
			}
		}

		template<typename Type>
		SparseVector<Type> SparseVector<Type>::operator *(const Type& value) const {
			SparseVector<Type> result(length, size());
			Operations::vct_mul(*this, value, result);

			return result;
		}

		template<typename Type>
		SparseVector<typename TypeTraits<Type>::RationalType> SparseVector<Type>::operator /(const Type& value) const {
			SparseVector<typename TypeTraits<Type>::RationalType> result(length, size());
			Operations::vct_div(*this, value, result);

			return result;
		}

		template<typename Type>
		Type SparseVector<Type>::operator *(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_mul(*this, vct);
		}

		template<typename Type>
		Type SparseVector<Type>::operator *(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_mul(*this, vct);
		}

		template<typename Type>
		Vector<Type> SparseVector<Type>::operator +(const Vector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_add(*this, vct, result);

			return result;
		}

		template<typename Type>
		SparseVector<Type> SparseVector<Type>::operator +(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			SparseVector<Type> result(length);
			Operations::vct_add(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type> SparseVector<Type>::operator -(const Vector<Type>& vct) const {
			assert(length == vct.length);

			Vector<Type> result(length, 00);
			Operations::vct_sub(*this, vct, result);

			return result;
		}

		template<typename Type>
		SparseVector<Type> SparseVector<Type>::operator -(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			SparseVector<Type> result(length);
			Operations::vct_sub(*this, vct, result);

			return result;
		}

		template<typename Type>
		bool SparseVector<Type>::operator ==(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator ==(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator !=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return !Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator !=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return !Operations::vct_eq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator <(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_lt(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator <(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_lt(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator <=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_leq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator <=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_leq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator >(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_gt(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator >(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_gt(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator >=(const Vector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_geq(*this, vct);
		}

		template<typename Type>
		bool SparseVector<Type>::operator >=(const SparseVector<Type>& vct) const {
			assert(length == vct.length);

			return Operations::vct_geq(*this, vct);
		}

	}
}
