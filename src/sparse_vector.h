#ifndef GENUSYS_LINALG_SPARSE_VECTOR_H_
#define GENUSYS_LINALG_SPARSE_VECTOR_H_

#include <iostream>
#include <vector>

#include "number_traits.h"

#include "vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		class Vector;

		template<typename Type>
		class Matrix;

		template<typename Type>
		class SparseMatrix;

		template<typename Type>
		class SparseVector {

			friend struct Traits;

			friend struct Algorithms;

			friend struct Operations;

			template<typename T>
			friend class Vector;

			template<typename T>
			friend class SparseVector;

			template<typename T>
			friend class Matrix;

			template<typename T>
			friend class SparseMatrix;

			template<unsigned int P>
			friend struct PNorm;

			template<typename T>
			friend class OperatorNorm;

		public:

			struct Entry {

				unsigned int idx;

				Type value;
	
				Entry();

				Entry(unsigned int idx, const Type& value);

			};

		private:

			unsigned int length;

			std::vector<Entry> elem;

			unsigned int size() const;

			unsigned int search(unsigned int idx) const;

			void push(unsigned int idx, const Type& value);

			SparseVector(unsigned int length, unsigned int size);

		public:

			SparseVector();
			
			SparseVector(unsigned int length);

			SparseVector(const SparseVector<Type>& vct);

			template<typename SourceType>
			SparseVector(const SparseVector<SourceType>& vct);

			template<typename SourceType>
			SparseVector(const Vector<SourceType>& vct);

			virtual ~SparseVector();

			SparseVector<Type>& operator =(const SparseVector<Type>& vct);

			template<typename SourceType>
			SparseVector<Type>& operator =(const SparseVector<SourceType>& vct);

			template<typename SourceType>
			SparseVector<Type>& operator =(const Vector<SourceType>& vct);

			// Component access

			unsigned int getLength() const;

			void set(unsigned int idx, const Type& value);

			const Type& operator [](unsigned int idx) const;

			// Vector - Value

			SparseVector<Type> operator *(const Type& value) const;

			SparseVector<typename TypeTraits<Type>::RationalType> operator /(const Type& value) const;

			// Vector - Vector

			Type operator *(const Vector<Type>& vct) const;

			Type operator *(const SparseVector<Type>& vct) const;

			Vector<Type> operator +(const Vector<Type>& vct) const;

			SparseVector<Type> operator +(const SparseVector<Type>& vct) const;

			Vector<Type> operator -(const Vector<Type>& vct) const;

			SparseVector<Type> operator -(const SparseVector<Type>& vct) const;

			bool operator ==(const Vector<Type>& vct) const;

			bool operator ==(const SparseVector<Type>& vct) const;

			bool operator !=(const Vector<Type>& vct) const;

			bool operator !=(const SparseVector<Type>& vct) const;

			bool operator <(const Vector<Type>& vct) const;

			bool operator <(const SparseVector<Type>& vct) const;

			bool operator <=(const Vector<Type>& vct) const;

			bool operator <=(const SparseVector<Type>& vct) const;

			bool operator >(const Vector<Type>& vct) const;

			bool operator >(const SparseVector<Type>& vct) const;

			bool operator >=(const Vector<Type>& vct) const;

			bool operator >=(const SparseVector<Type>& vct) const;

			// IO

			friend std::ostream& operator <<(std::ostream &os, const SparseVector<Type>& vct) {
				os << "[";
				for (unsigned int i = 0; i < vct.size(); ++i) {
					os << " (" << vct.elem[i].idx << ", " << vct.elem[i].value << ")";
				}
				os << " ]";
				return os;
			}
			
		};

	}
}

// Include implementation
#include "sparse_vector.hpp"

#endif // GENUSYS_LINALG_SPARSE_VECTOR_H_
