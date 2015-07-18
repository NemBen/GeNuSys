#ifndef GENUSYS_LINALG_VECTOR_H_
#define GENUSYS_LINALG_VECTOR_H_

#include <iostream>
#include <vector>

#include "number_traits.h"

#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		class SparseVector;

		template<typename Type>
		class Matrix;

		template<typename Type>
		class SparseMatrix;

		template<typename Type>
		class Vector {

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

		private:

			unsigned int length;

			std::vector<Type> elem;

			Vector(unsigned int length, int);

		public:

			Vector();

			Vector(unsigned int length);

			Vector(const Vector<Type>& vct);

			template<typename SourceType>
			Vector(const Vector<SourceType>& vct);

			template<typename SourceType>
			Vector(const SparseVector<SourceType>& vct);

			virtual ~Vector();

			Vector<Type>& operator =(const Vector<Type>& vct);

			template<typename SourceType>
			Vector<Type>& operator =(const Vector<SourceType>& vct);

			template<typename SourceType>
			Vector<Type>& operator =(const SparseVector<SourceType>& vct);

			// Component access

			unsigned int getLength() const;

			void set(unsigned int idx, const Type& value);

			Type& operator [](unsigned int idx);

			const Type& operator [](unsigned int idx) const;

			// Vector - Value

			Vector<Type> operator *(const Type& value) const;

			Vector<typename TypeTraits<Type>::RationalType> operator /(const Type& value) const;

			// Vector - Vector

			Type operator *(const Vector<Type>& vct) const;

			Type operator *(const SparseVector<Type>& vct) const;

			Vector<Type> operator +(const Vector<Type>& vct) const;

			Vector<Type> operator +(const SparseVector<Type>& vct) const;

			Vector<Type>& operator +=(const Vector<Type>& vct);

			Vector<Type>& operator +=(const SparseVector<Type>& vct);

			Vector<Type> operator -(const Vector<Type>& vct) const;

			Vector<Type> operator -(const SparseVector<Type>& vct) const;
			
			Vector<Type>& operator -=(const Vector<Type>& vct);

			Vector<Type>& operator -=(const SparseVector<Type>& vct);

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

			friend std::ostream& operator <<(std::ostream &os, const Vector<Type>& vct) {
				os << "[ ";
				for (unsigned int i = 0; i < vct.length; ++i) {
					os << vct.elem[i] << " ";
				}
				os << "]";
				return os;
			}

		};

	}
}

// Include implementation
#include "vector.hpp"

#endif // GENUSYS_LINALG_VECTOR_H_
