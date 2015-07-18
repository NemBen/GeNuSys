#ifndef GENUSYS_LINALG_SPARSE_MATRIX_H_
#define GENUSYS_LINALG_SPARSE_MATRIX_H_

#include <iostream>
#include <vector>

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		class Vector;

		template<typename Type>
		class SparseVector;

		template<typename Type>
		class Matrix;

		template<typename Type>
		class SparseMatrix {

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

			template<unsigned int p>
			friend struct PNorm;

			friend struct FrobeniusNorm;

			template<typename T>
			friend class OperatorNorm;

		public:

			static SparseMatrix<Type> diag(const Vector<Type>& vct);

			static SparseMatrix<Type> diag(const SparseVector<Type>& vct);

			static SparseMatrix<Type> identity(unsigned int rows, unsigned int cols);

		public:

			struct Entry {

				unsigned int col_idx;

				Type value;
	
				Entry();

				Entry(unsigned int col_idx, const Type& value);

			};

		private:

			unsigned int rows;

			unsigned int cols;

			std::vector<unsigned int> row_ptr;

			std::vector<Entry> elem;

			unsigned int size() const;

			unsigned int search(unsigned int row_idx, unsigned int col_idx) const;

			void push(unsigned int col_idx, const Type& value);

			SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity, unsigned int size);

			SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity);

		public:

			SparseMatrix();

			SparseMatrix(unsigned int rows, unsigned int cols);

			SparseMatrix(const SparseMatrix<Type>& mat);

			template<typename SourceType>
			SparseMatrix(const SparseMatrix<SourceType>& mat);

			template<typename SourceType>
			SparseMatrix(const Matrix<SourceType>& mat);

			virtual ~SparseMatrix();

			SparseMatrix<Type>& operator =(const SparseMatrix<Type>& mat);

			template<typename SourceType>
			SparseMatrix<Type>& operator =(const SparseMatrix<SourceType>& mat);

			template<typename SourceType>
			SparseMatrix<Type>& operator =(const Matrix<SourceType>& mat);

			// Component access

			unsigned int getRows() const;

			unsigned int getCols() const;

			void set(unsigned int row_idx, unsigned int col_idx, const Type& value);

			const Type& operator ()(unsigned int row_idx, unsigned int col_idx) const;

			// Matrix

			SparseMatrix<Type> transpose() const;

			SparseMatrix<Type> conjugateTranspose() const;

			// Matrix - Value

			SparseMatrix<Type> operator *(const Type& value) const;

			SparseMatrix<typename TypeTraits<Type>::RationalType> operator /(const Type& value) const;

			// Matrix - Vector

			Vector<Type> operator *(const Vector<Type>& vct) const;

			Vector<Type> operator *(const SparseVector<Type>& vct) const;

			// Matrix - Matrix

			Matrix<Type> operator +(const Matrix<Type>& mat) const;

			SparseMatrix<Type> operator +(const SparseMatrix<Type>& mat) const;

			Matrix<Type> operator -(const Matrix<Type>& mat) const;

			SparseMatrix<Type> operator -(const SparseMatrix<Type>& mat) const;

			Matrix<Type> operator *(const Matrix<Type>& mat) const;

			Matrix<Type> operator *(const SparseMatrix<Type>& mat) const;

			bool operator ==(const Matrix<Type>& mat);

			bool operator ==(const SparseMatrix<Type>& mat);

			bool operator !=(const Matrix<Type>& mat);

			bool operator !=(const SparseMatrix<Type>& mat);

			// IO

			friend std::ostream& operator <<(std::ostream &os, const SparseMatrix<Type>& mat) {
				os << "[ " << std::endl;
				for (unsigned int row = 0; row < mat.rows; ++row) {
					for (unsigned int i = mat.row_ptr[row]; i < mat.row_ptr[row + 1]; ++i) {
						os << " (" << row << ", " << mat.elem[i].col_idx << ") -> " << mat.elem[i].value << std::endl;
					}
				}
				os << "]";
				return os;
			}

		};

	}
}

// Include implementation
#include "sparse_matrix.hpp"

#endif // GENUSYS_LINALG_SPARSE_MATRIX_H_
