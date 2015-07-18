#ifndef GENUSYS_LINALG_MATRIX_H_
#define GENUSYS_LINALG_MATRIX_H_

#include <iostream>
#include <vector>

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		class Vector;

		template<typename Type>
		class SparseVector;

		template<typename Type>
		class SparseMatrix;

		template<typename Type>
		class Matrix {

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

			static Matrix<Type> diag(const Vector<Type>& vct);

			static Matrix<Type> diag(const SparseVector<Type>& vct);

			static Matrix<Type> identity(unsigned int rows, unsigned int cols);

		private:

			unsigned int rows;

			unsigned int cols;

			unsigned int size() const;

			std::vector<Type> elem;

			Matrix(unsigned int rows, unsigned int cols, int);

		public:

			Matrix();

			Matrix(unsigned int rows, unsigned int cols);

			Matrix(const Matrix<Type>& mat);

			template<typename SourceType>
			Matrix(const Matrix<SourceType>& mat);

			template<typename SourceType>
			Matrix(const SparseMatrix<SourceType>& mat);

			virtual ~Matrix();

			Matrix<Type>& operator =(const Matrix<Type>& mat);

			template<typename SourceType>
			Matrix<Type>& operator =(const Matrix<SourceType>& mat);

			template<typename SourceType>
			Matrix<Type>& operator =(const SparseMatrix<SourceType>& mat);

			// Component access

			unsigned int getRows() const;

			unsigned int getCols() const;

			void set(unsigned int row_idx, unsigned int col_idx, const Type& value);

			Type& operator ()(unsigned int row_idx, unsigned int col_idx);

			const Type& operator ()(unsigned int row_idx, unsigned int col_idx) const;

			// Matrix

			Matrix<Type> transpose() const;

			Matrix<Type> conjugateTranspose() const;

			// Matrix - Value

			Matrix<Type> operator *(const Type& value) const;

			Matrix<typename TypeTraits<Type>::RationalType> operator /(const Type& value) const;

			// Matrix - Vector

			Vector<Type> operator *(const Vector<Type>& vct) const;

			Vector<Type> operator *(const SparseVector<Type>& vct) const;

			// Matrix - Matrix

			Matrix<Type> operator +(const Matrix<Type>& mat) const;

			Matrix<Type> operator +(const SparseMatrix<Type>& mat) const;

			Matrix<Type>& operator +=(const Matrix<Type>& mat);

			Matrix<Type>& operator +=(const SparseMatrix<Type>& mat);

			Matrix<Type> operator -(const Matrix<Type>& mat) const;

			Matrix<Type> operator -(const SparseMatrix<Type>& mat) const;

			Matrix<Type>& operator -=(const Matrix<Type>& mat);

			Matrix<Type>& operator -=(const SparseMatrix<Type>& mat);

			Matrix<Type> operator *(const Matrix<Type>& mat) const;

			Matrix<Type> operator *(const SparseMatrix<Type>& mat) const;

			bool operator ==(const Matrix<Type>& mat);

			bool operator ==(const SparseMatrix<Type>& mat);

			bool operator !=(const Matrix<Type>& mat);

			bool operator !=(const SparseMatrix<Type>& mat);

			// IO

			friend std::ostream& operator <<(std::ostream &os, const Matrix<Type>& mat) {
				os << "[" << std::endl;
				for (unsigned int i = 0; i < mat.rows; ++i) {
					os << "  [ ";
					for (unsigned int j = 0; j < mat.cols; ++j) {
						os << " " << mat(i,j);
					}
					os << " ]" << std::endl;
				}
				os << "]";
				return os;
			}

		};

	}
}

// Include implementation
#include "matrix.hpp"

#endif // GENUSYS_LINALG_MATRIX_H_
