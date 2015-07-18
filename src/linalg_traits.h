#ifndef GENUSYS_LINALG_TRAITS_H_
#define GENUSYS_LINALG_TRAITS_H_

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		struct Traits {

			template<typename SourceType, typename ConvertedType>
			static Matrix<ConvertedType> convertUnsafe(const Matrix<SourceType>& mat);

			template<typename SourceType, typename ConvertedType>
			static SparseMatrix<ConvertedType> convertUnsafe(const SparseMatrix<SourceType>& mat);

			template<typename SourceType, typename ConvertedType>
			static Vector<ConvertedType> convertUnsafe(const Vector<SourceType>& vct);

			template<typename SourceType, typename ConvertedType>
			static SparseVector<ConvertedType> convertUnsafe(const SparseVector<SourceType>& vct);

			template<typename Type>
			static void swapRows(Matrix<Type>& mat, const int rowA, const int rowB);

			template<typename Type>
			static void swapCols(Matrix<Type>& mat, const int colA, const int colB);

			template<typename Type>
			static Vector<Type> getSubVector(const Vector<Type>& vct, unsigned int from, unsigned int to);
			
			template<typename Type>
			static SparseVector<Type> getSubVector(const SparseVector<Type>& vct, unsigned int from, unsigned int to);

			template<typename Type>
			static Vector<Type> getRow(const Matrix<Type>& mat, unsigned int row);

			template<typename Type>
			static SparseVector<Type> getRow(const SparseMatrix<Type>& mat, unsigned int row);

			template<typename Type>
			static Vector<Type> getRow(const Matrix<Type>& mat, unsigned int row, unsigned int from, unsigned int to);

			template<typename Type>
			static SparseVector<Type> getRow(const SparseMatrix<Type>& mat, unsigned int row, unsigned int from, unsigned int to);

			template<typename Type>
			static Matrix<Type> getRows(const Matrix<Type>& mat, unsigned int from, unsigned int to);

			template<typename Type>
			static SparseMatrix<Type> getRows(const SparseMatrix<Type>& mat, unsigned int from, unsigned int to);

			template<typename Type>
			static Vector<Type> getCol(const Matrix<Type>& mat, unsigned int col);

			template<typename Type>
			static SparseVector<Type> getCol(const SparseMatrix<Type>& mat, unsigned int col);

			template<typename Type>
			static Vector<Type> getCol(const Matrix<Type>& mat, unsigned int col, unsigned int from, unsigned int to);

			template<typename Type>
			static SparseVector<Type> getCol(const SparseMatrix<Type>& mat, unsigned int col, unsigned int from, unsigned int to);

			template<typename Type>
			static Matrix<Type> getCols(const Matrix<Type>& mat, unsigned int from, unsigned int to);

			template<typename Type>
			static SparseMatrix<Type> getCols(const SparseMatrix<Type>& mat, unsigned int from, unsigned int to);

			template<typename Type>
			static Matrix<Type> getSubMatrix(const Matrix<Type>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol);
			
			template<typename Type>
			static SparseMatrix<Type> getSubMatrix(const SparseMatrix<Type>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol);

		};

	}
}

// Include implementation
#include "linalg_traits.hpp"

#endif // GENUSYS_LINALG_TRAITS_H_
