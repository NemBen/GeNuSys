#include <cassert>
#include <algorithm>

#include "linalg_operations.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		Matrix<Type> Matrix<Type>::diag(const Vector<Type>& vct) {
			Matrix<Type> result(vct.length, vct.length);
			for (unsigned int i = 0, idx = 0; i < vct.length; ++i, idx += vct.length + 1) {
				result.elem[idx] = vct.elem[i];
			}

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::diag(const SparseVector<Type>& vct) {
			Matrix<Type> result(vct.length, vct.length);
			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[vct.elem[i].idx * (vct.length + 1)] = vct.elem[i].value;
			}

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::identity(unsigned int rows, unsigned int cols) {
			Matrix<Type> result(rows, cols);
			unsigned int min = cols < rows ? cols : rows;
			for (unsigned int i = 0, idx = 0; i < min; ++i, idx += cols + 1) {
				result.elem[idx] = NumberTraits<Type>::one;
			}

			return result;
		}

		template<typename Type>
		Matrix<Type>::Matrix():rows(0),cols(0),elem(0) {
		}

		template<typename Type>
		Matrix<Type>::Matrix(unsigned int rows, unsigned int cols):rows(rows),cols(cols),elem(rows * cols, NumberTraits<Type>::zero) {
		}

		template<typename Type>
		Matrix<Type>::Matrix(unsigned int rows, unsigned int cols, int):rows(rows),cols(cols),elem(rows * cols) {
		}

		template<typename Type>
		Matrix<Type>::Matrix(const Matrix<Type>& mat):rows(mat.rows),cols(mat.cols),elem(mat.elem) {
		}

		template<typename Type>
		template<typename SourceType>
		Matrix<Type>::Matrix(const Matrix<SourceType>& mat):rows(mat.rows),cols(mat.cols),elem(mat.rows * mat.cols) {
			for (unsigned int i = 0; i < mat.elem.size(); ++i) {
				elem[i] = TypeTraits<SourceType>::template asType<Type>(mat.elem[i]);
			}
		}

		template<typename Type>
		template<typename SourceType>
		Matrix<Type>::Matrix(const SparseMatrix<SourceType>& mat):rows(mat.rows),cols(mat.cols),elem(mat.rows * mat.cols, NumberTraits<SourceType>::zero) {
			for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i, idxA += cols) {
				for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i+1]; ++j) {
					elem[idxA + mat.elem[j].col_idx] = TypeTraits<SourceType>::template asType<Type>(mat.elem[j].value);
				}
			}
		}

		template<typename Type>
		Matrix<Type>::~Matrix() {
		}

		template<typename Type>
		Matrix<Type>& Matrix<Type>::operator =(const Matrix<Type>& mat) {
			if (this != &mat) {
				rows = mat.rows;
				cols = mat.cols;
				elem = mat.elem;
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		Matrix<Type>& Matrix<Type>::operator =(const Matrix<SourceType>& mat) {
			rows = mat.rows;
			cols = mat.cols;
			elem = std::vector<Type>(mat.rows * mat.cols);
			for (unsigned int i = 0; i < mat.elem.size(); ++i) {
				elem[i] = TypeTraits<SourceType>::template asType<Type>(mat.elem[i]);
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		Matrix<Type>& Matrix<Type>::operator =(const SparseMatrix<SourceType>& mat) {
			rows = mat.rows;
			cols = mat.cols;
			elem = std::vector<Type>(mat.rows * mat.cols, NumberTraits<SourceType>::zero);
			for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i, idxA += cols) {
				for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i+1]; ++j) {
					elem[idxA + mat.elem[j].col_idx] = TypeTraits<SourceType>::template asType<Type>(mat.elem[j].value);
				}
			}

			return *this;
		}

		template<typename Type>
		unsigned int Matrix<Type>::size() const {
			return elem.size();
		}

		template<typename Type>
		unsigned int Matrix<Type>::getRows() const {
			return rows;
		}

		template<typename Type>
		unsigned int Matrix<Type>::getCols() const {
			return cols;
		}

		template<typename Type>
		void Matrix<Type>::set(unsigned int row_idx, unsigned int col_idx, const Type& value) {
			assert(row_idx < rows && col_idx < cols);

			elem[row_idx * cols + col_idx] = value;
		}

		template<typename Type>
		Type& Matrix<Type>::operator ()(unsigned int row_idx, unsigned int col_idx) {
			assert(row_idx < rows && col_idx < cols);

			return elem[row_idx * cols + col_idx];
		}

		template<typename Type>
		const Type& Matrix<Type>::operator ()(unsigned int row_idx, unsigned int col_idx) const {
			assert(row_idx < rows && col_idx < cols);

			return elem[row_idx * cols + col_idx];
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::transpose() const {
			Matrix<Type> result(cols, rows, 00);
			Operations::mat_transpose(*this, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::conjugateTranspose() const {
			Matrix<Type> result(cols, rows, 00);
			Operations::mat_conjugate_transpose(*this, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator *(const Type& value) const {
			Matrix<Type> result(rows, cols, 00);
			Operations::mat_mul(*this, value, result);

			return result;
		}

		template<typename Type>
		Matrix<typename TypeTraits<Type>::RationalType> Matrix<Type>::operator /(const Type& value) const {
			Matrix<typename TypeTraits<Type>::RationalType> result(rows, cols, 00);
			Operations::mat_div(*this, value, result);

			return result;
		}

		template<typename Type>
		Vector<Type> Matrix<Type>::operator *(const Vector<Type>& vct) const {
			assert(cols == vct.length);

			Vector<Type> result(rows, 00);
			Operations::mat_mul(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type> Matrix<Type>::operator *(const SparseVector<Type>& vct) const {
			assert(cols == vct.length);

			Vector<Type> result(rows, 00);
			Operations::mat_mul(*this, vct, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator +(const Matrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_add(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator +(const SparseMatrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_add(*this, mat, result);
			
			return result;
		}

		template<typename Type>
		Matrix<Type>& Matrix<Type>::operator +=(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			Operations::mat_add(*this, mat);

			return *this;
		}

		template<typename Type>
		Matrix<Type>& Matrix<Type>::operator +=(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			Operations::mat_add(*this, mat);

			return *this;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator -(const Matrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_sub(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator -(const SparseMatrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_sub(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type>& Matrix<Type>::operator -=(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			Operations::mat_sub(*this, mat);

			return *this;
		}

		template<typename Type>
		Matrix<Type>& Matrix<Type>::operator -=(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			Operations::mat_sub(*this, mat);

			return *this;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator *(const Matrix<Type>& mat) const {
			assert(cols == mat.rows);

			Matrix<Type> result(rows, mat.cols, 00);
			Operations::mat_mul(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> Matrix<Type>::operator *(const SparseMatrix<Type>& mat) const {
			assert(cols == mat.rows);

			Matrix<Type> result(rows, mat.cols, 00);
			Operations::mat_mul(*this, mat, result);

			return result;
		}

		template<typename Type>
		bool Matrix<Type>::operator ==(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool Matrix<Type>::operator ==(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool Matrix<Type>::operator !=(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return !Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool Matrix<Type>::operator !=(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return !Operations::mat_eq(*this, mat);
		}

	}
}
