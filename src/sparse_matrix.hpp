#include <cassert>
#include <algorithm>

#include "linalg_operations.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		SparseMatrix<Type>::Entry::Entry() {
		}

		template<typename Type>
		SparseMatrix<Type>::Entry::Entry(unsigned int col_idx, const Type& value):col_idx(col_idx),value(value) {
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::diag(const Vector<Type>& vct) {
			SparseMatrix<Type> result(vct.length, vct.length);
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = Entry(i, vct.elem[i]);
				result.row_ptr[i + 1] = i + 1;
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::diag(const SparseVector<Type>& vct) {
			SparseMatrix<Type> result(vct.length, vct.length, vct.size());
			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = Entry(vct.elem[i].idx, vct.elem[i].elem);
				result.row_ptr[i + 1] = i + 1;
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::identity(unsigned int rows, unsigned int cols) {
			unsigned int min = cols < rows ? cols : rows;

			SparseMatrix<Type> result(rows, cols, min);
			for (unsigned int i = 0; i < min; ++i) {
				result.elem[i] = Entry(i, NumberTraits<Type>::one);
				result.row_ptr[i + 1] = i + 1;
			}
			for (unsigned int i = min; i < rows; ++i) {
				result.row_ptr[i + 1] = i;
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type>::SparseMatrix():rows(0),cols(0),row_ptr(0),elem() {
		}

		template<typename Type>
		SparseMatrix<Type>::SparseMatrix(unsigned int rows, unsigned int cols, unsigned int capacity):rows(rows),cols(cols),row_ptr(rows + 1, 0),elem(capacity) {
		}

		template<typename Type>
		SparseMatrix<Type>::SparseMatrix(unsigned int rows, unsigned int cols):rows(rows),cols(cols),row_ptr(rows + 1, 0),elem() {
		}

		template<typename Type>
		SparseMatrix<Type>::SparseMatrix(const SparseMatrix<Type>& mat):rows(mat.rows),cols(mat.cols),row_ptr(mat.row_ptr),elem(mat.elem) {
		}

		template<typename Type>
		template<typename SourceType>
		SparseMatrix<Type>::SparseMatrix(const SparseMatrix<SourceType>& mat):rows(mat.rows),cols(mat.cols),row_ptr(mat.row_ptr),elem(mat.elem.size()) {
			for (unsigned int i = 0; i < size(); ++i) {
				elem[i] = Entry(mat.elem[i].col_idx, TypeTraits<SourceType>::template asType<Type>(mat.elem[i].value));
			}
		}

		template<typename Type>
		template<typename SourceType>
		SparseMatrix<Type>::SparseMatrix(const Matrix<SourceType>& mat):rows(mat.rows),cols(mat.cols),row_ptr(rows + 1, 0) {
			unsigned int nnz = 0;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				if (mat.elem[i] != NumberTraits<SourceType>::zero) {
					++nnz;
				}
			}

			elem = std::vector<Entry>(nnz);
			unsigned int size = 0;
			for (unsigned int i = 0, idxB = 0; i < rows; ++i) {
				for (unsigned int j = 0; j < cols; ++j, ++idxB) {
					if (mat.elem[idxB] != NumberTraits<SourceType>::zero) {
						elem[size++] = Entry(j, TypeTraits<SourceType>::template asType<Type>(mat.elem[idxB]));
					}
				}
				row_ptr[i + 1] = size;
			}
		}

		template<typename Type>
		SparseMatrix<Type>::~SparseMatrix() {
		}

		template<typename Type>
		SparseMatrix<Type>& SparseMatrix<Type>::operator =(const SparseMatrix<Type>& mat) {
			if (this != &mat) {
				rows = mat.rows;
				cols = mat.cols;
				row_ptr = mat.row_ptr;
				elem = mat.elem;
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		SparseMatrix<Type>& SparseMatrix<Type>::operator =(const SparseMatrix<SourceType>& mat) {
			rows = mat.rows;
			cols = mat.cols;
			row_ptr = mat.row_ptr;
			elem = std::vector<Entry>(mat.elem.size());
			for (unsigned int i = 0; i < size(); ++i) {
				elem[i] = Entry(mat.elem[i].col_idx, TypeTraits<SourceType>::template asType<Type>(mat.elem[i].value));
			}

			return *this;
		}

		template<typename Type>
		template<typename SourceType>
		SparseMatrix<Type>& SparseMatrix<Type>::operator =(const Matrix<SourceType>& mat) {
			rows = mat.rows;
			cols = mat.cols;
			row_ptr = std::vector<unsigned int>(rows + 1, 0);

			unsigned int nnz = 0;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				if (mat.elem[i] != NumberTraits<SourceType>::zero) {
					++nnz;
				}
			}

			elem = std::vector<Entry>(nnz);
			unsigned int size = 0;
			for (unsigned int i = 0, idxB = 0; i < rows; ++i) {
				for (unsigned int j = 0; j < cols; ++j, ++idxB) {
					if (mat.elem[idxB] != NumberTraits<SourceType>::zero) {
						elem[size++] = Entry(j, TypeTraits<SourceType>::template asType<Type>(mat.elem[idxB]));
					}
				}
				row_ptr[i + 1] = size;
			}

			return *this;
		}

		template<typename Type>
		unsigned int SparseMatrix<Type>::size() const {
			return elem.size();
		}

		template<typename Type>
		unsigned int SparseMatrix<Type>::search(unsigned int row_idx, unsigned int col_idx) const {
			unsigned int min = row_ptr[row_idx];
			unsigned int max = row_ptr[row_idx + 1];
			unsigned int mid;
			while (min < max) {
				mid = (min + max) / 2;
				if (elem[mid].col_idx < col_idx) {
					min = mid + 1;
				} else {
					max = mid;
				}
			}
			return min;
		}

		template<typename Type>
		void SparseMatrix<Type>::push(unsigned int col_idx, const Type& value) {
			elem.push_back(Entry(col_idx, value));
		}

		template<typename Type>
		unsigned int SparseMatrix<Type>::getRows() const {
			return rows;
		}

		template<typename Type>
		unsigned int SparseMatrix<Type>::getCols() const {
			return cols;
		}

		template<typename Type>
		void SparseMatrix<Type>::set(unsigned int row_idx, unsigned int col_idx, const Type& value) {
			assert(row_idx < rows && col_idx < cols);

			unsigned int entry_idx = search(row_idx, col_idx);
			if (value != NumberTraits<Type>::zero) {
				if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx) {
					elem[entry_idx].value = value;
				} else {
					elem.insert(elem.begin() + entry_idx, Entry(col_idx, value));
					for (unsigned int i = row_idx + 1; i <= rows; ++i) {
						++row_ptr[i];
					}
				}
			} else {
				if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx) {
					elem.erase(elem.begin() + entry_idx);
					for (unsigned int i = row_idx + 1; i <= rows; ++i) {
						++row_ptr[i];
					}
				}
			}
		}

		template<typename Type>
		const Type& SparseMatrix<Type>::operator ()(unsigned int row_idx, unsigned int col_idx) const {
			assert(row_idx < rows && col_idx < cols);

			unsigned int entry_idx = search(row_idx, col_idx);
			if (entry_idx < row_ptr[row_idx + 1] && elem[entry_idx].col_idx == col_idx) {
				return elem[entry_idx].value;
			} else {
				return NumberTraits<Type>::zero;
			}
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::transpose() const {
			SparseMatrix<Type> result(cols, rows, size());
			Operations::mat_transpose(*this, result);

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::conjugateTranspose() const {
			SparseMatrix<Type> result(cols, rows, size());
			Operations::mat_conjugate_transpose(*this, result);

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::operator *(const Type& value) const {
			SparseMatrix<Type> result(rows, cols, size());
			Operations::mat_mul(*this, value, result);

			return result;
		}

		template<typename Type>
		SparseMatrix<typename TypeTraits<Type>::RationalType> SparseMatrix<Type>::operator /(const Type& value) const {
			typedef typename TypeTraits<Type>::RationalType RationalType;

			SparseMatrix<RationalType> result(rows, cols, size());
			Operations::mat_div(*this, value, result);

			return result;
		}

		template<typename Type>
		Vector<Type> SparseMatrix<Type>::operator *(const Vector<Type>& vct) const {
			assert(cols == vct.length);

			Vector<Type> result(rows, 00);
			Operations::mat_mul(*this, vct, result);

			return result;
		}

		template<typename Type>
		Vector<Type> SparseMatrix<Type>::operator *(const SparseVector<Type>& vct) const {
			assert(cols == vct.length);

			Vector<Type> result(rows, 00);
			Operations::mat_mul(*this, vct, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> SparseMatrix<Type>::operator +(const Matrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_add(*this, mat, result);

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::operator +(const SparseMatrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			SparseMatrix<Type> result(rows, cols);
			Operations::mat_add(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> SparseMatrix<Type>::operator -(const Matrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			Matrix<Type> result(rows, cols, 00);
			Operations::mat_sub(*this, mat, result);

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> SparseMatrix<Type>::operator -(const SparseMatrix<Type>& mat) const {
			assert(rows == mat.rows && cols == mat.cols);

			SparseMatrix<Type> result(rows, cols);
			Operations::mat_sub(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> SparseMatrix<Type>::operator *(const Matrix<Type>& mat) const {
			assert(cols == mat.rows);

			Matrix<Type> result(rows, mat.cols, 00);
			Operations::mat_mul(*this, mat, result);

			return result;
		}

		template<typename Type>
		Matrix<Type> SparseMatrix<Type>::operator *(const SparseMatrix<Type>& mat) const {
			assert(cols == mat.rows);

			Matrix<Type> result(rows, mat.cols, 00);
			Operations::mat_mul(*this, mat, result);

			return result;
		}

		template<typename Type>
		bool SparseMatrix<Type>::operator ==(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool SparseMatrix<Type>::operator ==(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool SparseMatrix<Type>::operator !=(const Matrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return !Operations::mat_eq(*this, mat);
		}

		template<typename Type>
		bool SparseMatrix<Type>::operator !=(const SparseMatrix<Type>& mat) {
			assert(rows == mat.rows && cols == mat.cols);

			return !Operations::mat_eq(*this, mat);
		}

	}
}
