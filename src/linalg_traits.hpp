#include <algorithm>
#include <complex>

namespace GeNuSys {
	namespace LinAlg {

		template<typename SourceType, typename ConvertedType>
		Matrix<ConvertedType> Traits::convertUnsafe(const Matrix<SourceType>& mat) {
			Matrix<ConvertedType> result(mat.rows, mat.cols, 00);
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = TypeTraits<SourceType>::template asTypeUnsafe<ConvertedType>(mat.elem[i]);
			}

			return result;
		}

		template<typename SourceType, typename ConvertedType>
		SparseMatrix<ConvertedType> Traits::convertUnsafe(const SparseMatrix<SourceType>& mat) {
			SparseMatrix<ConvertedType> result(mat.rows, mat.cols);
			result.row_ptr = mat.row_ptr;
			result.elem = std::vector<typename SparseMatrix<ConvertedType>::Entry>(mat.elem.size());
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<ConvertedType>::Entry(mat.elem[i].col_idx, TypeTraits<SourceType>::template asTypeUnsafe<ConvertedType>(mat.elem[i].value));
			}

			return result;
		}

		template<typename SourceType, typename ConvertedType>
		Vector<ConvertedType> Traits::convertUnsafe(const Vector<SourceType>& vct) {
			Vector<ConvertedType> result(vct.length, 00);
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = TypeTraits<SourceType>::template asTypeUnsafe<ConvertedType>(vct.elem[i]);
			}

			return result;
		}

		template<typename SourceType, typename ConvertedType>
		SparseVector<ConvertedType> Traits::convertUnsafe(const SparseVector<SourceType>& vct) {
			SparseVector<ConvertedType> result(vct.length);
			result.elem = std::vector<typename SparseVector<ConvertedType>::Entry>(vct.elem.size());
			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<ConvertedType>::Entry(vct.elem[i].idx, TypeTraits<SourceType>::template asTypeUnsafe<ConvertedType>(vct.elem[i].value));
			}

			return result;
		}
		
		template<typename Type>
		void Traits::swapRows(Matrix<Type>& mat, const int rowA, const int rowB) {
			for (unsigned int i = 0, idxA = rowA * mat.cols, idxB = rowB * mat.cols; i < mat.cols; ++i, ++idxA, ++idxB) {
				std::swap(mat.elem[idxA], mat.elem[idxB]);
			}
		}
		
		template<typename Type>
		void Traits::swapCols(Matrix<Type>& mat, const int colA, const int colB) {
			for (unsigned int i = 0, idxA = colA, idxB = colB; i < mat.rows; ++i, idxA += mat.cols, idxB += mat.cols) {
				std::swap(mat.elem[idxA], mat.elem[idxB]);
			}
		}

		template<typename Type>
		Vector<Type> Traits::getSubVector(const Vector<Type>& vct, unsigned int from, unsigned int to) {
			assert(from < to && to <= vct.length);

			Vector<Type> result(to - from, 00);
			for (unsigned int i = 0, idx = from; idx < to; ++i, ++idx) {
				result.elem[i] = vct.elem[idx];
			}

			return result;
		}

		template<typename Type>
		SparseVector<Type> Traits::getSubVector(const SparseVector<Type>& vct, unsigned int from, unsigned int to) {
			assert(from < to && to <= vct.length);

			SparseVector<Type> result(to - from);
			for (unsigned int i = vct.search(from); vct.elem[i].idx < to; ++i) {
				result.push(i - from, vct.elem[i].value);
			}

			return result;
		}

		template<typename Type>
		Vector<Type> Traits::getRow(const Matrix<Type>& mat, unsigned int row) {
			assert(row < mat.rows);

			Vector<Type> result(mat.cols, 00);
			for (unsigned int i = 0, idx = row * mat.cols; i < mat.cols; ++i, ++idx) {
				result.elem[i] = mat.elem[idx];
			}

			return result;
		}

		template<typename Type>
		SparseVector<Type> Traits::getRow(const SparseMatrix<Type>& mat, unsigned int row) {
			assert(row < mat.rows);

			SparseVector<Type> result(mat.cols);
			for (unsigned int i = mat.row_ptr[row]; i < mat.row_ptr[row + 1]; ++i) {
				result.push(mat.elem[i].col_idx, mat.elem[i].value);
			}

			return result;
		}

		template<typename Type>
		Vector<Type> Traits::getRow(const Matrix<Type>& mat, unsigned int row, unsigned int from, unsigned int to) {
			assert(row < mat.rows);
			assert(from < to && to <= mat.cols);

			Vector<Type> result(to - from, 00);
			for (unsigned int i = from, idxA = row * mat.cols + from, idxR = 0; i < to; ++i, ++idxA, ++idxR) {
				result.elem[idxR] = mat.elem[idxA];
			}

			return result;
		}

		template<typename Type>
		SparseVector<Type> Traits::getRow(const SparseMatrix<Type>& mat, unsigned int row, unsigned int from, unsigned int to) {
			assert(row < mat.rows);
			assert(from < to && to <= mat.cols);

			SparseVector<Type> result(to - from);
			for (unsigned int idx = mat.search(row, from); idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < to; ++idx) {
				result.push(mat.elem[idx].col_idx - from, mat.elem[idx].value);
			}

			return result;
		}

		template<typename Type>
		Matrix<Type> Traits::getRows(const Matrix<Type>& mat, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.rows);

			Matrix<Type> result(to - from, mat.cols, 00);
			for (unsigned int i = from, idxA = from * mat.cols, idxR = 0; i < to; ++i) {
				for (unsigned int j = 0; j < mat.cols; ++j, ++idxA, ++idxR) {
					result.elem[idxR] = mat.elem[idxA];
				}
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> Traits::getRows(const SparseMatrix<Type>& mat, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.rows);

			SparseMatrix<Type> result(to - from, mat.cols);
			for (unsigned int i = from, rowR = 0; i < to; ++i, ++rowR) {
				for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i + 1]; ++j) {
					result.push(mat.elem[j].idx, mat.elem[j].value);
				}
				result.row_ptr[rowR + 1] = result.size();
			}

			return result;
		}
		
		template<typename Type>
		Vector<Type> Traits::getCol(const Matrix<Type>& mat, unsigned int col) {
			assert(col < mat.cols);

			Vector<Type> result(mat.rows, 00);
			for (unsigned int i = 0, idxA = col; i < mat.rows; ++i, idxA += mat.cols) {
				result.elem[i] = mat.elem[idxA];
			}
		}
		
		template<typename Type>
		SparseVector<Type> Traits::getCol(const SparseMatrix<Type>& mat, unsigned int col) {
			assert(col < mat.cols);

			SparseVector<Type> result(mat.rows);
			for (unsigned int row = 0; row < mat.rows; ++row) {
				unsigned int idx = mat.search(row, col);
				if (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx == col) {
					result.push(row, mat.elem[idx].value);
				}
			}

			return result;
		}

		template<typename Type>
		Vector<Type> Traits::getCol(const Matrix<Type>& mat, unsigned int col, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.rows);
			assert(col < mat.cols);
			
			Vector<Type> result(to - from, 00);
			for (unsigned int i = from, idxA = from * mat.cols + col, idxR = 0; i < to; ++i, idxA += mat.cols, ++idxR) {
				result.elem[idxR] = mat.elem[idxA];
			}

			return result;
		}

		template<typename Type>
		SparseVector<Type> Traits::getCol(const SparseMatrix<Type>& mat, unsigned int col, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.rows);
			assert(col < mat.cols);

			SparseVector<Type> result(to - from);
			for (unsigned int row = from, rowR = 0; row < to; ++row, ++rowR) {
				unsigned int idx = mat.search(row, col);
				if (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx == col) {
					result.push(rowR, mat.elem[idx].value);
				}
				result.row_ptr[rowR + 1] = result.size();
			}

			return result;
		}

		template<typename Type>
		Matrix<Type> Traits::getCols(const Matrix<Type>& mat, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.cols);

			Matrix<Type> result(mat.rows, to - from, 00);
			for (unsigned int i = 0, idxA = from, idxR = 0; i < mat.rows; ++i, idxA += mat.cols - to + from) {
				for (unsigned int j = from; j < to; ++j, ++idxA, ++idxR) {
					result.elem[idxR] = mat.elem[idxA];
				}
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> Traits::getCols(const SparseMatrix<Type>& mat, unsigned int from, unsigned int to) {
			assert(from < to && to <= mat.cols);

			SparseMatrix<Type> result(mat.rows, to - from);
			for (unsigned int row = 0; row < mat.rows; ++row) {
				unsigned int idx = mat.search(row, from);
				while (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < to) {
					result.push(mat.elem[idx].col_idx - from, mat.elem[idx].value);
				}
				result.row_ptr[row + 1] = result.size();
			}

			return result;
		}

		template<typename Type>
		Matrix<Type> Traits::getSubMatrix(const Matrix<Type>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol) {
			assert(fromRow < toRow && toRow <= mat.rows);
			assert(fromCol < toCol && toCol <= mat.cols);
			
			Matrix<Type> result(toRow - fromRow, toCol - fromCol, 00);
			for (unsigned int i = fromRow, idxA = fromRow * mat.cols + fromCol, idxR = 0; i < toRow; ++i, idxA += mat.cols - toCol + fromCol) {
				for (unsigned int j = fromCol; j < toCol; ++j, ++idxA, ++idxR) {
					result.elem[idxR] = mat.elem[idxA];
				}
			}

			return result;
		}

		template<typename Type>
		SparseMatrix<Type> Traits::getSubMatrix(const SparseMatrix<Type>& mat, unsigned int fromRow, unsigned int fromCol, unsigned int toRow, unsigned int toCol) {
			assert(fromRow < toRow && toRow <= mat.rows);
			assert(fromCol < toCol && toCol <= mat.cols);

			SparseMatrix<Type> result(toRow - fromRow, toCol - fromCol);
			for (unsigned int row = fromRow, rowR = 0; row < toRow; ++row, ++rowR) {
				unsigned int idx = mat.search(row, fromCol);
				while (idx < mat.row_ptr[row + 1] && mat.elem[idx].col_idx < toCol) {
					result.push(mat.elem[idx].col_idx - fromCol, mat.elem[idx].value);
				}
				result.row_ptr[rowR + 1] = result.size();
			}

			return result;
		}

	}
}
