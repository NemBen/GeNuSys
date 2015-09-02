#include "number_traits.h"

namespace GeNuSys {
	namespace Tests {

		template<typename T>
		GeNuSys::LinAlg::Matrix<T> TestUtils::randomMatrix(unsigned int rows, unsigned int cols, unsigned int range) {
			GeNuSys::LinAlg::Matrix<T> result(rows, cols);
			for (unsigned int i = 0; i < rows; ++i) {
				for (unsigned int j = 0; j < cols; ++j) {
					result.set(i, j, GeNuSys::TypeTraits<int>::template asType<T>(rand() % (2 * range + 1) - range));
				}
			}

			return result;
		}

		template<typename T>
		GeNuSys::LinAlg::Vector<T> TestUtils::randomVector(unsigned int length, unsigned int range) {
			GeNuSys::LinAlg::Vector<T> result(length);
			for (unsigned int i = 0; i < length; ++i) {
				result.set(i, GeNuSys::TypeTraits<int>::template asType<T>(rand() % (2 * range + 1) - range));
			}

			return result;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = (vctA[i] == GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = (vctA[i] == GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = (vctA[i] == GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = (vctA[i] == GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = (matA(i,j) == GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = (matA(i,j) == GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = (matA(i,j) == GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = (matA(i,j) == TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = NumberTraits<T>::isEpsilon(vctA[i] - GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = NumberTraits<T>::isEpsilon(vctA[i] - GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = NumberTraits<T>::isEpsilon(vctA[i] - GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB) {
			if (vctA.getLength() != vctB.getLength()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < vctA.getLength() && equals; ++i) {
				equals = NumberTraits<T>::isEpsilon(vctA[i] - GeNuSys::TypeTraits<S>::template asType<T>(vctB[i]));
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = NumberTraits<T>::isEpsilon(matA(i,j) - GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = NumberTraits<T>::isEpsilon(matA(i,j) - GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = NumberTraits<T>::isEpsilon(matA(i,j) - GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

		template<typename T, typename S>
		bool TestUtils::nearEquals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB) {
			if (matA.getRows() != matB.getRows() || matA.getCols() != matB.getCols()) {
				return false;
			}

			bool equals = true;
			for (unsigned int i = 0; i < matA.getRows() && equals; ++i) {
				for (unsigned int j = 0; j < matA.getCols() && equals; ++j) {
					equals = NumberTraits<T>::isEpsilon(matA(i,j) - GeNuSys::TypeTraits<S>::template asType<T>(matB(i,j)));
				}
			}

			return equals;
		}

	}
}
