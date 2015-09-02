#ifndef GENUSYS_TESTS_TEST_UTILS_H_
#define GENUSYS_TESTS_TEST_UTILS_H_

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace Tests {

		struct TestUtils {

			template<typename T>
			static GeNuSys::LinAlg::Matrix<T> randomMatrix(unsigned int rows, unsigned int cols, unsigned int range);

			template<typename T>
			static GeNuSys::LinAlg::Vector<T> randomVector(unsigned int length, unsigned int range);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

			template<typename T, typename S>
			static bool equals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::Vector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::Vector<S>& vctB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::SparseVector<T>& vctA, const GeNuSys::LinAlg::SparseVector<S>& vctB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::Matrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::Matrix<S>& matB);

			template<typename T, typename S>
			static bool nearEquals(const GeNuSys::LinAlg::SparseMatrix<T>& matA, const GeNuSys::LinAlg::SparseMatrix<S>& matB);

		};

	}
}

// Include implementation
#include "test_utils.hpp"

#endif // GENUSYS_TESTS_TEST_UTILS_H_
