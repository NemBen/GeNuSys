#ifndef GENUSYS_LINALG_FROBENIUS_NORM_H_
#define GENUSYS_LINALG_FROBENIUS_NORM_H_

#include "number_traits.h"

#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		struct FrobeniusNorm {

			template<typename Type>
			struct NormType {

				typedef typename TypeTraits<typename TypeTraits<Type>::AbsSqrType>::RealType Value;

			};

			template<typename Type>
			static typename NormType<Type>::Value norm(const Matrix<Type>& mat);

			template<typename Type>
			static typename NormType<Type>::Value norm(const SparseMatrix<Type>& mat);

			template<typename Type>
			typename NormType<Type>::Value operator ()(const Matrix<Type>& mat) const;

			template<typename Type>
			typename NormType<Type>::Value operator ()(const SparseMatrix<Type>& mat) const;

		};

	}
}

// Include implementation
#include "frobenius_norm.hpp"

#endif // GENUSYS_LINALG_FROBENIUS_NORM_H_
