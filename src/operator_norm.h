#ifndef GENUSYS_LINALG_OPERATOR_NORM_H_
#define GENUSYS_LINALG_OPERATOR_NORM_H_

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

#include "linalg_algorithms.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename BaseType>
		class OperatorNorm {

		private:

			typedef typename TypeTraits<typename TypeTraits<BaseType>::ComplexType>::RealType ComplexRealType;

			Matrix<ComplexRealType> S;

			Matrix<ComplexRealType> invS;

		public:

			OperatorNorm(const Matrix<BaseType>& mat);

			template<typename Type>
			struct NormType {

				typedef typename PNorm<00>::template NormType<ComplexRealType>::Value Value;

			};

			template<typename Type>
			typename NormType<Type>::Value norm(const Vector<Type>& vct) const;

			template<typename Type>
			typename NormType<Type>::Value norm(const SparseVector<Type>& vct) const;

			template<typename Type>
			typename NormType<Type>::Value norm(const Matrix<Type>& mat) const;

			template<typename Type>
			typename NormType<Type>::Value norm(const SparseMatrix<Type>& mat) const;

			template<typename Type>
			typename NormType<Type>::Value operator ()(const Vector<Type>& vct) const;

			template<typename Type>
			typename NormType<Type>::Value operator ()(const SparseVector<Type>& vct) const;

			template<typename Type>
			typename NormType<Type>::Value operator ()(const Matrix<Type>& mat) const;

			template<typename Type>
			typename NormType<Type>::Value operator ()(const SparseMatrix<Type>& mat) const;

		};

	}
}

// Include implementation
#include "operator_norm.hpp"

#endif // GENUSYS_LINALG_OPERATOR_NORM_H_
