#ifndef GENUSYS_LINALG_P_NORM_H_
#define GENUSYS_LINALG_P_NORM_H_

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<unsigned int P>
		struct PNorm {

			template<typename Type>
			struct NormType {

				typedef typename TypeTraits<typename TypeTraits<Type>::AbsType>::RealType Value;

			};

			template<typename Type>
			static typename NormType<Type>::Value norm(const Vector<Type>& vct);

			template<typename Type>
			static typename NormType<Type>::Value norm(const SparseVector<Type>& vct);
			
			template<typename Type>
			static typename NormType<Type>::Value norm(const Matrix<Type>& mat);

			template<typename Type>
			static typename NormType<Type>::Value norm(const SparseMatrix<Type>& mat);

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
#include "p_norm.hpp"

#endif // GENUSYS_LINALG_P_NORM_H_
