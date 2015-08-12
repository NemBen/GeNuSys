#ifndef GENUSYS_LINALG_ALGORITHMS_H_
#define GENUSYS_LINALG_ALGORITHMS_H_

#include <vector>

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		struct LU {

			const Matrix<typename TypeTraits<Type>::RationalType> L;

			const Matrix<typename TypeTraits<Type>::RationalType> U;

			LU(const Matrix<typename TypeTraits<Type>::RationalType>& L,
				const Matrix<typename TypeTraits<Type>::RationalType>& U);

		};

		template<typename Type>
		struct QR {

			const Matrix<typename TypeTraits<Type>::RealType> Q;

			const Matrix<typename TypeTraits<Type>::RealType> R;

			QR(const Matrix<typename TypeTraits<Type>::RealType>& Q,
				const Matrix<typename TypeTraits<Type>::RealType>& R);

		};

		template<typename Type>
		struct HessenbergForm {

			const Matrix<typename TypeTraits<Type>::RealType> Q;

			const Matrix<typename TypeTraits<Type>::RealType> A;

			const Matrix<typename TypeTraits<Type>::RealType> QT;

			HessenbergForm(const Matrix<typename TypeTraits<Type>::RealType>& Q,
				const Matrix<typename TypeTraits<Type>::RealType>& A,
				const Matrix<typename TypeTraits<Type>::RealType>& QT);

		};

		template<typename Type>
		struct SchurForm {

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> Q;

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> U;

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> QT;

			SchurForm(const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& Q,
				const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& U,
				const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& QT);

		};

		template<typename Type>
		struct JordanForm {

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> P;

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> J;

			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType> invP;

			JordanForm(const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& P,
				const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& J,
				const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& invP);

		};
		
		template<typename Type>
		struct SmithNormalForm {

			const Matrix<Type> S;

			const Matrix<Type> U;

			const Matrix<Type> V;

			SmithNormalForm(const Matrix<Type>& S,
				const Matrix<Type>& U,
				const Matrix<Type>& V);

		};

		struct Algorithms {

			template<typename Type>
			static typename TypeTraits<Type>::RationalType det(const Matrix<Type>& mat);

			template<typename Type>
			static Matrix<typename TypeTraits<Type>::RationalType> invert(const Matrix<Type>& mat);

			template<typename Type>
			static Matrix<Type> getAdjoint(const Matrix<Type>& mat);

			template<typename Type>
			static LU<Type> decomposeLU(const Matrix<Type>& mat);

			template<typename Type>
			static Matrix<typename TypeTraits<Type>::RealType> getHouseholderMatrix(const Vector<Type>& vct, const unsigned int N);

			template<typename Type>
			static void decomposeQR(const Matrix<Type>& mat, QR<Type>& qr);

			template<typename Type>
			static QR<Type> decomposeQR(const Matrix<Type>& mat);

			template<typename Type>
			static HessenbergForm<Type> getHessenbergForm(const Matrix<Type>& mat);

			template<typename Type>
			static SchurForm<Type> getSchurForm(const Matrix<Type>& mat);

			template<typename Type>
			static typename TypeTraits<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>::AbsType getSpectralRadius(const Matrix<Type>& mat);
			
			template<typename Type>
			static unsigned int getRank(const Matrix<Type>& mat);

			template<typename Type>
			static std::vector<Vector<typename TypeTraits<Type>::RealType> > solveHomogeneous(const Matrix<Type>& mat);

			template<typename Type>
			static JordanForm<Type> getJordanForm(const Matrix<Type>& mat);

			template<typename Type>
			static SmithNormalForm<Type> getSmithNormalForm(const Matrix<Type>& mat);
			
		};
		
	}
}

// Include implementation
#include "linalg_algorithms.hpp"

#endif // GENUSYS_LINALG_ALGORITHMS_H_
