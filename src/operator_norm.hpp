#include "linalg_algorithms.h"
#include "p_norm.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename BaseType>
		OperatorNorm<BaseType>::OperatorNorm(const Matrix<BaseType>& mat) {
			typedef typename TypeTraits<ComplexRealType>::AbsType AbsType;
			typedef typename TypeTraits<BaseType>::RationalType RationalType;

			JordanForm<RationalType> jordanForm = Algorithms::getJordanForm(Algorithms::invert(mat));

			const unsigned int N = jordanForm.J.rows;

			Matrix<ComplexRealType> D = Matrix<ComplexRealType>::identity(N, N);

			bool first = true;
			ComplexRealType act;
			AbsType mu, muPow;
			for (int i = N - 1; i >= 0; --i) {
				if (first || jordanForm.J(i, i) != act) {
					mu = NumberTraits<AbsType>::one - NumberTraits<ComplexRealType>::abs(jordanForm.J(i, i));
					if (NumberTraits<AbsType>::isEpsilon(mu)) {
						mu = NumberTraits<AbsType>::one;
					}
					muPow = mu;
					first = false;
				} else {
					muPow *= mu;
				}
				D.set(i, i, muPow);
			}

			S = jordanForm.P * D;
			invS = Algorithms::invert(S);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::norm(const Vector<Type>& vct) const {
			return PNorm<00>::norm(S * Vector<ComplexRealType>(vct));
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::norm(const SparseVector<Type>& vct) const {
			return PNorm<00>::norm(S * SparseVector<ComplexRealType>(vct));
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::norm(const Matrix<Type>& mat) const {
			return PNorm<00>::norm(S * Matrix<ComplexRealType>(mat) * invS);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::norm(const SparseMatrix<Type>& mat) const {
			return PNorm<00>::norm(S * SparseMatrix<ComplexRealType>(mat) * invS);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::operator ()(const Vector<Type>& vct) const {
			return norm(vct);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::operator ()(const SparseVector<Type>& vct) const {
			return norm(vct);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::operator ()(const Matrix<Type>& mat) const {
			return norm(mat);
		}

		template<typename BaseType>
		template<typename Type>
		typename OperatorNorm<BaseType>::template NormType<Type>::Value OperatorNorm<BaseType>::operator ()(const SparseMatrix<Type>& mat) const {
			return norm(mat);
		}

	}
}