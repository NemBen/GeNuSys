namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		typename FrobeniusNorm::template NormType<Type>::Value FrobeniusNorm::norm(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsSqrType AbsSqrType;

			AbsSqrType sum = NumberTraits<AbsSqrType>::zero;
			for (unsigned int idx = 0; idx < mat.size(); ++idx) {
				sum += NumberTraits<Type>::absSqr(mat.elem[idx]);
			}

			return NumberTraits<AbsSqrType>::sqrt(sum);
		}

		template<typename Type>
		typename FrobeniusNorm::template NormType<Type>::Value FrobeniusNorm::norm(const SparseMatrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsSqrType AbsSqrType;

			AbsSqrType sum = NumberTraits<AbsSqrType>::zero;
			for (unsigned int idx = 0; idx < mat.size(); ++idx) {
				sum += NumberTraits<AbsSqrType>::absSqr(mat.elem[idx].value);
			}

			return NumberTraits<AbsSqrType>::sqrt(sum);
		}

		template<typename Type>
		typename FrobeniusNorm::template NormType<Type>::Value FrobeniusNorm::operator ()(const Matrix<Type>& mat) const {
			return norm(mat);
		}

		template<typename Type>
		typename FrobeniusNorm::template NormType<Type>::Value FrobeniusNorm::operator ()(const SparseMatrix<Type>& mat) const {
			return norm(mat);
		}

	}
}
