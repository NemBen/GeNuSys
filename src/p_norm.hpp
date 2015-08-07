#include <vector>

namespace GeNuSys {
	namespace LinAlg {

		template<unsigned int P>
		template<typename Type>
		typename PNorm<P>::template NormType<Type>::Value PNorm<P>::norm(const Vector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType sum = NumberTraits<AbsType>::zero;
			for (unsigned int i = 0; i < vct.length; ++i) {
				sum += NumberTraits<AbsType>::pow(NumberTraits<Type>::abs(vct.elem[i]), P);
			}

			return NumberTraits<AbsType>::root(sum, P);
		}

		template<unsigned int P>
		template<typename Type>
		typename PNorm<P>::template NormType<Type>::Value PNorm<P>::norm(const SparseVector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType sum = NumberTraits<AbsType>::zero;
			for (unsigned int i = 0; i < vct.size(); ++i) {
				sum += NumberTraits<AbsType>::pow(NumberTraits<Type>::abs(vct.elem[i].value), P);
			}

			return NumberTraits<AbsType>::root(sum, P);
		}

		template<unsigned int P>
		template<typename Type>
		typename PNorm<P>::template NormType<Type>::Value PNorm<P>::operator ()(const Vector<Type>& vct) const {
			return norm(vct);
		}

		template<unsigned int P>
		template<typename Type>
		typename PNorm<P>::template NormType<Type>::Value PNorm<P>::operator ()(const SparseVector<Type>& vct) const {
			return norm(vct);
		}

		template<>
		struct PNorm<1> {

			template<typename Type>
			struct NormType {
				
				typedef typename TypeTraits<Type>::AbsType Value;

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

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::norm(const Vector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType sum = NumberTraits<AbsType>::zero;
			for (unsigned int i = 0; i < vct.length; ++i) {
				sum += NumberTraits<Type>::abs(vct.elem[i]);
			}

			return sum;
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::norm(const SparseVector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType sum = NumberTraits<AbsType>::zero;
			for (unsigned int i = 0; i < vct.size(); ++i) {
				sum += NumberTraits<Type>::abs(vct.elem[i].value);
			}

			return sum;
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::norm(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType max = NumberTraits<AbsType>::zero;
			for (unsigned int row = 0, idx = 0; row < mat.rows; ++row, idx += mat.cols) {
				max += NumberTraits<Type>::abs(mat.elem[idx]);
			}
			for (unsigned int col = 1; col < mat.cols; ++col) {
				AbsType sum = NumberTraits<AbsType>::zero;
				for (unsigned int row = 0, idx = col; row < mat.rows; ++row, idx += mat.cols) {
					sum += NumberTraits<Type>::abs(mat.elem[idx]);
				}
				if (max < sum) {
					max = sum;
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::norm(const SparseMatrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			std::vector<AbsType> cache(mat.cols, NumberTraits<AbsType>::zero);
			for (unsigned int i = 0; i < mat.size(); ++i) {
				cache[mat.elem[i].col_idx] += NumberTraits<Type>::abs(mat.elem[i].value);
			}

			AbsType max = cache[0];
			for (unsigned int i = 1; i < mat.cols; ++i) {
				if (max < cache[i]) {
					max = cache[i];
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::operator ()(const Vector<Type>& vct) const {
			return norm(vct);
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::operator ()(const SparseVector<Type>& vct) const {
			return norm(vct);
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::operator ()(const Matrix<Type>& mat) const {
			return norm(mat);
		}

		template<typename Type>
		typename PNorm<1>::template NormType<Type>::Value PNorm<1>::operator ()(const SparseMatrix<Type>& mat) const {
			return norm(mat);
		}

		template<>
		struct PNorm<2> {

			template<typename Type>
			struct NormType {

				typedef typename TypeTraits<Type>::AbsSqrType SqrType;

				typedef typename TypeTraits<SqrType>::RealType Value;

			};
			
			template<typename Type>
			static typename NormType<Type>::SqrType normSqr(const Vector<Type>& vct);
			
			template<typename Type>
			static typename NormType<Type>::SqrType normSqr(const SparseVector<Type>& vct);
			
			template<typename Type>
			static typename NormType<Type>::Value norm(const Vector<Type>& vct);
			
			template<typename Type>
			static typename NormType<Type>::Value norm(const SparseVector<Type>& vct);
			
			template<typename Type>
			typename NormType<Type>::Value operator ()(const Vector<Type>& vct) const;
			
			template<typename Type>
			typename NormType<Type>::Value operator ()(const SparseVector<Type>& vct) const;

		};

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::SqrType PNorm<2>::normSqr(const Vector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsSqrType AbsSqrType;

			AbsSqrType sum = NumberTraits<AbsSqrType>::zero;
			for (unsigned int i = 0; i < vct.length; ++i) {
				sum += NumberTraits<Type>::absSqr(vct.elem[i]);
			}

			return sum;
		}

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::SqrType PNorm<2>::normSqr(const SparseVector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsSqrType AbsSqrType;

			AbsSqrType sum = NumberTraits<AbsSqrType>::zero;
			for (unsigned int i = 0; i < vct.size(); ++i) {
				sum += NumberTraits<Type>::absSqr(vct.elem[i].value);
			}

			return sum;
		}

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::Value PNorm<2>::norm(const Vector<Type>& vct) {
			return NumberTraits<typename NormType<Type>::SqrType>::sqrt(normSqr(vct));
		}

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::Value PNorm<2>::norm(const SparseVector<Type>& vct) {
			return NumberTraits<typename NormType<Type>::SqrType>::sqrt(normSqr(vct));
		}

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::Value PNorm<2>::operator ()(const Vector<Type>& vct) const {
			return norm(vct);
		}

		template<typename Type>
		typename PNorm<2>::template NormType<Type>::Value PNorm<2>::operator ()(const SparseVector<Type>& vct) const {
			return norm(vct);
		}

		template<>
		struct PNorm<00> {

			template<typename Type>
			struct NormType {

				typedef typename TypeTraits<Type>::AbsType Value;

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

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::norm(const Vector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType max = NumberTraits<Type>::abs(vct.elem[0]);
			for (unsigned int i = 1; i < vct.length; ++i) {
				AbsType abs = NumberTraits<Type>::abs(vct.elem[i]);
				if (abs > max) {
					max = abs;
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::norm(const SparseVector<Type>& vct) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			if (vct.size() == 0) {
				return NumberTraits<AbsType>::zero;
			}

			AbsType max = NumberTraits<Type>::abs(vct.elem[0].value);
			for (unsigned int i = 1; i < vct.size(); ++i) {
				AbsType abs = NumberTraits<Type>::abs(vct.elem[i].value);
				if (abs > max) {
					max = abs;
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::norm(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			unsigned int idx = 0;
			AbsType max = NumberTraits<AbsType>::zero;
			for (unsigned int col = 0; col < mat.cols; ++col) {
				max += NumberTraits<Type>::abs(mat.elem[idx++]);
			}
			for (unsigned int row = 1; row < mat.rows; ++row) {
				AbsType sum = NumberTraits<AbsType>::zero;
				for (unsigned int col = 0; col < mat.cols; ++col) {
					sum += NumberTraits<Type>::abs(mat.elem[idx++]);
				}
				if (max < sum) {
					max = sum;
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::norm(const SparseMatrix<Type>& mat) {
			typedef typename TypeTraits<Type>::AbsType AbsType;

			AbsType max = NumberTraits<AbsType>::zero;
			for (unsigned int idx = mat.row_ptr[0]; idx < mat.row_ptr[1]; ++idx) {
				max += NumberTraits<Type>::abs(mat.elem[idx].value);
			}
			for (unsigned int row = 1; row < mat.rows; ++row) {
				AbsType sum = NumberTraits<AbsType>::zero;
				for (unsigned int idx = mat.row_ptr[row]; idx < mat.row_ptr[row + 1]; ++idx) {
					sum += NumberTraits<Type>::abs(mat.elem[idx].value);
				}
				if (max < sum) {
					max = sum;
				}
			}

			return max;
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::operator ()(const Vector<Type>& vct) const {
			return norm(vct);
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::operator ()(const SparseVector<Type>& vct) const {
			return norm(vct);
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::operator ()(const Matrix<Type>& mat) const {
			return norm(mat);
		}

		template<typename Type>
		typename PNorm<00>::template NormType<Type>::Value PNorm<00>::operator ()(const SparseMatrix<Type>& mat) const {
			return norm(mat);
		}

	}
}