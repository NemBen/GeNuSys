#include <algorithm>

#include "number_traits.h"

#include "linalg_traits.h"

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		LU<Type>::LU(const Matrix<typename TypeTraits<Type>::RationalType>& L,
			const Matrix<typename TypeTraits<Type>::RationalType>& U):L(L),U(U) {
		}

		template<typename Type>
		QR<Type>::QR(const Matrix<typename TypeTraits<Type>::RealType>& Q,
			const Matrix<typename TypeTraits<Type>::RealType>& R):Q(Q),R(R) {
		}

		template<typename Type>
		HessenbergForm<Type>::HessenbergForm(const Matrix<typename TypeTraits<Type>::RealType>& Q,
			const Matrix<typename TypeTraits<Type>::RealType>& A,
			const Matrix<typename TypeTraits<Type>::RealType>& QT):Q(Q),A(A),QT(QT) {
		}

		template<typename Type>
		SchurForm<Type>::SchurForm(const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& Q,
			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& U,
			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& QT):Q(Q),U(U),QT(QT) {
		}

		template<typename Type>
		JordanForm<Type>::JordanForm(const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& P,
			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& J,
			const Matrix<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>& invP):P(P),J(J),invP(invP) {
		}

		template<typename Type>
		SmithNormalForm<Type>::SmithNormalForm(const Matrix<Type>& S,
			const Matrix<Type>& U,
			const Matrix<Type>& V):S(S),U(U),V(V) {
		}

		template<typename Type>
		typename TypeTraits<Type>::RationalType Algorithms::det(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RationalType RationalType;

			const unsigned int N = mat.rows;

			Matrix<RationalType> U = mat;
			int perm = 1;

			for (unsigned int i = 0, idxPivot = 0; i < N - 1; ++i, idxPivot += N + 1) {
				bool pivotFound = false;
				RationalType pivot = U.elem[idxPivot];
				unsigned int pivotRow = i;
				for (unsigned int j = i, idx = idxPivot; j < N; ++j, idx += N) {
					if (U.elem[idx] == NumberTraits<RationalType>::zero) {
						continue;
					}
					RationalType val = NumberTraits<RationalType>::abs(U.elem[idx]);
					if (!pivotFound || pivot < val) {
						pivotFound = true;

						pivot = val;
						pivotRow = j;
					}
				}

				if (!pivotFound) {
					return NumberTraits<RationalType>::zero;
				}

				if (pivotRow != i) {
					unsigned int idxRow1 = i * N + i;
					unsigned int idxRow2 = pivotRow * N + i;
					for (unsigned int j = i; j < N; ++j, ++idxRow1, ++idxRow2) {
						std::swap(U.elem[idxRow1], U.elem[idxRow2]);
					}
					perm = -perm;
				}

				for (unsigned int j = i + 1, idxERow = idxPivot + N; j < N; ++j, idxERow += N) {
					RationalType coef = U.elem[idxERow] / U.elem[idxPivot];
					for (unsigned int k = i, idxP = idxPivot, idxE = idxERow; k < N; ++k, ++idxP, ++idxE) {
						U.elem[idxE] -= U.elem[idxP] * coef;
					}
				}
			}

			RationalType prod = NumberTraits<RationalType>::one;
			for (unsigned int i = 0, idxU = 0; i < mat.rows; ++i, idxU += mat.cols + 1) {
				prod *= U.elem[idxU];
			}
			return perm * prod;
		}

		template<typename Type>
		Matrix<typename TypeTraits<Type>::RationalType> Algorithms::invert(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RationalType RationalType;

			const unsigned int N = mat.rows;
			Matrix<RationalType> A = mat;
			Matrix<RationalType> I = Matrix<RationalType>::identity(N, N);
			for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += N) {
				bool pivotFound = false;
				typename TypeTraits<RationalType>::AbsType pivotMax;
				unsigned int pivotRow;
				for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += N) {
					if (A.elem[idxA] == NumberTraits<RationalType>::zero) {
						continue;
					}

					typename TypeTraits<RationalType>::AbsType val = NumberTraits<RationalType>::abs(A.elem[idxA]);
					if (!pivotFound || pivotMax < val) {
						pivotFound = true;
						pivotMax = val;
						pivotRow = j;
					}
				}

				if (!pivotFound) {
					assert(false);
					break; // NOT INVERTIBLE
				}

				if (pivotRow != i) {
					unsigned int idxP = idxPivotRow;
					unsigned int idxE = pivotRow * N;
					for (unsigned int j = 0; j < i; ++j, ++idxP, ++idxE) {
						std::swap(I.elem[idxP], I.elem[idxE]);
					}
					for (unsigned int j = i; j < N; ++j, ++idxP, ++idxE) {
						std::swap(I.elem[idxP], I.elem[idxE]);
						std::swap(A.elem[idxP], A.elem[idxE]);
					}
				}

				for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + N; j < N; ++j, idxEliminateRow += N) {
					RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
					unsigned int idxP = idxPivotRow;
					unsigned int idxE = idxEliminateRow;
					for (unsigned int k = 0; k < i; ++k, ++idxP, ++idxE) {
						I.elem[idxE] -= I.elem[idxP] * coef;
					}
					for (unsigned int k = i; k < N; ++k, ++idxP, ++idxE) {
						A.elem[idxE] -= A.elem[idxP] * coef;
						I.elem[idxE] -= I.elem[idxP] * coef;
					}
				}
			}
			for (unsigned int i = N - 1, idxPivotRow = (N - 1) * N; i > 0; --i, idxPivotRow -= N) {
				for (unsigned int j = 0, idxP = idxPivotRow; j < N; ++j, ++idxP) {
					I.elem[idxP] /= A.elem[idxPivotRow + i];
				}
				for (unsigned int j = 0, idxEliminateRow = 0; j < i; ++j, idxEliminateRow += N) {
					for (unsigned int k = 0, idxP = idxPivotRow, idxE = idxEliminateRow; k < N; ++k, ++idxP, ++idxE) {
						I.elem[idxE] -= I.elem[idxP] * A.elem[idxEliminateRow + i];
					}
				}
			}
			for (unsigned int j = 0, idxP = 0; j < N; ++j, ++idxP) {
				I.elem[idxP] /= A.elem[0];
			}

			return I;
		}

#ifdef __unix__

		template<typename Type>
		Matrix<Type> Algorithms::getAdjoint(const Matrix<Type>& mat) {
			Matrix<mpz_class> mpz_mat = mat;
			Matrix<mpz_class> mpz_adj = Traits::convertUnsafe<mpq_class, mpz_class>(Algorithms::invert(mpz_mat) * Algorithms::det(mpz_mat));

			Matrix<Type> adj(mpz_adj.getRows(), mpz_adj.getCols());
			for (unsigned int i = 0; i < mpz_adj.size(); ++i) {
				adj.elem[i] = mpz_get_si(mpz_adj.elem[i].get_mpz_t());
			}

			return adj;
		}

#else

		template<typename Type>
		Matrix<Type> Algorithms::getAdjoint(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::RationalType RationalType;

			return Traits::convertUnsafe<RationalType, Type>(Algorithms::invert(mat) * Algorithms::det(mat));
		}

#endif // __unix__

		template<typename Type>
		LU<Type> Algorithms::decomposeLU(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RationalType RationalType;

			const unsigned int N = mat.rows;
			Matrix<RationalType> L = Matrix<RationalType>::identity(N, N);
			Matrix<RationalType> U(N, N);

			for (unsigned int row = 0; row < N; ++row) {
				for (unsigned int i = row; i < N; ++i) {
					RationalType sum = NumberTraits<RationalType>::zero;
					for (unsigned int j = 0; j < row; ++j) {
						sum += L.elem[row * N + j] * U.elem[j * N + i];
					}
					U.elem[row * N + i] = TypeTraits<Type>::template asType<RationalType>(mat.elem[row * N + i]) - sum;
				}
				for (unsigned int i = row + 1; i < N; ++i) {
					RationalType sum = NumberTraits<RationalType>::zero;
					for (unsigned int j = 0; j < row; ++j) {
						sum += L.elem[i * N + j] * U.elem[j * N + row];
					}
					assert(U.elem[row * N + row] != NumberTraits<RationalType>::zero);
					L.elem[i * N + row] = (TypeTraits<Type>::template asType<RationalType>(mat.elem[i * N + row]) - sum) / U.elem[row * N + row];
				}
			}

			return LU<Type>(L, U);
		}

		template<typename Type>
		Matrix<typename TypeTraits<Type>::RealType> Algorithms::getHouseholderMatrix(const Vector<Type>& vct, const unsigned int N) {
			typedef typename TypeTraits<Type>::RealType RealType;

			const unsigned int M = vct.getLength();

			Vector<RealType> houseVct = vct;
			RealType sgn = NumberTraits<RealType>::sgn(houseVct.elem[0]);
			if (sgn == NumberTraits<RealType>::zero) {
				sgn = NumberTraits<RealType>::one;
			}
			houseVct.elem[0] += sgn * PNorm<2>::norm(houseVct);
			RealType normSqr = PNorm<2>::normSqr(houseVct);

			Matrix<RealType> result = Matrix<RealType>::identity(N, N);
			if (NumberTraits<RealType>::isEpsilon(normSqr)) {
				return result;
			}
			for (unsigned int i = 0, idx = (N - M) * (N + 1); i < M; ++i, idx += N - M) {
				for (unsigned int j = 0; j < M; ++j, ++idx) {
					result.elem[idx] -= TypeTraits<int>::template asType<RealType>(2) * houseVct[i] * NumberTraits<RealType>::conj(houseVct[j]) / normSqr;
				}
			}

			return result;
		}

		template<typename Type>
		QR<Type> Algorithms::decomposeQR(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RealType RealType;

			const unsigned int N = mat.cols;
			Matrix<RealType> Q = Matrix<RealType>::identity(N, N);
			Matrix<RealType> R = mat;

			for (unsigned int col = 0; col < N - 1; ++col) {
				Matrix<RealType> transform = getHouseholderMatrix(Traits::getCol(R, col, col, N), N);

				Q = transform * Q;
				R = transform * R;
			}
			
			return QR<Type>(Q.conjugateTranspose(), R);
		}

		template<typename Type>
		HessenbergForm<Type> Algorithms::getHessenbergForm(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RealType RealType;

			const unsigned int N = mat.cols;
			Matrix<RealType> Q = Matrix<RealType>::identity(N, N);
			Matrix<RealType> A = mat;

			for (unsigned int col = 0; col < N - 1; ++col) {
				Matrix<RealType> transform = getHouseholderMatrix(Traits::getCol(A, col, col + 1, N), N);

				Q = transform * Q;
				A = transform * A * transform.conjugateTranspose();
			}
			
			return HessenbergForm<Type>(Q, A, Q.conjugateTranspose());
		}

		template<typename Type>
		SchurForm<Type> Algorithms::getSchurForm(const Matrix<Type>& mat) {
			// TODO: This method should be improved
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType ComplexRealType;
			typedef typename TypeTraits<ComplexRealType>::AbsType AbsComplexRealType;

			HessenbergForm<Type> hf = getHessenbergForm(mat);

			Matrix<ComplexRealType> Q = hf.QT;
			Matrix<ComplexRealType> H = hf.A;
			Matrix<ComplexRealType> I = Matrix<ComplexRealType>::identity(H.getRows(), H.getRows());
			
			bool flag = true;
			do {
				unsigned int shiftRow = H.getRows() - 1;
				for (unsigned int round = 0; flag && round < 100; ++round) {
					////////////////////////////////
					// Chose o near an eigenvalue //
					////////////////////////////////
					ComplexRealType h11 = H(shiftRow - 1, shiftRow - 1);
					ComplexRealType h12 = H(shiftRow - 1, shiftRow);
					ComplexRealType h21 = H(shiftRow, shiftRow - 1);
					ComplexRealType h22 = H(shiftRow, shiftRow);

					ComplexRealType sqrt = NumberTraits<ComplexRealType>::sqrt((h11 + h22) * (h11 + h22) - ComplexRealType(4, 0) * (h11 * h22 - h12 * h21));
					ComplexRealType o = NumberTraits<ComplexRealType>::div((h11 + h22) + sqrt, ComplexRealType(2,0));
					////////////////////////////////

					Matrix<ComplexRealType> shift = I * o;
					QR<ComplexRealType> qr = Algorithms::decomposeQR(H - shift);

					Q = Q * qr.Q;
					H = qr.R * qr.Q + shift;
				
					bool l = true;
					for (unsigned int j = 0, idxH = H.cols * H.rows - 2; j < H.cols - 1 && l; ++j, idxH -= H.cols + 1) {
						l = (NumberTraits<ComplexRealType>::abs(H.elem[idxH]) <= NumberTraits<AbsComplexRealType>::epsilon * NumberTraits<AbsComplexRealType>::epsilon);
						shiftRow = H.getRows() - (j + 1);
					}
					flag = !l;
				}
				
				
				for (unsigned int round = 0; flag && round < 10; ++round) {
					// Try a random shift
					unsigned int r = rand() % (H.rows - 1);
					Matrix<ComplexRealType> shift = I * H(r,r);
					QR<ComplexRealType> qr = Algorithms::decomposeQR(H - shift);

					Q = Q * qr.Q;
					H = qr.R * qr.Q + shift;

					bool l = true;
					for (unsigned int j = 0, idxH = H.cols * H.rows - 2; j < H.cols - 1 && l; ++j, idxH -= H.cols + 1) {
						l = (NumberTraits<ComplexRealType>::abs(H.elem[idxH]) <= NumberTraits<AbsComplexRealType>::epsilon * NumberTraits<AbsComplexRealType>::epsilon);
						shiftRow = H.getRows() - (j + 1);
					}
					flag = !l;
				}
			} while (flag);

			return SchurForm<Type>(Q, H, Q.conjugateTranspose());
		}

		template<typename Type>
		typename TypeTraits<typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType>::AbsType Algorithms::getSpectralRadius(const Matrix<Type>& mat) {
			typedef typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType ComplexRealType;
			typedef typename TypeTraits<ComplexRealType>::AbsType AbsType;

			SchurForm<Type> schurForm = getSchurForm(mat);
			AbsType max = NumberTraits<ComplexRealType>::abs(schurForm.U.elem[0]);
			for (unsigned int i = 1, idx = schurForm.U.cols + 2; i < schurForm.U.rows; ++i, idx += schurForm.U.cols + 1) {
				AbsType val = NumberTraits<ComplexRealType>::abs(schurForm.U.elem[idx]);
				if (val > max) {
					max = val;
				}
			}

			return max;
		}

		template<typename Type>
		unsigned int Algorithms::getRank(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::RationalType RationalType;
			typedef typename TypeTraits<RationalType>::AbsType AbsType;

			const unsigned int N = mat.rows;
			const unsigned int M = mat.cols;

			Matrix<RationalType> A = mat;
			Matrix<RationalType> P = Matrix<RationalType>::identity(N, N);

			unsigned int rank = 0;
			for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += M) {
				bool pivotFound = false;
				AbsType pivotMax = NumberTraits<AbsType>::zero;
				unsigned int pivotRow = i;
				unsigned int pivotCol = i;
				for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += i) {
					for (unsigned int k = i; k < M; ++k, ++idxA) {
						AbsType val = NumberTraits<RationalType>::abs(A.elem[idxA]);
						if (val <= NumberTraits<AbsType>::epsilon) {
							continue;
						}

						if (!pivotFound || pivotMax < val) {
							pivotFound = true;
							pivotMax = val;
							pivotRow = j;
							pivotCol = k;
						}
					}
				}
				if (!pivotFound) {
					break;
				}
				if (pivotRow != i) {
					unsigned int idxP = idxPivotRow + i;
					unsigned int idxE = pivotRow * M + i;
					for (unsigned int j = i; j < M; ++j, ++idxP, ++idxE) {
						std::swap(A.elem[idxP], A.elem[idxE]);
					}
				}
				if (pivotCol != i) {
					unsigned int idxP = i;
					unsigned int idxE = pivotCol;
					for (unsigned int j = 0; j < N; ++j, idxP += M, idxE += M) {
						std::swap(P.elem[idxP], P.elem[idxE]);
						std::swap(A.elem[idxP], A.elem[idxE]);
					}
				}

				for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + M; j < N; ++j, idxEliminateRow += M) {
					RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
					unsigned int idxP = idxPivotRow + i;
					unsigned int idxE = idxEliminateRow + i;
					for (unsigned int k = i; k < M; ++k, ++idxP, ++idxE) {
						A.elem[idxE] -= A.elem[idxP] * coef;
					}
				}

				++rank;
			}

			return rank;
		}

		template<typename Type>
		std::vector<Vector<typename TypeTraits<Type>::RealType> > Algorithms::solveHomogeneous(const Matrix<Type>& mat) {
			typedef typename TypeTraits<Type>::RationalType RationalType;
			typedef typename TypeTraits<Type>::RealType RealType;
			typedef typename TypeTraits<RationalType>::AbsType AbsType;

			const unsigned int N = mat.rows;
			const unsigned int M = mat.cols;

			Matrix<RationalType> A = mat;
			Matrix<RationalType> P = Matrix<RationalType>::identity(M, M);

			unsigned int rank = 0;
			for (unsigned int i = 0, idxPivotRow = 0; i < N; ++i, idxPivotRow += M) {
				bool pivotFound = false;
				AbsType pivotMax = NumberTraits<AbsType>::zero;
				unsigned int pivotRow = i;
				unsigned int pivotCol = i;
				for (unsigned int j = i, idxA = idxPivotRow + i; j < N; ++j, idxA += i) {
					for (unsigned int k = i; k < M; ++k, ++idxA) {
						AbsType val = NumberTraits<RationalType>::abs(A.elem[idxA]);
						if (val <= NumberTraits<AbsType>::epsilon) {
							continue;
						}

						if (!pivotFound || pivotMax < val) {
							pivotFound = true;
							pivotMax = val;
							pivotRow = j;
							pivotCol = k;
						}
					}
				}
				if (!pivotFound) {
					break;
				}
				if (pivotRow != i) {
					unsigned int idxP = idxPivotRow + i;
					unsigned int idxE = pivotRow * M + i;
					for (unsigned int j = i; j < M; ++j, ++idxP, ++idxE) {
						std::swap(A.elem[idxP], A.elem[idxE]);
					}
				}
				if (pivotCol != i) {
					unsigned int idxP = i;
					unsigned int idxE = pivotCol;
					for (unsigned int j = 0; j < N; ++j, idxP += M, idxE += M) {
						std::swap(P.elem[idxP], P.elem[idxE]);
						std::swap(A.elem[idxP], A.elem[idxE]);
					}
				}

				for (unsigned int j = i + 1, idxEliminateRow = idxPivotRow + M; j < N; ++j, idxEliminateRow += M) {
					RationalType coef = A.elem[idxEliminateRow + i] / A.elem[idxPivotRow + i];
					unsigned int idxP = idxPivotRow + i;
					unsigned int idxE = idxEliminateRow + i;
					for (unsigned int k = i; k < M; ++k, ++idxP, ++idxE) {
						A.elem[idxE] -= A.elem[idxP] * coef;
					}
				}

				++rank;
			}

			for (unsigned int i = rank - 1, idxPivotRow = (rank - 1) * M; i > 0; --i, idxPivotRow -= M) {
				RationalType coef = A.elem[idxPivotRow + i];
				for (unsigned int j = i, idxA = idxPivotRow + i; j < M; ++j, ++idxA) {
					A.elem[idxA] /= coef;
				}
				for (unsigned int j = 0, idxEliminateRow = 0; j < i; ++j, idxEliminateRow += M) {
					RationalType coef = A.elem[idxEliminateRow + i];
					for (unsigned int k = i, idxP = idxPivotRow + i, idxE = idxEliminateRow + i; k < M; ++k, ++idxP, ++idxE) {
						A.elem[idxE] -= A.elem[idxP] * coef;
					}
				}
			}
			RationalType coef = A.elem[0];
			for (unsigned int idxA = 0; idxA < M; ++idxA) {
				A.elem[idxA] /= coef;
			}

			std::vector<Vector<RealType> > result(N - rank);
			for (unsigned int i = rank; i < N; ++i) {
				Vector<RealType> eigenVector(M);
				for (unsigned int j = 0, idxA = i; j < rank; ++j, idxA += M) {
					eigenVector.elem[j] = -A.elem[idxA];
				}
				eigenVector.elem[i] = NumberTraits<RealType>::one;

				result[i - rank] = P * eigenVector;

				typename PNorm<2>::NormType<RealType>::Value norm = PNorm<2>::norm(result[i - rank]);
				Operations::vct_div(result[i - rank], TypeTraits<typename PNorm<2>::NormType<RealType>::Value>::template asType<RealType>(norm));
			}

			return result;
		}

		template<typename Type>
		JordanForm<Type> Algorithms::getJordanForm(const Matrix<Type>& mat) {
			assert(mat.cols == mat.rows);

			typedef typename TypeTraits<Type>::RealType RealType;
			typedef typename TypeTraits<typename TypeTraits<Type>::ComplexType>::RealType ComplexRealType;
			typedef typename TypeTraits<ComplexRealType>::AbsType AbsComplexRealType;

			const unsigned int N = mat.rows;

			if (N == 1) {
				return JordanForm<Type>(Matrix<RealType>::identity(1,1), mat, Matrix<RealType>::identity(1,1));
			}

			SchurForm<Type> schurForm = getSchurForm(mat);

			Matrix<ComplexRealType> M = mat;
			Matrix<ComplexRealType> I = Matrix<ComplexRealType>::identity(N, N);

			Matrix<ComplexRealType> P(N, N);
			Matrix<ComplexRealType> J(N, N);

			unsigned int idxJ = 0;
			unsigned int colP = 0;

			std::vector<ComplexRealType> eigenvalues;
			for (unsigned int i = 0; i < N; ++i) {
				eigenvalues.push_back(schurForm.U(i, i));
			}
			std::sort(eigenvalues.begin(), eigenvalues.end(), complex_comparator<RealType>);

			ComplexRealType act = eigenvalues[0];
			ComplexRealType actSum = eigenvalues[0];
			unsigned int mul = 1;
			for (unsigned int i = 1; i <= eigenvalues.size(); ++i) {
				if (i == eigenvalues.size() || NumberTraits<ComplexRealType>::abs(eigenvalues[i] - act) > 2 * NumberTraits<AbsComplexRealType>::epsilon) {
					act = NumberTraits<ComplexRealType>::div(actSum, mul);

					Matrix<ComplexRealType> base = M - I * act;
					Matrix<ComplexRealType> pow = I;

					std::vector<unsigned int> blocks;
					std::vector<Matrix<ComplexRealType> > gevBases;
					std::vector<std::vector<Vector<ComplexRealType> > > gevs;

					do {
						pow = base * pow;
						
						std::vector<Vector<ComplexRealType> > solutions = solveHomogeneous(pow);
						blocks.push_back(solutions.size());

						gevBases.push_back(pow);
						gevs.push_back(solutions);
					} while (blocks[blocks.size() - 1] < mul);

					// Calculate block size -> number of blocks
					for (unsigned int j = blocks.size() - 1; j > 0; --j) {
						blocks[j] -= blocks[j - 1];
					}
					for (unsigned int j = 0; j < blocks.size() - 1; ++j) {
						blocks[j] -= blocks[j + 1];
					}

					// Fill Jordan form with blocks
					for (unsigned int j = 0; j < blocks.size(); ++j) {
						for (unsigned int k = 0; k < blocks[j]; ++k) {
							J.elem[idxJ] = act;
							++idxJ;
							for (unsigned int l = 1; l < j+1; ++l) {
								J.elem[idxJ] = NumberTraits<ComplexRealType>::one;
								idxJ += N;
								J.elem[idxJ] = act;
								++idxJ;
							}
							idxJ += N;
						}
					}

					Matrix<ComplexRealType> T(N, mul);
					int colT = mul;
					unsigned int rankT = 0;
					for (int j = blocks.size() - 1; j >= 0; --j) {
						for (unsigned int k = 0; k < blocks[j]; ++k) {
							unsigned int selected;
							for (unsigned int l = 0; l < gevs[j].size(); ++l) {
								// NOT ZERO WITH PREVIOUS
								if (j > 0 && PNorm<00>::norm(gevBases[j - 1] * gevs[j][l]) <= NumberTraits<typename PNorm<00>::template NormType<ComplexRealType>::Value>::epsilon) {
									continue;
								}

								--colT;
								for (unsigned int m = 0, idxT = colT; m < N; ++m, idxT += mul) {
									T.elem[idxT] = gevs[j][l].elem[m];
								}
								if (getRank(T) > rankT) {
									selected = l;
									break;
								} else {
									++colT;
								}
							}

							Vector<ComplexRealType> evp = gevs[j][selected];
							for (int l = 0; l < j; ++l) {
								--colT;
								evp = base * evp;
								for (unsigned int m = 0, idxT = colT; m < N; ++m, idxT += mul) {
									T.elem[idxT] = evp.elem[m];
								}
							}

							rankT = getRank(T);
						}
					}

					for (unsigned int j = 0; j < mul; ++j) {
						for (unsigned int k = 0; k < N; ++k) {
							P.set(k, colP + j, T(k, j));
						}
					}
					colP += mul;

					// NEXT
					if (i != eigenvalues.size()) {
						act = eigenvalues[i];
						actSum = eigenvalues[i];
						mul = 1;
					}
				} else {
					++mul;
					actSum += eigenvalues[i];
				}
			}

			return JordanForm<Type>(invert(P), J, P);
		}

		template<typename Type>
		SmithNormalForm<Type> Algorithms::getSmithNormalForm(const Matrix<Type>& mat) {
			Matrix<Type> S = mat;
			Matrix<Type> U = Matrix<Type>::identity(mat.rows, mat.cols);
			Matrix<Type> V = Matrix<Type>::identity(mat.rows, mat.cols);

			for (unsigned int i = 0, idxPivot = 0; i < mat.rows - 1; ++i, idxPivot += mat.cols + 1) {
				//
				// STEP 1: PIVOT
				//
				bool pivotFound = false;
				Type min = NumberTraits<Type>::zero;
				unsigned int pivotRow = i;
				unsigned int pivotCol = i;
				for (unsigned int j = i, idx = idxPivot; j < mat.rows; ++j, idx += i) {
					for (unsigned int k = i; k < mat.cols; ++k, ++idx) {
						if (S.elem[idx] == NumberTraits<Type>::zero) {
							continue;
						}
						Type val = NumberTraits<Type>::abs(S.elem[idx]);
						if (!pivotFound || min > val) {
							pivotFound = true;

							min = val;
							pivotRow = j;
							pivotCol = k;
						}
					}
				}
				if (!pivotFound) {
					break;
				}
				if (pivotRow != i) {
					Traits::swapRows(S, pivotRow, i);
					Traits::swapRows(U, pivotRow, i);
				}
				if (pivotCol != i) {
					Traits::swapCols(S, pivotCol, i);
					Traits::swapCols(V, pivotCol, i);
				}

				//
				// STEP 2: TRANSFORM COLS / ROWS
				//
				bool flag;
				do {
					flag = false;

					// ROW TRANSFORM
					for (unsigned int j = i + 1, idxPCol = idxPivot + mat.cols; j < mat.rows; ++j, idxPCol += mat.cols) {
						if (!IntegerTraits<Type>::divisible(S.elem[idxPCol], S.elem[idxPivot])) {
							ExtendedGCD<Type> gcd = IntegerTraits<Type>::egcd(S.elem[idxPivot], S.elem[idxPCol]);

							Type coef1 = -gcd.b / gcd.gcd;
							Type coef2 = gcd.a / gcd.gcd;

							unsigned int idxR1 = i * mat.cols;
							unsigned int idxR2 = j * mat.cols;
							for (unsigned int k = 0; k < i; ++k, ++idxR1, ++idxR2) {
								Type ui = gcd.coefA * U.elem[idxR1] + gcd.coefB * U.elem[idxR2];
								Type uj = coef1 * U.elem[idxR1] + coef2 * U.elem[idxR2];
								U.elem[idxR1] = ui;
								U.elem[idxR2] = uj;
							}
							for (unsigned int k = i; k < mat.cols; ++k, ++idxR1, ++idxR2) {
								Type ui = gcd.coefA * U.elem[idxR1] + gcd.coefB * U.elem[idxR2];
								Type uj = coef1 * U.elem[idxR1] + coef2 * U.elem[idxR2];
								U.elem[idxR1] = ui;
								U.elem[idxR2] = uj;

								Type si = gcd.coefA * S.elem[idxR1] + gcd.coefB * S.elem[idxR2];
								Type sj = coef1 * S.elem[idxR1] + coef2 * S.elem[idxR2];
								S.elem[idxR1] = si;
								S.elem[idxR2] = sj;
							}

							flag = true;
						}
					}
					
					// COL TRANSFORM
					for (unsigned int j = i + 1, idxPRow = idxPivot + 1; j < mat.cols; ++j, idxPRow += 1) {
						if (!IntegerTraits<Type>::divisible(S.elem[idxPRow], S.elem[idxPivot])) {
							ExtendedGCD<Type> gcd = IntegerTraits<Type>::egcd(S.elem[idxPivot], S.elem[idxPRow]);

							Type coef1 = -gcd.b / gcd.gcd;
							Type coef2 = gcd.a / gcd.gcd;

							unsigned int idxC1 = i;
							unsigned int idxC2 = j;
							for (unsigned int k = 0; k < i; ++k, idxC1 += mat.cols, idxC2 += mat.cols) {
								Type vi = gcd.coefA * V.elem[idxC1] + gcd.coefB * V.elem[idxC2];
								Type vj = coef1 * V.elem[idxC1] + coef2 * V.elem[idxC2];
								V.elem[idxC1] = vi;
								V.elem[idxC2] = vj;
							}
							for (unsigned int k = i; k < mat.rows; ++k, idxC1 += mat.cols, idxC2 += mat.cols) {
								Type vi = gcd.coefA * V.elem[idxC1] + gcd.coefB * V.elem[idxC2];
								Type vj = coef1 * V.elem[idxC1] + coef2 * V.elem[idxC2];
								V.elem[idxC1] = vi;
								V.elem[idxC2] = vj;
								
								Type si = gcd.coefA * S.elem[idxC1] + gcd.coefB * S.elem[idxC2];
								Type sj = coef1 * S.elem[idxC1] + coef2 * S.elem[idxC2];
								S.elem[idxC1] = si;
								S.elem[idxC2] = sj;
							}

							flag = true;
						}
					}
				} while (flag);

				//
				// STEP 3: ELIMINATE COLS / ROWS
				//

				// ROW ELIMINATE
				for (unsigned int k = i + 1, idxCol = idxPivot - i + mat.cols; k < mat.rows; ++k, idxCol += mat.cols) {
					Type coef = S.elem[idxCol + i] / S.elem[idxPivot];
					unsigned int idxPRow = idxPivot - i;
					unsigned int idxERow = idxCol;
					for (unsigned int j = 0; j < i; ++j, ++idxPRow, ++idxERow) {
						U.elem[idxERow] -= U.elem[idxPRow] * coef;
					}
					for (unsigned int j = i; j < mat.cols; ++j, ++idxPRow, ++idxERow) {
						U.elem[idxERow] -= U.elem[idxPRow] * coef;
						S.elem[idxERow] -= S.elem[idxPRow] * coef;
					}
				}

				// COL ELIMINATE
				for (unsigned int k = i + 1, idxRow = i + 1; k < mat.rows; ++k, ++idxRow) {
					Type coef = S.elem[idxRow + idxPivot - i] / S.elem[idxPivot];
					unsigned int idxPCol = i;
					unsigned int idxECol = idxRow;
					for (unsigned int j = 0; j < i; ++j, idxPCol += mat.cols, idxECol += mat.cols) {
						V.elem[idxECol] -= V.elem[idxPCol] * coef;
					}
					for (unsigned int j = i; j < mat.rows; ++j, idxPCol += mat.cols, idxECol += mat.cols) {
						V.elem[idxECol] -= V.elem[idxPCol] * coef;
						S.elem[idxECol] -= S.elem[idxPCol] * coef;
					}
				}
			}

			//
			// FINAL STEP: CORRECT
			//
			bool corrected;
			do {
				corrected = false;

				for (unsigned int i = 0, idxPivot = 0; i < mat.rows - 1; ++i, idxPivot += mat.cols + 1) {
					if (S.elem[idxPivot] == NumberTraits<Type>::zero) {
						break;
					}

					if (IntegerTraits<Type>::divisible(S.elem[idxPivot + mat.cols + 1], S.elem[idxPivot])) {
						continue;
					}

					// ADD COLS (i. := i. + (i + 1).)
					for (unsigned int j = 0, idxC1 = i, idxC2 = i + 1; j < mat.rows; ++j, idxC1 += mat.cols, idxC2 += mat.cols) {
						V.elem[idxC1] += V.elem[idxC2];
					}
					S.elem[idxPivot + mat.cols] = S.elem[idxPivot + mat.cols + 1];

					// TRANSFORM
					bool flag;
					do {
						flag = false;

						// TRANSFORM ROW
						if (!IntegerTraits<Type>::divisible(S.elem[idxPivot + mat.cols], S.elem[idxPivot])) {
							ExtendedGCD<Type> gcd = IntegerTraits<Type>::egcd(S.elem[idxPivot], S.elem[idxPivot + mat.cols]);

							Type coef1 = -gcd.b / gcd.gcd;
							Type coef2 = gcd.a / gcd.gcd;

							unsigned int idxR1 = i * mat.cols;
							unsigned int idxR2 = idxR1 + mat.cols;
							for (unsigned int k = 0; k < mat.cols; ++k, ++idxR1, ++idxR2) {
								Type ui = gcd.coefA * U.elem[idxR1] + gcd.coefB * U.elem[idxR2];
								Type uj = coef1 * U.elem[idxR1] + coef2 * U.elem[idxR2];
								U.elem[idxR1] = ui;
								U.elem[idxR2] = uj;
							}
						
							Type si1 = gcd.coefA * S.elem[idxPivot] + gcd.coefB * S.elem[idxPivot + mat.cols];
							Type sj1 = coef1 * S.elem[idxPivot] + coef2 * S.elem[idxPivot + mat.cols];
							S.elem[idxPivot] = si1;
							S.elem[idxPivot + mat.cols] = sj1;

							Type si2 = gcd.coefA * S.elem[idxPivot + 1] + gcd.coefB * S.elem[idxPivot + mat.cols + 1];
							Type sj2 = coef1 * S.elem[idxPivot + 1] + coef2 * S.elem[idxPivot + mat.cols + 1];
							S.elem[idxPivot + 1] = si2;
							S.elem[idxPivot + mat.cols + 1] = sj2;

							flag = true;
						}

						// TRANSFORM COL
						if (!IntegerTraits<Type>::divisible(S.elem[idxPivot + 1], S.elem[idxPivot])) {
							ExtendedGCD<Type> gcd = IntegerTraits<Type>::egcd(S.elem[idxPivot], S.elem[idxPivot + 1]);

							Type coef1 = -gcd.b / gcd.gcd;
							Type coef2 = gcd.a / gcd.gcd;

							unsigned int idxC1 = i;
							unsigned int idxC2 = idxC1 + 1;
							for (unsigned int k = 0; k < mat.rows; ++k, idxC1 += mat.cols, idxC2 += mat.cols) {
								Type vi = gcd.coefA * V.elem[idxC1] + gcd.coefB * V.elem[idxC2];
								Type vj = coef1 * V.elem[idxC1] + coef2 * V.elem[idxC2];
								V.elem[idxC1] = vi;
								V.elem[idxC2] = vj;
							}
														
							Type si1 = gcd.coefA * S.elem[idxPivot] + gcd.coefB * S.elem[idxPivot + 1];
							Type sj1 = coef1 * S.elem[idxPivot] + coef2 * S.elem[idxPivot + 1];
							S.elem[idxPivot] = si1;
							S.elem[idxPivot + 1] = sj1;

							Type si2 = gcd.coefA * S.elem[idxPivot + mat.cols] + gcd.coefB * S.elem[idxPivot + mat.cols + 1];
							Type sj2 = coef1 * S.elem[idxPivot + mat.cols] + coef2 * S.elem[idxPivot + mat.cols + 1];
							S.elem[idxPivot + mat.cols] = si2;
							S.elem[idxPivot + mat.cols + 1] = sj2;

							flag = true;
						}

					} while (flag);

					// ELIMINATE ROW
					Type coefRow = S.elem[idxPivot + mat.cols] / S.elem[idxPivot];
					unsigned int idxPRow = i * mat.cols;
					unsigned int idxERow = idxPRow + mat.cols;
					for (unsigned int j = 0; j < mat.cols; ++j, ++idxPRow, ++idxERow) {
						U.elem[idxERow] -= U.elem[idxPRow] * coefRow;
					}
					S.elem[idxPivot + mat.cols] = 0;
					S.elem[idxPivot + mat.cols + 1] -= S.elem[idxPivot + 1] * coefRow;

					// ELIMINATE COL
					Type coefCol = S.elem[idxPivot + 1] / S.elem[idxPivot];
					unsigned int idxPCol = i;
					unsigned int idxECol = idxPCol + 1;
					for (unsigned int j = 0; j < mat.rows; ++j, idxPCol += mat.cols, idxECol += mat.cols) {
						V.elem[idxECol] -= V.elem[idxPCol] * coefCol;
					}
					S.elem[idxPivot + 1] = 0;
					S.elem[idxPivot + mat.cols + 1] -= S.elem[idxPivot + mat.cols] * coefCol;

					corrected = true;
				}

			} while (corrected);

			for (unsigned int i = 0, idxPivot = 0; i < mat.rows; ++i, idxPivot += mat.cols + 1) {
				if (S.elem[idxPivot] < 0) {
					S.elem[idxPivot] = -S.elem[idxPivot];
					for (unsigned int j = 0, idxU = i * mat.rows; j < mat.rows; ++j, ++idxU) {
						U.elem[idxU] = -U.elem[idxU];
					}
				}
			}

			return SmithNormalForm<Type>(S, U, V);
		}

	}
}
