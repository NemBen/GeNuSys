#ifndef GENUSYS_LINALG_OPERATIONS_H_
#define GENUSYS_LINALG_OPERATIONS_H_

#include "number_traits.h"

#include "vector.h"
#include "sparse_vector.h"
#include "matrix.h"
#include "sparse_matrix.h"

namespace GeNuSys {
	namespace LinAlg {

		struct Operations {

			template<typename Type>
			static void vct_mul(const Vector<Type>& vct, const Type& value, Vector<Type>& result);

			template<typename Type>
			static void vct_mul(Vector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_mul(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result);

			template<typename Type>
			static void vct_mul(SparseVector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_idiv(const Vector<Type>& vct, const Type& value, Vector<Type>& result);

			template<typename Type>
			static void vct_idiv(Vector<Type>& vct, const Type& value);
			
			template<typename Type>
			static void vct_idiv(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result);
			
			template<typename Type>
			static void vct_idiv(SparseVector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_mod(const Vector<Type>& vct, const Type& value, Vector<Type>& result);

			template<typename Type>
			static void vct_mod(Vector<Type>& vct, const Type& value);
			
			template<typename Type>
			static void vct_mod(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result);
			
			template<typename Type>
			static void vct_mod(SparseVector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_mods(const Vector<Type>& vct, const Type& value, Vector<Type>& result);

			template<typename Type>
			static void vct_mods(Vector<Type>& vct, const Type& value);
			
			template<typename Type>
			static void vct_mods(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result);
			
			template<typename Type>
			static void vct_mods(SparseVector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_div(const Vector<Type>& vct, const Type& value, Vector<typename TypeTraits<Type>::RationalType>& result);

			template<typename Type>
			static void vct_div(Vector<Type>& vct, const Type& value);
			
			template<typename Type>
			static void vct_div(const SparseVector<Type>& vct, const Type& value, SparseVector<typename TypeTraits<Type>::RationalType>& result);
			
			template<typename Type>
			static void vct_div(SparseVector<Type>& vct, const Type& value);

			template<typename Type>
			static void vct_add(const Vector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result);

			template<typename Type>
			static void vct_add(Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static void vct_add(const Vector<Type>& op1, const SparseVector<Type>& op2, Vector<Type>& result);

			template<typename Type>
			static void vct_add(Vector<Type>& op1, const SparseVector<Type>& op2);
			
			template<typename Type>
			static void vct_add(const SparseVector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result);

			template<typename Type>
			static void vct_add(const SparseVector<Type>& op1, const SparseVector<Type>& op2, SparseVector<Type>& result);

			template<typename Type>
			static void vct_sub(const Vector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result);

			template<typename Type>
			static void vct_sub(Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static void vct_sub(const Vector<Type>& op1, const SparseVector<Type>& op2, Vector<Type>& result);

			template<typename Type>
			static void vct_sub(Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static void vct_sub(const SparseVector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result);
			
			template<typename Type>
			static void vct_sub(const SparseVector<Type>& op1, const SparseVector<Type>& op2, SparseVector<Type>& result);

			template<typename Type>
			static Type vct_mul(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static Type vct_mul(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static Type vct_mul(const SparseVector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static Type vct_mul(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_eq(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_eq(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_eq(const SparseVector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_eq(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_gt(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_gt(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_gt(const SparseVector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_gt(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_geq(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_geq(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_geq(const SparseVector<Type>& op1, const Vector<Type>& op2);
			
			template<typename Type>
			static bool vct_geq(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_lt(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_lt(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_lt(const SparseVector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_lt(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_leq(const Vector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_leq(const Vector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static bool vct_leq(const SparseVector<Type>& op1, const Vector<Type>& op2);

			template<typename Type>
			static bool vct_leq(const SparseVector<Type>& op1, const SparseVector<Type>& op2);

			template<typename Type>
			static void mat_transpose(const Matrix<Type>& mat, Matrix<Type>& result);

			template<typename Type>
			static void mat_transpose(const SparseMatrix<Type>& mat, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_conjugate_transpose(const Matrix<Type>& mat, Matrix<Type>& result);

			template<typename Type>
			static void mat_conjugate_transpose(const SparseMatrix<Type>& mat, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result);

			template<typename Type>
			static void mat_mul(Matrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(SparseMatrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_idiv(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result);

			template<typename Type>
			static void mat_idiv(Matrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_idiv(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_idiv(SparseMatrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mod(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result);

			template<typename Type>
			static void mat_mod(Matrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mod(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mod(SparseMatrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mods(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result);

			template<typename Type>
			static void mat_mods(Matrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mods(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mods(SparseMatrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_div(const Matrix<Type>& mat, const Type& value, Matrix<typename TypeTraits<Type>::RationalType>& result);

			template<typename Type>
			static void mat_div(Matrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_div(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<typename TypeTraits<Type>::RationalType>& result);

			template<typename Type>
			static void mat_div(SparseMatrix<Type>& mat, const Type& value);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& mat, const Vector<Type>& vct, Vector<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& mat, const Vector<Type>& vct, SparseVector<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& mat, const SparseVector<Type>& vct, Vector<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& mat, const SparseVector<Type>& vct, SparseVector<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& mat, const Vector<Type>& vct, Vector<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& mat, const Vector<Type>& vct, SparseVector<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& mat, const SparseVector<Type>& vct, Vector<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& mat, const SparseVector<Type>& vct, SparseVector<Type>& result);

			template<typename Type>
			static void mat_add(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_add(Matrix<Type>& op1, const Matrix<Type>& op2);

			template<typename Type>
			static void mat_add(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_add(Matrix<Type>& op1, const SparseMatrix<Type>& op2);

			template<typename Type>
			static void mat_add(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_add(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_sub(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_sub(Matrix<Type>& op1, const Matrix<Type>& op2);

			template<typename Type>
			static void mat_sub(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_sub(Matrix<Type>& op1, const SparseMatrix<Type>& op2);

			template<typename Type>
			static void mat_sub(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_sub(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& op1, const Matrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_mul(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result);

			template<typename Type>
			static void mat_mul(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result);

			template<typename Type>
			static bool mat_eq(const Matrix<Type>& op1, const Matrix<Type>& op2);

			template<typename Type>
			static bool mat_eq(const Matrix<Type>& op1, const SparseMatrix<Type>& op2);

			template<typename Type>
			static bool mat_eq(const SparseMatrix<Type>& op1, const Matrix<Type>& op2);

			template<typename Type>
			static bool mat_eq(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2);

		};

	}
}

#include "linalg_operations.hpp"

#endif // GENUSYS_LINALG_OPERATIONS_H_
