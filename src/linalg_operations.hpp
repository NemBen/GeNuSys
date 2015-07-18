#include <cassert>
#include <algorithm>
#include <vector>

namespace GeNuSys {
	namespace LinAlg {

		template<typename Type>
		void Operations::vct_mul(const Vector<Type>& vct, const Type& value, Vector<Type>& result) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = vct.elem[i] * value;
			}
		}

		template<typename Type>
		void Operations::vct_mul(Vector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				vct.elem[i] *= value;
			}
		}

		template<typename Type>
		void Operations::vct_mul(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result) {
			result.elem.resize(vct.size());

			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<Type>::Entry(vct.elem[i].idx, vct.elem[i].value * value);
			}
		}

		template<typename Type>
		void Operations::vct_mul(SparseVector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				vct.elem[i].value *= value;
			}
		}

		template<typename Type>
		void Operations::vct_idiv(const Vector<Type>& vct, const Type& value, Vector<Type>& result) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = IntegerTraits<Type>::div(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_idiv(Vector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				vct.elem[i] = IntegerTraits<Type>::div(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_idiv(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result) {
			result.elem.resize(vct.size());

			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<Type>::Entry(vct.elem[i].idx, IntegerTraits<Type>::div(vct.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::vct_idiv(SparseVector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				vct.elem[i].value = IntegerTraits<Type>::div(vct.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::vct_mod(const Vector<Type>& vct, const Type& value, Vector<Type>& result) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = IntegerTraits<Type>::mod(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_mod(Vector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				vct.elem[i] = IntegerTraits<Type>::mod(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_mod(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result) {
			result.elem.resize(vct.size());

			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<Type>::Entry(vct.elem[i].idx, IntegerTraits<Type>::mod(vct.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::vct_mod(SparseVector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				vct.elem[i] = IntegerTraits<Type>::mod(vct.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::vct_mods(const Vector<Type>& vct, const Type& value, Vector<Type>& result) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = IntegerTraits<Type>::mods(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_mods(Vector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				vct.elem[i] = IntegerTraits<Type>::mods(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_mods(const SparseVector<Type>& vct, const Type& value, SparseVector<Type>& result) {
			result.elem.resize(vct.size());

			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<Type>::Entry(vct.elem[i].idx, IntegerTraits<Type>::mods(vct.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::vct_mods(SparseVector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				vct.elem[i].value = IntegerTraits<Type>::mods(vct.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::vct_div(const Vector<Type>& vct, const Type& value, Vector<typename TypeTraits<Type>::RationalType>& result) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				result.elem[i] = NumberTraits<Type>::div(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_div(Vector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.length; ++i) {
				vct.elem[i] = NumberTraits<Type>::div(vct.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::vct_div(const SparseVector<Type>& vct, const Type& value, SparseVector<typename TypeTraits<Type>::RationalType>& result) {
			result.elem.resize(vct.size());

			for (unsigned int i = 0; i < vct.size(); ++i) {
				result.elem[i] = typename SparseVector<typename TypeTraits<Type>::RationalType>::Entry(vct.elem[i].idx, NumberTraits<Type>::div(vct.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::vct_div(SparseVector<Type>& vct, const Type& value) {
			for (unsigned int i = 0; i < vct.size(); ++i) {
				vct.elem[i].value = NumberTraits<Type>::div(vct.elem[i].value, value);
			}
		}
			
		template<typename Type>
		void Operations::vct_add(const Vector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				result.elem[i] = op1.elem[i] + op2.elem[i];
			}
		}
			
		template<typename Type>
		void Operations::vct_add(Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				op1.elem[i] += op2.elem[i];
			}
		}

		template<typename Type>
		void Operations::vct_add(const Vector<Type>& op1, const SparseVector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				result.elem[i] = op1.elem[i];
			}
			for (unsigned int i = 0; i < op2.size(); ++i) {
				result.elem[op2.elem[i].idx] += op2.elem[i].value;
			}
		}

		template<typename Type>
		void Operations::vct_add(Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op2.size(); ++i) {
				op1.elem[op2.elem[i].idx] += op2.elem[i].value;
			}
		}

		template<typename Type>
		void Operations::vct_add(const SparseVector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op2.length; ++i) {
				result.elem[i] = op2.elem[i];
			}
			for (unsigned int i = 0; i < op1.size(); ++i) {
				result.elem[op1.elem[i].idx] += op1.elem[i].value;
			}
		}
			
		template<typename Type>
		void Operations::vct_add(const SparseVector<Type>& op1, const SparseVector<Type>& op2, SparseVector<Type>& result) {
			assert(op1.length == op2.length);

			result.elem.clear();

			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size()) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					result.push(op1.elem[i].idx, op1.elem[i].value);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					result.push(op2.elem[j].idx, op2.elem[j].value);
					++j;
				} else {
					result.push(op1.elem[i].idx, op1.elem[i].value + op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size()) {
				result.push(op1.elem[i].idx, op1.elem[i].value);
				++i;
			}
			while (j < op2.size()) {
				result.push(op2.elem[j].idx, op2.elem[j].value);
				++j;
			}
		}

		template<typename Type>
		void Operations::vct_sub(const Vector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				result.elem[i] = op1.elem[i] - op2.elem[i];
			}
		}

		template<typename Type>
		void Operations::vct_sub(Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				op1.elem[i] -= op2.elem[i];
			}
		}

		template<typename Type>
		void Operations::vct_sub(const Vector<Type>& op1, const SparseVector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op1.length; ++i) {
				result.elem[i] = op1.elem[i];
			}
			for (unsigned int i = 0; i < op2.size(); ++i) {
				result.elem[op2.elem[i].idx] -= op2.elem[i].value;
			}
		}

		template<typename Type>
		void Operations::vct_sub(Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op2.size(); ++i) {
				op1.elem[op2.elem[i].idx] -= op2.elem[i].value;
			}
		}

		template<typename Type>
		void Operations::vct_sub(const SparseVector<Type>& op1, const Vector<Type>& op2, Vector<Type>& result) {
			assert(op1.length == op2.length);

			for (unsigned int i = 0; i < op2.length; ++i) {
				result.elem[i] = -op2.elem[i];
			}
			for (unsigned int i = 0; i < op1.size(); ++i) {
				result.elem[op1.elem[i].idx] += op1.elem[i].value;
			}
		}

		template<typename Type>
		void Operations::vct_sub(const SparseVector<Type>& op1, const SparseVector<Type>& op2, SparseVector<Type>& result) {
			assert(op1.length == op2.length);

			result.elem.clear();

			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size()) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					result.push(op1.elem[i].idx, op1.elem[i].value);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					result.push(op2.elem[j].idx, -op2.elem[j].value);
					++j;
				} else {
					result.push(op1.elem[i].idx, op1.elem[i].value - op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size()) {
				result.push(op1.elem[i].idx, op1.elem[i].value);
				++i;
			}
			while (j < op2.size()) {
				result.push(op2.elem[j].idx, -op2.elem[j].value);
				++j;
			}
		}

		template<typename Type>
		Type Operations::vct_mul(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			Type result = NumberTraits<Type>::zero;
			for (unsigned int i = 0; i < op1.length; ++i) {
				result += op1.elem[i] * op2.elem[i];
			}

			return result;
		}

		template<typename Type>
		Type Operations::vct_mul(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			Type result = NumberTraits<Type>::zero;
			for (unsigned int i = 0; i < op2.size(); ++i) {
				result += op1.elem[op2.elem[i].idx] * op2.elem[i].value;
			}

			return result;
		}

		template<typename Type>
		Type Operations::vct_mul(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			Type result = NumberTraits<Type>::zero;
			for (unsigned int i = 0; i < op1.size(); ++i) {
				result += op2.elem[op1.elem[i].idx] * op1.elem[i].value;
			}

			return result;
		}

		template<typename Type>
		Type Operations::vct_mul(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			Type result = NumberTraits<Type>::zero;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size()) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					++j;
				} else {
					result += op1.elem[i].value * op2.elem[j].value;
					++i;
					++j;
				}
			}

			return result;
		}

		template<typename Type>
		bool Operations::vct_eq(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			for (unsigned int i = 0; i < op1.length && l; ++i) {
				l = (op1.elem[i] == op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_eq(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.length && j < op2.size() && l) {
				if (i < op2.elem[j].idx) {
					l = (op1.elem[i] == NumberTraits<Type>::zero);
					++i;
				} else {
					l = (op1.elem[i] == op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.length && l) {
				l = (op1.elem[i] == NumberTraits<Type>::zero);
				++i;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_eq(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.length && l) {
				if (j < op1.elem[i].idx) {
					l = (NumberTraits<Type>::zero == op2.elem[j]);
					++j;
				} else {
					l = (op1.elem[i].value == op2.elem[j]);
					++i;
					++j;
				}
			}
			while (j < op2.length && l) {
				l = (NumberTraits<Type>::zero == op2.elem[j]);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_eq(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = (op1.size() == op2.size());
			for (unsigned int i = 0; i < op1.size() && l; ++i) {
				l = (op1.elem[i].idx == op2.elem[i].idx && op1.elem[i].value == op2.elem[i].value);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_gt(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			for (unsigned int i = 0; i < op1.length && l; ++i) {
				l = (op1.elem[i] > op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_gt(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.length && j < op2.size() && l) {
				if (i < op2.elem[j].idx) {
					l = (op1.elem[i] > NumberTraits<Type>::zero);
					++i;
				} else {
					l = (op1.elem[i] > op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.length && l) {
				l = (op1.elem[i] > NumberTraits<Type>::zero);
				++i;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_gt(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.length && l) {
				if (j < op1.elem[i].idx) {
					l = (NumberTraits<Type>::zero > op2.elem[j]);
					++j;
				} else {
					l = (op1.elem[i].value > op2.elem[j]);
					++i;
					++j;
				}
			}
			while (j < op2.length && l) {
				l = (NumberTraits<Type>::zero > op2.elem[j]);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_gt(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size() && l) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					l = (op1.elem[i].value > NumberTraits<Type>::zero);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					l = (NumberTraits<Type>::zero > op2.elem[j].value);
					++j;
				} else {
					l = (op1.elem[i].value > op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size() && l) {
				l = (op1.elem[i].value > NumberTraits<Type>::zero);
				++i;
			}
			while (j < op2.size() && l) {
				l = (NumberTraits<Type>::zero > op2.elem[j].value);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_geq(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			for (unsigned int i = 0; i < op1.length && l; ++i) {
				l = (op1.elem[i] >= op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_geq(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.length && j < op2.size() && l) {
				if (i < op2.elem[j].idx) {
					l = (op1.elem[i] >= NumberTraits<Type>::zero);
					++i;
				} else {
					l = (op1.elem[i] >= op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.length && l) {
				l = (op1.elem[i] >= NumberTraits<Type>::zero);
				++i;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_geq(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.length && l) {
				if (j < op1.elem[i].idx) {
					l = (NumberTraits<Type>::zero >= op2.elem[j]);
					++j;
				} else {
					l = (op1.elem[i].value >= op2.elem[j]);
					++i;
					++j;
				}
			}
			while (j < op2.length && l) {
				l = (NumberTraits<Type>::zero >= op2.elem[j]);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_geq(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size() && l) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					l = (op1.elem[i].value >= NumberTraits<Type>::zero);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					l = (NumberTraits<Type>::zero >= op2.elem[j].value);
					++j;
				} else {
					l = (op1.elem[i].value >= op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size() && l) {
				l = (op1.elem[i].value >= NumberTraits<Type>::zero);
				++i;
			}
			while (j < op2.size() && l) {
				l = (NumberTraits<Type>::zero >= op2.elem[j].value);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_lt(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			for (unsigned int i = 0; i < op1.length && l; ++i) {
				l = (op1.elem[i] < op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_lt(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.length && j < op2.size() && l) {
				if (i < op2.elem[j].idx) {
					l = (op1.elem[i] < NumberTraits<Type>::zero);
					++i;
				} else {
					l = (op1.elem[i] < op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.length && l) {
				l = (op1.elem[i] < NumberTraits<Type>::zero);
				++i;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_lt(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.length && l) {
				if (j < op1.elem[i].idx) {
					l = (NumberTraits<Type>::zero < op2.elem[j]);
					++j;
				} else {
					l = (op1.elem[i].value < op2.elem[j]);
					++i;
					++j;
				}
			}
			while (j < op2.length && l) {
				l = (NumberTraits<Type>::zero < op2.elem[j]);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_lt(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size() && l) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					l = (op1.elem[i].value < NumberTraits<Type>::zero);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					l = (NumberTraits<Type>::zero < op2.elem[j].value);
					++j;
				} else {
					l = (op1.elem[i].value < op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size() && l) {
				l = (op1.elem[i].value < NumberTraits<Type>::zero);
				++i;
			}
			while (j < op2.size() && l) {
				l = (NumberTraits<Type>::zero < op2.elem[j].value);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_leq(const Vector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			for (unsigned int i = 0; i < op1.length && l; ++i) {
				l = (op1.elem[i] <= op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_leq(const Vector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.length && j < op2.size() && l) {
				if (i < op2.elem[j].idx) {
					l = (op1.elem[i] <= NumberTraits<Type>::zero);
					++i;
				} else {
					l = (op1.elem[i] <= op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.length && l) {
				l = (op1.elem[i] <= NumberTraits<Type>::zero);
				++i;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_leq(const SparseVector<Type>& op1, const Vector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.length && l) {
				if (j < op1.elem[i].idx) {
					l = (NumberTraits<Type>::zero <= op2.elem[j]);
					++j;
				} else {
					l = (op1.elem[i].value <= op2.elem[j]);
					++i;
					++j;
				}
			}
			while (j < op2.length && l) {
				l = (NumberTraits<Type>::zero <= op2.elem[j]);
				++j;
			}

			return l;
		}

		template<typename Type>
		bool Operations::vct_leq(const SparseVector<Type>& op1, const SparseVector<Type>& op2) {
			assert(op1.length == op2.length);

			bool l = true;
			unsigned int i = 0, j = 0;
			while (i < op1.size() && j < op2.size() && l) {
				if (op1.elem[i].idx < op2.elem[j].idx) {
					l = (op1.elem[i].value <= NumberTraits<Type>::zero);
					++i;
				} else if (op1.elem[i].idx > op2.elem[j].idx) {
					l = (NumberTraits<Type>::zero <= op2.elem[j].value);
					++j;
				} else {
					l = (op1.elem[i].value <= op2.elem[j].value);
					++i;
					++j;
				}
			}
			while (i < op1.size() && l) {
				l = (op1.elem[i].value <= NumberTraits<Type>::zero);
				++i;
			}
			while (j < op2.size() && l) {
				l = (NumberTraits<Type>::zero <= op2.elem[j].value);
				++j;
			}

			return l;
		}
		
		template<typename Type>
		void Operations::mat_transpose(const Matrix<Type>& mat, Matrix<Type>& result) {
			for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i) {
				for (unsigned int j = 0, idxR = i; j < mat.cols; ++j, ++idxA, idxR += mat.rows) {
					result.elem[idxR] = mat.elem[idxA];
				}
			}
		}
		
		template<typename Type>
		void Operations::mat_transpose(const SparseMatrix<Type>& mat, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());
			std::fill(result.row_ptr.begin(), result.row_ptr.end(), NumberTraits<Type>::zero);

			for (unsigned int i = 0; i < mat.size(); ++i) {
				++result.row_ptr[mat.elem[i].col_idx + 1];
			}
			for (unsigned int i = 1; i <= result.rows; ++i) {
				result.row_ptr[i] += result.row_ptr[i-1];
			}
			std::vector<unsigned int> inserted(result.rows, 0);
			for (unsigned int i = 0; i < mat.rows; ++i) {
				for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i+1]; ++j) {
					result.elem[result.row_ptr[mat.elem[j].col_idx] + inserted[mat.elem[j].col_idx]++] =
						typename SparseMatrix<Type>::Entry(i, mat.elem[j].value);
				}
			}
		}

		template<typename Type>
		void Operations::mat_conjugate_transpose(const Matrix<Type>& mat, Matrix<Type>& result) {
			for (unsigned int i = 0, idxA = 0; i < mat.rows; ++i) {
				for (unsigned int j = 0, idxR = i; j < mat.cols; ++j, ++idxA, idxR += mat.rows) {
					result.elem[idxR] = NumberTraits<Type>::conj(mat.elem[idxA]);
				}
			}
		}
		
		template<typename Type>
		void Operations::mat_conjugate_transpose(const SparseMatrix<Type>& mat, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());
			std::fill(result.row_ptr.begin(), result.row_ptr.end(), NumberTraits<Type>::zero);

			for (unsigned int i = 0; i < mat.size(); ++i) {
				++result.row_ptr[mat.elem[i].col_idx + 1];
			}
			for (unsigned int i = 1; i <= result.rows; ++i) {
				result.row_ptr[i] += result.row_ptr[i-1];
			}
			std::vector<unsigned int> inserted(result.rows, 0);
			for (unsigned int i = 0; i < mat.rows; ++i) {
				for (unsigned int j = mat.row_ptr[i]; j < mat.row_ptr[i+1]; ++j) {
					result.elem[result.row_ptr[mat.elem[j].col_idx] + inserted[mat.elem[j].col_idx]++] =
						typename SparseMatrix<Type>::Entry(i, NumberTraits<Type>::conj(mat.elem[j].value));
				}
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = mat.elem[i] * value;
			}
		}

		template<typename Type>
		void Operations::mat_mul(Matrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i] *= value;
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());

			result.row_ptr = mat.row_ptr;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<Type>::Entry(mat.elem[i].col_idx, mat.elem[i].value * value);
			}
		}

		template<typename Type>
		void Operations::mat_mul(SparseMatrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i].value *= value;
			}
		}

		template<typename Type>
		void Operations::mat_idiv(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = IntegerTraits<Type>::div(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_idiv(Matrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i] = IntegerTraits<Type>::div(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_idiv(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());

			result.row_ptr = mat.row_ptr;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<Type>::Entry(mat.elem[i].col_idx, IntegerTraits<Type>::div(mat.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::mat_idiv(SparseMatrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i].value = IntegerTraits<Type>::div(mat.elem[i].value, value);
			}
		}
		
		template<typename Type>
		void Operations::mat_mod(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = IntegerTraits<Type>::mod(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_mod(Matrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i] = IntegerTraits<Type>::mod(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_mod(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());

			result.row_ptr = mat.row_ptr;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<Type>::Entry(mat.elem[i].col_idx, IntegerTraits<Type>::mod(mat.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::mat_mod(SparseMatrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i].value = IntegerTraits<Type>::mod(mat.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::mat_mods(const Matrix<Type>& mat, const Type& value, Matrix<Type>& result) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = IntegerTraits<Type>::mods(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_mods(Matrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i] = IntegerTraits<Type>::mods(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_mods(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<Type>& result) {
			result.elem.resize(mat.size());

			result.row_ptr = mat.row_ptr;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<Type>::Entry(mat.elem[i].col_idx, IntegerTraits<Type>::mods(mat.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::mat_mods(SparseMatrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i].value = IntegerTraits<Type>::mods(mat.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::mat_div(const Matrix<Type>& mat, const Type& value, Matrix<typename TypeTraits<Type>::RationalType>& result) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = NumberTraits<Type>::div(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_div(Matrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem[i] = NumberTraits<Type>::div(mat.elem[i], value);
			}
		}

		template<typename Type>
		void Operations::mat_div(const SparseMatrix<Type>& mat, const Type& value, SparseMatrix<typename TypeTraits<Type>::RationalType>& result) {
			result.elem.resize(mat.size());

			result.row_ptr = mat.row_ptr;
			for (unsigned int i = 0; i < mat.size(); ++i) {
				result.elem[i] = typename SparseMatrix<Type>::Entry(mat.elem[i].col_idx, NumberTraits<Type>::div(mat.elem[i].value, value));
			}
		}

		template<typename Type>
		void Operations::mat_div(SparseMatrix<Type>& mat, const Type& value) {
			for (unsigned int i = 0; i < mat.size(); ++i) {
				mat.elem.value[i] = NumberTraits<Type>::div(mat.elem[i].value, value);
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& mat, const Vector<Type>& vct, Vector<Type>& result) {
			assert(mat.cols == vct.length);

			for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxB = 0; idxB < vct.length; ++idxB, ++idxA) {
					prod += mat.elem[idxA] * vct.elem[idxB];
				}
				result.elem[idxR] = prod;
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& mat, const Vector<Type>& vct, SparseVector<Type>& result) {
			assert(mat.cols == vct.length);

			result.elem.clear();

			for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxB = 0; idxB < vct.length; ++idxB, ++idxA) {
					prod += mat.elem[idxA] * vct.elem[idxB];
				}
				result.push(idxR, prod);
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& mat, const SparseVector<Type>& vct, Vector<Type>& result) {
			assert(mat.cols == vct.length);

			for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR, idxA += mat.cols) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxB = 0; idxB < vct.size(); ++idxB) {
					prod += mat.elem[idxA + vct.elem[idxB].idx] * vct.elem[idxB].value;
				}
				result.elem[idxR] = prod;
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& mat, const SparseVector<Type>& vct, SparseVector<Type>& result) {
			assert(mat.cols == vct.length);

			result.elem.clear();

			for (unsigned int idxA = 0, idxR = 0; idxA < mat.size(); ++idxR, idxA += mat.cols) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxB = 0; idxB < vct.size(); ++idxB) {
					prod += mat.elem[idxA + vct.elem[idxB].idx] * vct.elem[idxB].value;
				}
				result.push(idxR, prod);
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& mat, const Vector<Type>& vct, Vector<Type>& result) {
			assert(mat.cols == vct.length);

			for (unsigned int row = 0; row < mat.rows; ++row) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxA = mat.row_ptr[row]; idxA < mat.row_ptr[row + 1]; ++idxA) {
					prod += mat.elem[idxA].value * vct.elem[mat.elem[idxA].col_idx];
				}
				result.elem[row] = prod;
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& mat, const Vector<Type>& vct, SparseVector<Type>& result) {
			assert(mat.cols == vct.length);

			result.elem.clear();

			for (unsigned int row = 0; row < mat.rows; ++row) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxA = mat.row_ptr[row]; idxA < mat.row_ptr[row + 1]; ++idxA) {
					prod += mat.elem[idxA].value * vct.elem[mat.elem[idxA].col_idx];
				}
				result.push(row, prod);
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& mat, const SparseVector<Type>& vct, Vector<Type>& result) {
			assert(mat.cols == vct.length);

			for (unsigned int row = 0; row < mat.rows; ++row) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxA = mat.row_ptr[row], idxB = 0; idxA < mat.row_ptr[row + 1] && idxB < vct.size(); ) {
					if (mat.elem[idxA].col_idx < vct.elem[idxB].idx) {
						++idxA;
					} else if (mat.elem[idxA].col_idx > vct.elem[idxB].idx) {
						++idxB;
					} else {
						prod += mat.elem[idxA].value * vct.elem[idxB].value;
						++idxA;
						++idxB;
					}
				}
				result.elem[row] = prod;
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& mat, const SparseVector<Type>& vct, SparseVector<Type>& result) {
			assert(mat.cols == vct.length);

			result.elem.clear();

			for (unsigned int row = 0; row < mat.rows; ++row) {
				Type prod = NumberTraits<Type>::zero;
				for (unsigned int idxA = mat.row_ptr[row], idxB = 0; idxA < mat.row_ptr[row + 1] && idxB < vct.size(); ) {
					if (mat.elem[idxA].col_idx < vct.elem[idxB].idx) {
						++idxA;
					} else if (mat.elem[idxA].col_idx > vct.elem[idxB].idx) {
						++idxB;
					} else {
						prod += mat.elem[idxA].value * vct.elem[idxB].value;
						++idxA;
						++idxB;
					}
				}
				result.push(row, prod);
			}
		}

		template<typename Type>
		void Operations::mat_add(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				result.elem[idxR] = op1.elem[idxR] + op2.elem[idxR];
			}
		}

		template<typename Type>
		void Operations::mat_add(Matrix<Type>& op1, const Matrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				op1.elem[idxR] += op2.elem[idxR];
			}
		}

		template<typename Type>
		void Operations::mat_add(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				result.elem[idxR] = op1.elem[idxR];
			}
			for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += result.cols) {
				for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB) {
					result.elem[idxR + op2.elem[idxB].col_idx] += op2.elem[idxB].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_add(Matrix<Type>& op1, const SparseMatrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += op1.cols) {
				for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB) {
					op1.elem[idxR + op2.elem[idxB].col_idx] += op2.elem[idxB].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_add(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op2.size(); ++idxR) {
				result.elem[idxR] = op2.elem[idxR];
			}
			for (unsigned int rowA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += result.cols) {
				for (unsigned int idxA = op1.row_ptr[rowA]; idxA < op1.row_ptr[rowA + 1]; ++idxA) {
					result.elem[idxR + op1.elem[idxA].col_idx] += op1.elem[idxA].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_add(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			result.elem.clear();

			for (unsigned int row = 0; row < op1.rows; ++row) {
				unsigned int idxA = op1.row_ptr[row];
				unsigned int idxB = op2.row_ptr[row];
				while (idxA < op1.row_ptr[row + 1] && idxB < op2.row_ptr[row + 1]) {
					if (op1.elem[idxA].col_idx < op2.elem[idxB].col_idx) {
						result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
						++idxA;
					} else if (op1.elem[idxA].col_idx > op2.elem[idxB].col_idx) {
						result.push(op2.elem[idxB].col_idx, op2.elem[idxB].value);
						++idxB;
					} else {
						result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value + op2.elem[idxB].value);
						++idxA;
						++idxB;
					}
				}
				while (idxA < op1.row_ptr[row + 1]) {
					result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
					++idxA;
				}
				while (idxB < op2.row_ptr[row + 1]) {
					result.push(op2.elem[idxB].col_idx, op2.elem[idxB].value);
					++idxB;
				}
				result.row_ptr[row + 1] = result.size();
			}
		}

		template<typename Type>
		void Operations::mat_sub(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				result.elem[idxR] = op1.elem[idxR] - op2.elem[idxR];
			}
		}

		template<typename Type>
		void Operations::mat_sub(Matrix<Type>& op1, const Matrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				op1.elem[idxR] -= op2.elem[idxR];
			}
		}

		template<typename Type>
		void Operations::mat_sub(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op1.size(); ++idxR) {
				result.elem[idxR] = op1.elem[idxR];
			}
			for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += result.cols) {
				for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB) {
					result.elem[idxR + op2.elem[idxB].col_idx] -= op2.elem[idxB].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_sub(Matrix<Type>& op1, const SparseMatrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int rowB = 0, idxR = 0; rowB < op2.rows; ++rowB, idxR += op1.cols) {
				for (unsigned int idxB = op2.row_ptr[rowB]; idxB < op2.row_ptr[rowB + 1]; ++idxB) {
					op1.elem[idxR + op2.elem[idxB].col_idx] -= op2.elem[idxB].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_sub(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			for (unsigned int idxR = 0; idxR < op2.size(); ++idxR) {
				result.elem[idxR] = -op2.elem[idxR];
			}
			for (unsigned int rowA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += result.cols) {
				for (unsigned int idxA = op1.row_ptr[rowA]; idxA < op1.row_ptr[rowA + 1]; ++idxA) {
					result.elem[idxR + op1.elem[idxA].col_idx] += op1.elem[idxA].value;
				}
			}
		}

		template<typename Type>
		void Operations::mat_sub(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			result.elem.clear();

			for (unsigned int row = 0; row < op1.rows; ++row) {
				unsigned int idxA = op1.row_ptr[row];
				unsigned int idxB = op2.row_ptr[row];
				while (idxA < op1.row_ptr[row + 1] && idxB < op2.row_ptr[row + 1]) {
					if (op1.elem[idxA].col_idx < op2.elem[idxB].col_idx) {
						result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
						++idxA;
					} else if (op1.elem[idxA].col_idx > op2.elem[idxB].col_idx) {
						result.push(op2.elem[idxB].col_idx, -op2.elem[idxB].value);
						++idxB;
					} else {
						result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value - op2.elem[idxB].value);
						++idxA;
						++idxB;
					}
				}
				while (idxA < op1.row_ptr[row + 1]) {
					result.push(op1.elem[idxA].col_idx, op1.elem[idxA].value);
					++idxA;
				}
				while (idxB < op2.row_ptr[row + 1]) {
					result.push(op2.elem[idxB].col_idx, -op2.elem[idxB].value);
					++idxB;
				}
				result.row_ptr[row + 1] = result.size();
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.cols == op2.rows);

			// TODO: Improve indexing
			for (unsigned int rowR = 0, idxR = 0; rowR < op1.rows; ++rowR) {
				for (unsigned int colR = 0; colR < op2.cols; ++colR, ++idxR) {
					Type prod = NumberTraits<Type>::zero;
					for (unsigned int idx = 0, idxA = rowR * op1.cols, idxB = colR; idx < op1.cols; ++idx, ++idxA, idxB += op2.cols) {
						prod += op1.elem[idxA] * op2.elem[idxB];
					}
					result.elem[idxR] = prod;
				}
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& op1, const Matrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.cols == op2.rows);

			result.elem.clear();

			for (unsigned int rowR = 0; op1.rowR < op1.rows; ++rowR) {
				for (unsigned int colR = 0; colR < op2.cols; ++colR) {
					Type prod = NumberTraits<Type>::zero;
					for (unsigned int idx = 0, idxA = rowR * op1.cols, idxB = colR; idx < op1.cols; ++idx, ++idxA, idxB += op2.cols) {
						prod += op1.elem[idxA] * op2.elem[idxB];
					}
					result.push(colR, prod);
				}
				result.row_ptr[rowR + 1] = result.size();
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.cols == op2.rows);

			std::fill(result.begin(), result.end(), NumberTraits<Type>::zero);

			for (unsigned int rowA = 0, idxA = 0, idxR = 0; rowA < op1.rows; ++rowA, idxR += op2.cols) {
				for (unsigned int colA = 0; colA < op1.cols; ++colA, ++idxA) {
					for (unsigned int idxB = op2.row_ptr[colA]; idxB < op2.row_ptr[colA + 1]; ++idxB) {
						result.elem[idxR + op2.elem[idxB].col_idx] += op1.elem[idxA] * op2.elem[idxB].value;
					}
				}
			}
		}

		template<typename Type>
		void Operations::mat_mul(const Matrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.cols == op2.rows);

			result.elem.clear();

			std::vector<Type> cache(op2.cols);
			for (unsigned int rowA = 0, idxA = 0; rowA < op1.rows; ++rowA) {
				std::fill(cache.begin(), cache.end(), NumberTraits<Type>::zero);
				for (unsigned int colA = 0; colA < op1.cols; ++colA, ++idxA) {
					for (unsigned int idxB = op2.row_ptr[colA]; idxB < op2.row_ptr[colA + 1]; ++idxB) {
						cache[op2.elem[idxB].col_idx] += op1.elem[idxA] * op2.elem[idxB].value;
					}
				}
				for (unsigned int idxC = 0; idxC < op2.cols; ++idxC) {
					result.push(idxC, cache[idxC]);
				}
				result.row_ptr[rowA + 1] = result.size();
			}

		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.cols == op2.rows);
			
			std::fill(result.elem.begin(), result.elem.end(), NumberTraits<Type>::zero);

			for (unsigned int row = 0, idxR = 0; row < op1.rows; ++row, idxR += op2.cols) {
				for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA) {
					for (unsigned int idxRR = idxR, idxB = op1.elem[idxA].col_idx * op2.cols; idxRR < idxR + op2.cols; ++idxRR, ++idxB) {
						result.elem[idxRR] += op1.elem[idxA].value * op2.elem[idxB];
					}
				}
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& op1, const Matrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.cols == op2.rows);

			result.elem.clear();

			std::vector<Type> cache(op2.cols);
			for (unsigned int row = 0; row < op1.rows; ++row) {
				std::fill(cache.begin(), cache.end(), NumberTraits<Type>::zero);
				for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA) {
					unsigned int idxB = op1.elem[idxA].col_idx * op2.cols;
					for (unsigned int idxC = 0; idxC < op2.cols; ++idxC, ++idxB) {
						cache[idxC] += op1.elem[idxA].value * op2.elem[idxB];
					}
				}
				for (unsigned int idxC = 0; idxC < op2.cols; ++idxC) {
					result.push(idxC, cache[idxC]);
				}
				result.row_ptr[row + 1] = result.size();
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, Matrix<Type>& result) {
			assert(op1.cols == op2.rows);

			std::fill(result.elem.begin(), result.elem.end(), NumberTraits<Type>::zero);

			for (unsigned int row = 0, idxR = 0; row < op1.rows; ++row, idxR += op2.cols) {
				for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA) {
					for (unsigned int idxB = op2.row_ptr[op1.elem[idxA].col_idx]; idxB < op2.row_ptr[op1.elem[idxA].col_idx + 1]; ++idxB) {
						result.elem[idxR + op2.elem[idxB].col_idx] += op1.elem[idxA].value * op2.elem[idxB].value;
					}
				}
			}
		}

		template<typename Type>
		void Operations::mat_mul(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2, SparseMatrix<Type>& result) {
			assert(op1.cols == op2.rows);

			std::vector<Type> cache(op2.cols);
			for (unsigned int row = 0; row < op1.rows; ++row) {
				std::fill(cache.begin(), cache.end(), NumberTraits<Type>::zero);
				for (unsigned int idxA = op1.row_ptr[row]; idxA < op1.row_ptr[row + 1]; ++idxA) {
					for (unsigned int idxB = op2.row_ptr[op1.elem[idxA].col_idx]; idxB < op2.row_ptr[op1.elem[idxA].col_idx + 1]; ++idxB) {
						cache[op2.elem[idxB].col_idx] += op1.elem[idxA].value * op2.elem[idxB].value;
					}
				}
				for (unsigned int idxC = 0; idxC < op2.cols; ++idxC) {
					result.push(idxC, cache[idxC]);
				}
				result.row_ptr[row + 1] = result.size();
			}
		}

		template<typename Type>
		bool Operations::mat_eq(const Matrix<Type>& op1, const Matrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			bool l = true;
			for (unsigned int i = 0; i < op1.size() && l; ++i) {
				l = (op1.elem[i] == op2.elem[i]);
			}

			return l;
		}

		template<typename Type>
		bool Operations::mat_eq(const Matrix<Type>& op1, const SparseMatrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			bool l = true;
			for (unsigned int row = 0, idxA = 0; op1.row < op1.rows && l; ++row) {
				unsigned int colA = 0;
				unsigned int idxB = op2.row_ptr[row];
				while (colA < op1.cols && idxB < op2.row_ptr[row + 1] && l) {
					if (colA == op2.elem[idxB].col_idx) {
						l = (op1.elem[idxA] == op2.elem[idxB].value);
						++colA;
						++idxA;
						++idxB;
					} else {
						l = (op1.elem[idxA] == NumberTraits<Type>::zero);
						++colA;
						++idxA;
					}
				}
				while (colA < op1.cols && l) {
					l = (op1.elem[idxA] == NumberTraits<Type>::zero);
					++colA;
					++idxA;
				}
			}

			return l;
		}

		template<typename Type>
		bool Operations::mat_eq(const SparseMatrix<Type>& op1, const Matrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			bool l = true;
			for (unsigned int row = 0, idxB = 0; row < op1.rows && l; ++row) {
				unsigned int idxA = op2.row_ptr[row];
				unsigned int colB = 0;
				while (idxA < op1.row_ptr[row + 1] && colB < op2.cols && l) {
					if (op1.elem[idxA].col_idx == colB) {
						l = (op1.elem[idxA].value == op2.elem[idxB]);
						++idxA;
						++colB;
						++idxB;
					} else {
						l = (NumberTraits<Type>::zero == op2.elem[idxB]);
						++colB;
						++idxB;
					}
				}
				while (colB < op2.cols && l) {
					l = (NumberTraits<Type>::zero == op2.elem[idxB]);
					++colB;
					++idxB;
				}
			}

			return l;
		}

		template<typename Type>
		bool Operations::mat_eq(const SparseMatrix<Type>& op1, const SparseMatrix<Type>& op2) {
			assert(op1.rows == op2.rows && op1.cols == op2.cols);

			bool l = (op1.row_ptr[0] == op2.row_ptr[0]);
			for (unsigned int row = 0; row < op1.rows && l; ++row) {
				l = (op1.row_ptr[row + 1] == op2.row_ptr[row + 1]);
				for (unsigned int idx = op1.row_ptr[row]; idx < op1.row_ptr[row + 1] && l; ++idx) {
					l = ((op1.elem[idx].col_idx == op2.elem[idx].col_idx) && (op1.elem[idx].value == op2.elem[idx].value));
				}
			}

			return l;
		}

	}
}
