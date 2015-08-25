#include "bit_vector.h"

namespace GeNuSys {
	namespace Utils {
		
		BitVector::BitVector(unsigned long long size):table((unsigned int)(size / (sizeof(unsigned long long) * 8)) + 1, 0) {
		}

		bool BitVector::operator [](unsigned long long idx) const {
			unsigned int d = (int) (idx / (sizeof(unsigned long long) * 8));
			unsigned int m = (int) (idx % (sizeof(unsigned long long) * 8));

			return (table[d] & ((unsigned long long) 1 << m)) != 0;
		}

		void BitVector::set(unsigned long long idx) {
			unsigned int d = (int) (idx / (sizeof(unsigned long long) * 8));
			unsigned int m = (int) (idx % (sizeof(unsigned long long) * 8));

			table[d] |= ((unsigned long long) 1 << m);
		}

		void BitVector::unset(unsigned long long idx) {
			unsigned int d = (int) (idx / (sizeof(unsigned long long) * 8));
			unsigned int m = (int) (idx % (sizeof(unsigned long long) * 8));

			table[d] &= ~((unsigned long long) 1 << m);
		}

	}
}
