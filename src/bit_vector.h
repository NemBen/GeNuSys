#ifndef GENUSYS_UTILS_BIT_VECTOR_H_
#define GENUSYS_UTILS_BIT_VECTOR_H_

#include <vector>

namespace GeNuSys {
	namespace Utils {

		class BitVector {

		private:

			std::vector<unsigned long long> table;

		public:

			BitVector(unsigned long long size);

			bool operator [](unsigned long long idx) const;

			void set(unsigned long long idx);

			void unset(unsigned long long idx);

		};

	}
}

#endif // GENUSYS_UTILS_BIT_VECTOR_H_
