#include <iostream>

namespace GeNuSys {
	namespace Tests {

		template<typename T>
		void TestSuite::assertEqual(const T& expected, const T& value, const std::string& message) {
			if (expected == value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertNotEqual(const T& expected, const T& value, const std::string& message) {
			if (expected != value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: NOT " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertLess(const T& expected, const T& value, const std::string& message) {
			if (expected > value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: LESS THAN " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertLessOrEqual(const T& expected, const T& value, const std::string& message) {
			if (expected >= value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: LESS THAN OR EQUAL TO " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertMore(const T& expected, const T& value, const std::string& message) {
			if (expected < value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: MORE THAN " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertMoreOrEqual(const T& expected, const T& value, const std::string& message) {
			if (expected <= value) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: MORE THAN OR EQUAL TO " << expected << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertBetween(const T& low, const T& high, const T& value, const std::string& message) {
			if (low <= value && value <= high) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: BETWEEN " << low << " AND " << high << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertNotBetween(const T& low, const T& high, const T& value, const std::string& message) {
			if (!(low <= value && value <= high)) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: NOT BETWEEN " << low << " AND " << high << " - Actual value: " << value << std::endl;
			}
		}

		template<typename T>
		void TestSuite::assertDifference(const T& expected, const T& threshold, const T& value, const std::string& message) {
			T diff = expected - value;
			if (-threshold <= diff && diff <= threshold) {
				pass(message);
			} else {
				fail(message);
				std::cout << "Expected: CLOSER TO " << expected << " THAN " << threshold << " - Actual value: " << value << " (DIFFERENCE: " << diff << ")" << std::endl;
			}
		}

	}
}
