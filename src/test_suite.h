#ifndef GENUSYS_TESTS_TEST_SUITE_H_
#define GENUSYS_TESTS_TEST_SUITE_H_

#include <string>

namespace GeNuSys {
	namespace Tests {

		class TestSuite {

			friend class TestRunner;

		private:

			std::string name;

			int cnt;

			int passed;

			int failed;

			void fail(const std::string& message);

			void pass(const std::string& message);

		public:

			TestSuite(const std::string& name);

			virtual ~TestSuite();

			void runTestSuite();

			virtual void run() = 0;

			void assertTrue(bool value, const std::string& message = "");

			void assertFalse(bool value, const std::string& message = "");

			template<typename T>
			void assertEqual(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertNotEqual(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertLess(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertLessOrEqual(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertMore(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertMoreOrEqual(const T& expected, const T& value, const std::string& message = "");

			template<typename T>
			void assertBetween(const T& low, const T& high, const T& value, const std::string& message = "");

			template<typename T>
			void assertNotBetween(const T& low, const T& high, const T& value, const std::string& message = "");

			template<typename T>
			void assertDifference(const T& expected, const T& threshold, const T& value, const std::string& message = "");

		};

	}
}

#include "test_suite.hpp"

#endif // GENUSYS_TESTS_TEST_SUITE_H_
