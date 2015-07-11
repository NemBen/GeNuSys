namespace GeNuSys {

	template<>
	bool IntegerTraits<int>::divisible(const int& a, const int& b) {
		return a % b == 0;
	}

	template<>
	int IntegerTraits<int>::idiv(const int& a, const int& b) {
		return a / b;
	}

	template<>
	int IntegerTraits<int>::mod(const int& a, const int& b) {
		return (a % b + b) % b;
	}

	template<>
	int IntegerTraits<int>::mods(const int& a, const int& b) {
		int m = mod(a, b);
		return (m > b / 2) ?  m - b : m;
	}

}
