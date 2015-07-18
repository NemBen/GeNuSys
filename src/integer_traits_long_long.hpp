namespace GeNuSys {

	template<>
	bool IntegerTraits<long long>::divisible(const long long& a, const long long& b) {
		return a % b == 0;
	}
	
	template<>
	long long IntegerTraits<long long>::div(const long long& a, const long long& b) {
		return a / b;
	}
	
	template<>
	long long IntegerTraits<long long>::mod(const long long& a, const long long& b) {
		return (a % b + b) % b;
	}
	
	template<>
	long long IntegerTraits<long long>::mods(const long long& a, const long long& b) {
		long long m = mod(a, b);
		return (m > b / 2) ? m - b : m;
	}

}
