#ifndef LZ_DEF_LZ_H_basic_src_basic_std
#define LZ_DEF_LZ_H_basic_src_basic_std 202105L

//C++11 or newer:
#if __cplusplus<201103L
	#error ISO C++ 2011 standard is required.
#endif

//Other C++ standard libraries:
//language feature:
#include<cstddef>
#include<cstdlib>
#include<memory>
//traits:
#include<type_traits>
#include<limits>
//iterator and container:
#include<iterator>
#include<vector>
//math:
#include<cmath>
#include<complex>
#include<algorithm>
#include<numeric>
#include<functional>
#include<utility>

//extending <type_traits>:
namespace std{
	template<typename _Tp>
	struct is_unsigned_integral:
	public is_integral<typename conditional<is_unsigned<_Tp>::value,_Tp,void>::type>
	{};
} //namespace std;

#endif //#ifndef LZ_DEF_LZ_H_basic_src_basic_std