#ifndef LZ_DEF_LZ_H_basic_basic
#define LZ_DEF_LZ_H_basic_basic 202105L

//external library files:
#define LZ_DEF_extLIB_Eigen 30309L //version 3.3.9;

//including files:
#include"LZ_H/basic/src/basic_std.hpp"
#ifdef LZ_DEF_extLIB_Eigen
	#include"LZ_H/basic/src/basic_Eigen.hpp"
#endif

//#define: deleted codes:
#define LZ_DEF_compile_ignore false

//checking traits:
#define LZ_DEF_func_check_traits(boolexpr) static_assert(boolexpr,"Error in LZ_DEF_func_check_traits.")

//typedef: default types:
namespace Liuze{
	typedef bool Type_Bool; //bool (true or false);
	typedef unsigned char Type_Stat; //state;
	typedef std::size_t Type_Size; //number of bytes of objects;
	typedef unsigned int Type_UInt; //unsigned integral;
	typedef int Type_Int; //signed integral;
	typedef double Type_Real; //floating point;
	template<typename RealType>
	using Type_CompTemp=std::complex<RealType>; //template: complex number;
	typedef Type_CompTemp<Type_Real> Type_Comp; //complex number;
	typedef char Type_Char; //character;
	template<typename Tp>
	using Type_ArrayTemp=std::vector<Tp>; //template: array, sequence and so on;
	#ifdef LZ_DEF_extLIB_Eigen
		template<typename NumType>
		using Type_MatTemp=Eigen::MatrixX<NumType>; //template: matrix;
	#endif //#ifdef LZ_DEF_extLIB_Eigen
} //namespace Liuze;

//#define: LZ_DEF_const_default_...
#ifndef LZ_DEF_const_default_precision
	#define LZ_DEF_const_default_precision 1e-5
#endif

#include"LZ_H/basic/src/basic_utility.hpp"

#endif //#ifndef LZ_DEF_LZ_H_basic_basic