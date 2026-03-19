template<int fun_id=-1>
class Fun_main{
	public:
	int operator ()() const {return 0;}
};
#define LZ_DEF_compile_part_val 2
#define LZ_DEF_compile_part(x) LZ_DEF_compile_part_val==x
#define LZ_DEF_compile_part_setfun(id,fun) \
	template<>class Fun_main<id>{public:int operator()()const{return fun();}}
#define LZ_DEF_compile_part_setfun_auto(fun) LZ_DEF_compile_part_setfun(LZ_DEF_compile_part_val,fun)
#define LZ_DEF_run_mainfun(x) Fun_main<x> fun_main;return fun_main()

#define DEF_debug_pos std::cout<<__FILE__<<"("<<__LINE__<<")"<<std::endl;
#define DEF_debug_out_head std::cout<<__FILE__<<"("<<__LINE__<<"): "
#define DEF_debug_out(x) DEF_debug_out_head<<(x)<<std::endl

#include<vector>
#include<iostream>
using namespace std;
//`#include"stats/test_DOE_crit.hpp"
//`#include"stats/test_DOE_crit_aberration.hpp"
//`#include"math/test_function.hpp"
//`#include"math/test_differentiation.hpp"
//`#include"math/test_minimisation.hpp"
#include"math/test_prime.hpp"
//`#include"math/test_polynomial.hpp"
//`#include"math/test_GaloisField.hpp"

int main(){
	#if false
	std::vector<int(*)(void)> fun_main={
		/*0*/&fun_test_Discrepancy_Leb_infty,&fun_test_Iterator_design_canon
		};
	return (*(fun_main[1]))();
	#endif
	LZ_DEF_run_mainfun(LZ_DEF_compile_part_val);
}