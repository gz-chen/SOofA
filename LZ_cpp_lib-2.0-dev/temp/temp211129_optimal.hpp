#ifndef LZ_DEF_LZ_H_math_optimal
#define LZ_DEF_LZ_H_math_optimal 202105L

#include"LZ_H/math/distribution.hpp"

namespace Liuze{
namespace math{
/****************************************************************************************************/
//ThresholdAccepting{
//fun: minimiser_ThresholdAccepting:
template<typename Tp,typename TargetType,typename FuncTranType,typename ValType=Type_Real,
typename SizeType=Type_Size,
typename RandUIntType=typename Type_RandEngine::result_type>
Tp minimiser_ThresholdAccepting
(TargetType & func_tar,Tp const & result_init,FuncTranType & func_tran,
RandomEngine<RandUIntType> & engine,
Type_ArrayTemp<ValType> const & threshold,SizeType const & num_iter_inner){
	typedef Tp type_res; //the result;
	typedef RandomEngine<RandUIntType> type_rand_eng; //the random engine;
	typedef SizeType type_size; //the integral;
	typedef ValType type_val; //the output value of the target function;
	//checks: types:
	LZ_DEF_func_check_traits((std::is_integral<RandUIntType>::value));
	LZ_DEF_func_check_traits((std::is_integral<type_size>::value));
	LZ_DEF_func_check_traits((std::is_arithmetic<type_val>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::result_of<typename std::decay<TargetType>::type(type_res const &)>::type,
		type_val>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::result_of<typename std::decay<FuncTranType>::type(
		type_res const &,type_rand_eng &)>::type,
		type_res>::value));
	//checks: values:
	if(num_iter_inner<=(type_size)0 || threshold.size()<=(Type_Size)0) return result_init;
	//checks end.
	Type_Size i_out;
	type_size i_in;
	type_res res0=result_init,res1=result_init;
	type_val val0=func_tar(res0),val1=val0;
	for(i_out=(Type_Size)0;i_out<threshold.size();++i_out){
		for(i_in=(type_size)0;i_in<num_iter_inner;++i_in){
			res1=func_tran(res0,engine);
			val1=func_tar(res1);
			if(val1-val0<=threshold[i_out]){
				val0=val1;
				res0=res1;
			}
		}
	}
	return res0;
} //(fun: minimiser_ThresholdAccepting)

//fun: minimiser_helper_ThresholdAccepting_threshold_adaptive:
//(See Lin, Sharpe and Winker (2010). Optimized U-type designs on flexible regions. Section 4)
template<typename Tp,typename TargetType,typename FuncTranType,typename ValType=Type_Real,
typename RandUIntType=typename Type_RandEngine::result_type>
Type_ArrayTemp<ValType> minimiser_helper_ThresholdAccepting_threshold_adaptive
(TargetType & func_tar,Tp const & result_init,FuncTranType & func_tran,
RandomEngine<RandUIntType> & engine,Type_Size const & num_threshold){
	typedef Tp type_res; //the result;
	typedef RandomEngine<RandUIntType> type_rand_eng; //the random engine;
	typedef ValType type_val; //the output value of the target function;
	//checks: types:
	LZ_DEF_func_check_traits((std::is_integral<RandUIntType>::value));
	LZ_DEF_func_check_traits((std::is_arithmetic<type_val>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::result_of<typename std::decay<TargetType>::type(type_res const &)>::type,
		type_val>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::result_of<typename std::decay<FuncTranType>::type(
		type_res const &,type_rand_eng &)>::type,
		type_res>::value));
	//checks: values:
	if(num_threshold<=(type_size)0) return Type_ArrayTemp<type_val>((Type_Size)0);
	//checks end.
	Type_Size i;
	type_res res1=result_init;
	type_val val0=func_tar(result_init);
	Type_ArrayTemp<type_val> threshold(num_threshold);
	for(i=(Type_Size)0;i<threshold.size();++i){
		res1=func_tran(result_init,engine);
		threshold[i]=abs(func_tar(res1)-val0);
	}
	auto order_desc=[](type_val const & y0,type_val const & y1)->Type_Bool{return y0>y1;};
	std::sort(threshold.begin(),threshold.end(),order_desc);
	return threshold;
} //(fun: minimiser_helper_ThresholdAccepting_threshold_adaptive)
//}(ThresholdAccepting)
/****************************************************************************************************/
} //namespace Liuze;
} //namespace math;

#endif //#ifndef LZ_DEF_LZ_H_math_optimal