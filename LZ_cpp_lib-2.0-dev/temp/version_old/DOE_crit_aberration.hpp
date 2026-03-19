#ifndef LZ_DEF_LZ_H_stats_DOE_crit_aberration
#define LZ_DEF_LZ_H_stats_DOE_crit_aberration 202105L

#include<Eigen/QR>
#include"LZ_H/stats/basic.hpp"
//`#include"LZ_H/math/polynomial.hpp"

#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace stats{
namespace crit{ //a namespace for criteria;

/****************************************************************************************************/
template<typename RealType=Type_Real,typename IntType=Type_Int,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value &&
	std::is_integral<IntType>::value>::type>
class WordlengthPattern{
	public:
		//types:
		typedef WordlengthPattern<RealType,IntType> type_this;
		typedef RealType type_real;
		typedef IntType type_int;
		typedef Type_MatTemp<type_int> type_des;
		typedef Type_ArrayTemp<type_real> type_res;
		
		//Constructor:
		explicit WordlengthPattern(type_int const & num_level);
		WordlengthPattern(type_this const & crit)=delete;
		WordlengthPattern(type_this && crit)=delete;
		
		//Destructor:
		~WordlengthPattern(void)=default;
		
		//other function:
		//beta: for beta-wordlength pattern:
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		type_res beta
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const;
		//enumerator:
		//See Tang, Yu and Xu, Hongquan (2019). Wordlength Enumerator for Fractional Factorial Designs.
		//The input 'enumval_begin' is the beginning iterator for {y_1,...,y_{'m_num_level'-1}},
		//and y_0 is always 1.
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename ForwardIterType_enumval,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		typename std::decay<typename std::iterator_traits<ForwardIterType_enumval>::value_type>::type
		enumerator
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		ForwardIterType_enumval const & enumval_begin) const;
	protected:
		//types:
		class type_tab_powval;
		
		//member values:
		type_int m_num_level;
		Type_ArrayTemp<Type_ArrayTemp<type_real> > m_tab_poly;
		Type_ArrayTemp<Type_ArrayTemp<type_real> > m_tab_polyval;
		
		//static member values:
		static type_tab_powval m_tab_powval;
}; //class WordlengthPattern;
/****************************************************************************************************/

} //namespace crit;
} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/DOE/DOE_crit_aberration.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_stats_DOE_crit_aberration