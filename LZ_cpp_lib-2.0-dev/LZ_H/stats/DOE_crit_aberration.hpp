#ifndef LZ_DEF_LZ_H_stats_DOE_crit_aberration
#define LZ_DEF_LZ_H_stats_DOE_crit_aberration 202105L

#include<array>
#include<map>
#include<Eigen/QR>
#include"LZ_H/stats/basic.hpp"

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
	protected:
		//types:
		class type_tab_polyval;
		//the type for calculating and storing the values of the orthogonal polynoimal contrasts;
		class type_tab_haarval;
		//the type for calculating and storing the values of the Haar-like contrasts;
		
	public:
		//Constructor:
		WordlengthPattern(void)=default;
		WordlengthPattern(type_this const & crit)=delete;
		WordlengthPattern(type_this && crit)=delete;
		
		//Destructor:
		~WordlengthPattern(void)=default;
		
		//operator:
		type_this & operator =(type_this const & crit)=delete;
		type_this & operator =(type_this && crit)=delete;
		
		//other function:
		//Alphabeta: for qualitative-quantitative mixed-level unequally-weighted designs.
		//Each design point is a column vector with integer entries,
		//whose upper and lower parts are respectively
		//the settings for qualitative and quantitative factors.
		//Weights should be positive real numbers and are standardised automatically.
		//`num_level_qual' and `num_level_quan' respectively list
		//the numbers of levels of qualitative and quantitative factors.
		//The numbers of levels should be integers no less than 1.
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		type_res Alphabeta
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		Type_ArrayTemp<type_int> const & num_level_qual,
		Type_ArrayTemp<type_int> const & num_level_quan) const;
		//enumerator_Alphabeta:
		//See Tang, Yu and Xu, Hongquan (2019). Wordlength Enumerator for Fractional Factorial Designs.
		//This is a modified version for
		//qualitative-quantitative mixed-level unequally-weighted designs.
		//Each design point is a column vector with integer entries,
		//whose upper and lower parts are respectively
		//the settings for qualitative and quantitative factors.
		//Weights should be positive real numbers and are standardised automatically.
		//`num_level_qual' and `num_level_quan' respectively list
		//the numbers of levels of qualitative and quantitative factors.
		//The numbers of levels should be integers no less than 1.
		//The enumerator is evaluated at `val'.
		//`type_real' should be able to be converted to `ValType'.
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename ValType,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		ValType enumerator_Alphabeta
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		Type_ArrayTemp<type_int> const & num_level_qual,
		Type_ArrayTemp<type_int> const & num_level_quan,
		ValType const & val) const;
		
		//static function:
		static typename type_tab_polyval::type_mat_real
		get_contrast(type_int const & num_level);
		static typename type_tab_haarval::type_mat_real
		get_contrast_Haar_value(type_int const & num_level,type_int const & base);
		static typename type_tab_haarval::type_arr_int
		get_contrast_Haar_order(type_int const & num_level,type_int const & base);
		
	protected:
		//static member values:
		static type_tab_polyval m_tab_polyval;
		//it stores the values of the orthogonal polynoimal contrasts;
		static type_tab_haarval m_tab_haarval;
		//it stores the values of the Haar-like contrasts;
}; //class WordlengthPattern;
/****************************************************************************************************/
template<typename RealType=Type_Real,typename IntType=Type_Int,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value &&
	std::is_integral<IntType>::value>::type>
class StratificationPattern{
	public:
		//types:
		typedef StratificationPattern<RealType,IntType> type_this;
		typedef RealType type_real;
		typedef IntType type_int;
		typedef Type_MatTemp<type_int> type_des;
		typedef Type_ArrayTemp<type_real> type_res;
	protected:
		//types:
		class type_tab_haarval;
		//the type for calculating and storing the values of the Haar-like contrasts;
		
	public:
		//Constructor:
		StratificationPattern(void)=default;
		StratificationPattern(type_this const & crit)=delete;
		StratificationPattern(type_this && crit)=delete;
		
		//Destructor:
		~StratificationPattern(void)=default;
		
		//operator:
		type_this & operator =(type_this const & crit)=delete;
		type_this & operator =(type_this && crit)=delete;
		
		//other function:
		//Stratification: for quantitative mixed-level unequally-weighted designs.
		//Each design point is a column vector with integer entries.
		//Weights should be positive real numbers and are standardised automatically.
		//`num_level' lists the numbers of levels of the quantitative factors.
		//The numbers of levels should be integers no less than 1.
		//`base' lists the base numbers for assessing the stratification property.
		//The base numbers should be integers no less than 2.
		//`order_max' is the maximum order of the pattern to be calculated.
		//When `order_max' is 0, the full pattern is calculated,
		//and when `order_max' is positive,
		//the pattern of orders up to `order_max' is calculated
		//with an error proportional to `err'.
		//For `err', -1 means the default setting.
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		type_res Stratification
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base) const;
		//Stratification: for quantitative mixed-level unequally-weighted designs.
		//This function is basically the same as the former one,
		//but has more arguments.
		//`order_max' is the maximum order of the pattern to be calculated.
		//When `order_max' is 0, the full pattern is calculated,
		//and when `order_max' is positive,
		//the pattern of orders up to `order_max' is calculated using an approximate algorithm
		//with a relative error no more than `err'.
		//For `err', -1 means the default setting.
		//For the approximate algorithm,
		//see Theorem 4 of
		//(Tian, Ye and Xu, Hongquan (2023). Stratification Pattern Enumerator and its Applications).
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		type_res Stratification
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base,
		type_int order_max,type_real err=type_real(-1)) const;
		//enumerator_Stratification:
		//See Tian, Ye and Xu, Hongquan (2023). Stratification Pattern Enumerator and its Applications.
		//This is a modified version for
		//quantitative mixed-level unequally-weighted designs.
		//Each design point is a column vector with integer entries.
		//Weights should be positive real numbers and are standardised automatically.
		//`num_level' lists the numbers of levels of the quantitative factors.
		//The numbers of levels should be integers no less than 1.
		//`base' lists the base numbers for assessing the stratification property.
		//The base numbers should be integers no less than 2.
		//The enumerator is evaluated at `val'.
		//`type_real' should be able to be converted to `ValType'.
		template<typename ForwardIterType_des,typename ForwardIterType_weight,
			typename ValType,
			typename=typename std::enable_if<
				std::is_convertible<typename std::iterator_traits<ForwardIterType_des>::value_type,
				type_des>::value &&
				std::is_convertible<typename std::iterator_traits<ForwardIterType_weight>::value_type,
				type_real>::value
			>::type>
		ValType enumerator_Stratification
		(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
		Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base,
		ValType const & val) const;
		
		//static function:
		static typename type_tab_haarval::type_mat_real
		get_contrast_Haar_value(type_int const & num_level,type_int const & base);
		static typename type_tab_haarval::type_arr_int
		get_contrast_Haar_order(type_int const & num_level,type_int const & base);
		
	protected:
		//static member values:
		static type_tab_haarval m_tab_haarval;
		//it stores the values of the Haar-like contrasts;
}; //class StratificationPattern;
/****************************************************************************************************/

} //namespace crit;
} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/DOE/DOE_crit_aberration.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_stats_DOE_crit_aberration