#ifndef LZ_DEF_LZ_H_math_prime
#define LZ_DEF_LZ_H_math_prime 202105L

#include"LZ_H/math/basic.hpp"

namespace Liuze{
namespace math{
/****************************************************************************************************/
//fun: Is_prime:
//returns `true' if `number' is a prime integer:
template<typename IntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Type_Bool Is_prime(IntType number);

//fun: Get_prime_all:
//returns all prime numbers in {2,...,`upper'}:
template<typename IntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Type_ArrayTemp<IntType> Get_prime_all
(IntType const & upper);

//fun: Get_prime_all:
//returns all prime numbers in {2,...,`upper'}, based on prime numbers no more than `prior_upper':
template<typename IntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Type_ArrayTemp<IntType> Get_prime_all
(IntType const & upper,Type_ArrayTemp<IntType> const & prior_prime,IntType prior_upper=0);

//fun: Get_prime_atleast:
//Iterators `least_begin' and `least_end' defines a sequence of natural numbers,
//say (a_{0},...,a_{n-1});
//this function finds the smallest prime number p_{i} that is no less than a_{i} for each i in n,
//and puts these primes numbers into the iterator starting at `prime_begin';
//the return value is the list of all prime numbers that is no more than p_{i*},
//where a_{i*} is the largest of (a_{0},...,a_{n-1});
//users can use `prime_list' and `prime_list_upper' to
//tell this function the increasing sequence of all primes that are no more than `prime_list_upper':
template<typename ForwardIterType,typename OutputIterType,
	typename IntType=typename std::remove_cv<typename std::remove_reference<
		typename ForwardIterType::value_type>::type>::type,
	typename=typename std::enable_if<
		std::is_integral<IntType>::value &&
		std::is_convertible<typename ForwardIterType::value_type,IntType>::value &&
		std::is_convertible<IntType,typename OutputIterType::value_type>::value
	>::type>
Type_ArrayTemp<IntType> Get_prime_atleast
(ForwardIterType const & least_begin,ForwardIterType const & least_end,
OutputIterType prime_begin,
Type_ArrayTemp<IntType> const prime_list=Type_ArrayTemp<IntType>((Type_Size)0),
IntType prime_list_upper=0);

//fun: integer_prime_decomposition:
//number: an integer;
//returns `decomp' such that,
//when the absolute value of `number' is at least 2,
//abs(`number') = `decomp[0][0]'^`decomp[0][1]' * ... , and `decomp[0][0]', ... are prime numbers,
//and otherwise, the size of `decomp' is 0.
template<typename IntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Type_ArrayTemp<Type_ArrayTemp<IntType> > integer_prime_decomposition(IntType const & number);
/****************************************************************************************************/
//class: PrimePowerQuerier:
template<typename IntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
class PrimePowerQuerier{
	public:
		//type:
		typedef PrimePowerQuerier<IntType> type_this;
		typedef IntType type_int;
		typedef Type_ArrayTemp<type_int> type_subtab; //subtable of numbers;
		typedef Type_ArrayTemp<type_subtab> type_tab; //table of numbers;
		
	public:
		//constructor:
		explicit PrimePowerQuerier(type_int upper=type_int(0));
		PrimePowerQuerier(type_this const & querier);
		PrimePowerQuerier(type_this && querier);
		
		//destructor:
		~PrimePowerQuerier(void)=default;
		
		//operator:
		type_this & operator =(type_this const & querier);
		type_this & operator =(type_this && querier);
		
		//other function:
		type_int const & get_upper() const;
		type_tab const & get_table_ascending() const;
		type_tab const & get_table_base() const;
		type_int set_upper(type_int const & upper);
	protected:
		//other function:
		void Cal_init();
		
	protected:
		//member value:
		type_int m_upper; //this querier stores prime powers in {2,...,`upper'};
		type_tab m_tab_asc;
			//this table stores prime powers in ascending order;
			//each subtable stores {value,base,exponent};
		type_tab m_tab_bas;
			//this table stores prime powers classified by the prime base number;
			//each subtable stores ascending powers of a prime base number;
}; //class: PrimePowerQuerier;
//(class: PrimePowerQuerier)
/****************************************************************************************************/
} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/prime.cpp"

#endif //#ifndef LZ_DEF_LZ_H_math_prime