#ifndef LZ_DEF_LZ_H_math_function
#define LZ_DEF_LZ_H_math_function 202105L

#include"LZ_H/math/basic.hpp"

namespace Liuze{
namespace math{

/****************************************************************************************************/
//general arithmetic{

//Get zero:
template<typename ValType>
ValType zero();
template<typename ValType>
ValType zero(ValType const &);

//Get one:
template<typename ValType>
ValType one();
template<typename ValType>
ValType one(ValType const &);

//Check whether an object is (approximately) zero:
template<typename ValType>
Type_Bool Is_approx_zero(ValType const & x);

//Inverse for multiplication:
template<typename ValType>
ValType inverse(ValType const & x);

//Take integer power:
template<typename ValType,typename IntType=Type_Int>
ValType powint(ValType const & base,IntType const & exp);

//Division with Remainder:
template<typename FactorType>
Type_Bool div
(FactorType const & numer,FactorType const & denom,
FactorType * ptr_quot=NULL,FactorType * ptr_rem=NULL);

//Division algorithm (Euclidean algorithm):
//this function is also valid for other types with operators % and /:
//return gcd(fac0,fac1)=coef0*fac0+coef1*fac1:
template<typename FactorType>
FactorType EuclideanDivision
(FactorType const & fac0,FactorType const & fac1,FactorType * coef0=NULL,FactorType * coef1=NULL);

//}(general arithmetic)
/****************************************************************************************************/
//functions for the integral number types{

//Take modulo:
template<typename IntType=Type_Int,typename ModIntType=Type_UInt>
ModIntType mod(IntType const & number,ModIntType const & modulo);
//Take division for integral:
template<typename IntType=Type_Int,typename ModIntType=Type_UInt>
IntType divint(IntType const & number,ModIntType const & modulo);

//Is co-prime:
template<typename IntType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Type_Bool Is_coprime(IntType const & n0,IntType const & n1);

//}(functions for the integral number types)
/****************************************************************************************************/
//functions for the real number types{

//Check: a real number is an integer:
template<typename RealType>
inline Type_Bool constexpr Is_integer(RealType const & number){
	LZ_DEF_func_check_traits(std::is_arithmetic<RealType>::value);
	return std::is_integral<RealType>::value ? true : number==floor(number);
}

//Absolute value function for the real number types:
template<typename RealType>
inline RealType constexpr abs(RealType const & number){
	LZ_DEF_func_check_traits(std::is_arithmetic<RealType>::value);
	return number>=(RealType)0 ? number : -number ;
}

//Sign function for the real number types:
template<typename RealType>
inline RealType constexpr sign(RealType const & number){
	LZ_DEF_func_check_traits(std::is_arithmetic<RealType>::value);
	return number>(RealType)0 ? (RealType)1 :
		(number<(RealType)0 ? (RealType)(-1) : (RealType)0);
}

#if LZ_DEF_compile_ignore //unrealised;
//Check: two objects are approximately equal:
template<typename Tp0,typename Tp1,typename RealType=Type_Real>
Type_Bool Is_approx
(Tp0 const & num0,Tp1 const & num1,RealType pre=(RealType)LZ_DEF_const_default_precision)
#endif //#if LZ_DEF_compile_ignore

//}(functions for the real number types)
/****************************************************************************************************/
//functions for matrices{
#ifdef LZ_DEF_extLIB_Eigen

//Vectorisation for a matrix:
template<typename EigenDerived>
Type_MatTemp<typename Eigen::DenseBase<EigenDerived>::Scalar>
Vec(Eigen::DenseBase<EigenDerived> const & X);

//Converting the `Scalar' type of a matrix:
template<typename NewScalar=Liuze::Type_Real,typename EigenDerived>
Type_MatTemp<NewScalar>
Matrix_Scalar_cast(Eigen::DenseBase<EigenDerived> const & X);

#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(functions for matrices)
/****************************************************************************************************/
//functions to calculate all kinds of combination numbers{
/*
template parameters:
IntType: integral type.
UIntType: unsigned integral type.
ResType: the type of results, should be able to represent more larger numbers then UIntType.
*/
namespace Comb_num{
	
	//the combination number:
	template<typename ResType=Type_Real,typename UIntType=Type_UInt>
	ResType comb(UIntType const & num_total,UIntType const & num_select);
	//the permutation or arrange number:
	template<typename ResType=Type_Real,typename UIntType=Type_UInt,
		typename=typename std::enable_if<std::is_arithmetic<ResType>::value &&
		std::is_integral<UIntType>::value>::type>
	ResType perm(UIntType const & num_total,UIntType const & num_select);
	
} //namespace Comb_num;
//}(functions to calculate all kinds of combination numbers)
/****************************************************************************************************/
//class Comb_enum{
/*
Base class to enumerate values of all kinds of combinations.
Itself is its iterator.
UIntType: unsigned integral type.
*/
template<typename UIntType=Type_UInt,
	typename=typename std::enable_if<std::is_integral<UIntType>::value,void>::type>
class Comb_enum{
	public:
		//typedef: this class:
		typedef Comb_enum<UIntType> type_this;
		//typedef: the base class:
		typedef void type_base;
		//typedef: the type of nature number:
		typedef UIntType type_nat;
		//typedef: the type of state:
		typedef Type_Stat type_stat;
		//typedef: the type of value of the combination:
		typedef Type_ArrayTemp<type_nat> type_val;
		//typedef: its iterator:
		typedef type_this iterator;
		//typedefs for iterator:
		typedef std::iterator<std::forward_iterator_tag,type_val const> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//values:
		static type_stat constexpr val_stat_deref=0; //de-reference-able;
		static type_stat constexpr val_stat_end=1; //at after the last entry;
		static type_stat constexpr val_stat_rend=2; //at before the first entry;
		static type_stat constexpr val_stat_null=10; //the set has no entry;
		
		Comb_enum();
		Comb_enum(type_this const & comb_enum);
		Comb_enum(type_this && comb_enum);
		
		virtual ~Comb_enum()=default;
		
		//operator:
		type_this & operator =(type_this const & comb_enum);
		type_this & operator =(type_this && comb_enum);
		Type_Bool operator ==(type_this const & comb_enum) const;
		explicit operator Type_Bool() const;
		virtual value_type & operator *() const final;
		virtual pointer operator ->() const final;
		virtual type_this & operator ++()=0;
		
		/*
		//virtual: deleted in type_this because of the return type{
		//operator ++(int):
		type_this operator ++(int);
		//Get: the iterator to begin:
		iterator begin() const;
		//Get: the iterator past to end:
		iterator end() const;
		//}(virtual: deleted in type_this because of the return type)
		*/
		
		//Is:
		Type_Bool Is_stat_null() const;
		Type_Bool Is_stat_deref() const;
		Type_Bool Is_stat_end() const;
		Type_Bool Is_stat_rend() const;
		
		//Set: the state:
		virtual type_stat Set_stat_null();
		virtual type_stat Set_stat_deref();
		virtual type_stat Set_stat_end();
		virtual type_stat Set_stat_rend();
		
		//Get: the state:
		type_stat Get_stat() const;
		//Get: get the value:
		type_val Get_val() const;
		
		//Get: get the next value:
		virtual type_val Get_val_next() final;
		//Get: get the value after num_step steps:
		virtual type_val Get_val_next(type_nat const & num_step) final;
		
		//Set: re-initialise:
		virtual type_stat Set_begin()=0;
	protected:
		type_stat stat; //state: see this->val_stat_...;
		type_val val;
};
//}(class Comb_enum)
/****************************************************************************************************/
//class Comb_enum_pow{
//enumerate all elements in the set of the form n^m:
template<typename UIntType=Type_UInt>
class Comb_enum_pow: public Comb_enum<UIntType>{
	public:
		//typedef: this class:
		typedef Comb_enum_pow<UIntType> type_this;
		//typedef: the base class:
		typedef Comb_enum<UIntType> type_base;
		//typedef: the type of nature number:
		typedef typename type_base::type_nat type_nat;
		//typedef: the type of state:
		typedef typename type_base::type_stat type_stat;
		//typedef: the type of value of the combination:
		typedef typename type_base::type_val type_val;
		//typedef: its iterator:
		typedef type_this iterator;
		//typedefs for iterator:
		typedef std::forward_iterator_tag iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		Comb_enum_pow();
		Comb_enum_pow(type_nat const & num_base,type_nat const & num_power);
		Comb_enum_pow(type_this const & comb_enum_pow);
		Comb_enum_pow(type_this && comb_enum_pow);
		virtual ~Comb_enum_pow()=default;
		
		//operator and virtual: from type_base:
		type_this& operator =(type_this const & comb_enum_pow);
		type_this& operator =(type_this && comb_enum_pow);
		Type_Bool operator ==(type_this const & comb_enum_pow) const;
		virtual type_this& operator ++();
		virtual type_stat Set_begin();
		
		//virtual: deleted in type_base because of the return type:
		iterator begin() const;
		iterator end() const;
		
		//Set: re-initialise with new parameters:
		type_stat Set_begin(type_nat const & num_base,type_nat const & num_power);
		type_stat Set_begin(type_this const & comb_enum_pow);
		
		//Get: the parameters:
		type_nat Get_base() const;
		type_nat Get_power() const;
		
		//static:
		template<typename IndexType=type_nat,
			typename=typename std::enable_if<std::is_arithmetic<IndexType>::value>::type>
		static type_stat value_to_index
		(type_nat const & num_base,type_val const & val_pow,IndexType & index);
		static type_stat index_to_value
		(type_nat const & num_base,type_nat const & num_power,type_nat const & index,
		type_val & val_pow);
	protected:
		type_nat base,power;
};
//}(class Comb_enum_pow)
/****************************************************************************************************/
//class Comb_enum_cart{
//enumerate all elements in the Cartesian product of several nature numbers:
template<typename UIntType=Type_UInt>
class Comb_enum_cart: public Comb_enum<UIntType>{
	public:
		//typedef: this class:
		typedef Comb_enum_cart<UIntType> type_this;
		//typedef: the base class:
		typedef Comb_enum<UIntType> type_base;
		//typedef: the type of nature number:
		typedef typename type_base::type_nat type_nat;
		//typedef: the type of a sequence of nature numbers:
		typedef Type_ArrayTemp<type_nat> type_nat_seq;
		//typedef: the type of state:
		typedef typename type_base::type_stat type_stat;
		//typedef: the type of value of the combination:
		typedef typename type_base::type_val type_val;
		//typedef: its iterator:
		typedef type_this iterator;
		//typedefs for iterator:
		typedef std::forward_iterator_tag iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		Comb_enum_cart();
		Comb_enum_cart(type_nat_seq const & num_base);
		Comb_enum_cart(type_this const & comb_enum_cart);
		Comb_enum_cart(type_this && comb_enum_cart);
		virtual ~Comb_enum_cart()=default;
		
		//operator and virtual: from type_base:
		type_this& operator =(type_this const & comb_enum_cart);
		type_this& operator =(type_this && comb_enum_cart);
		Type_Bool operator ==(type_this const & comb_enum_cart) const;
		virtual type_this& operator ++();
		virtual type_stat Set_begin();
		
		//virtual: deleted in type_base because of the return type:
		iterator begin() const;
		iterator end() const;
		
		//Set: re-initialise with new parameters:
		type_stat Set_begin(type_nat_seq const & num_base);
		type_stat Set_begin(type_this const & comb_enum_cart);
		
		//Get: the parameters:
		type_nat Get_base() const;
	protected:
		type_nat_seq base;
};
//}(class Comb_enum_cart)
/****************************************************************************************************/
//class Comb_enum_comb{
//enumerate all elements in the set of the form {s : s is a subset of n and Card(s)=m}:
template<typename UIntType=Type_UInt>
class Comb_enum_comb: public Comb_enum<UIntType>{
	public:
		//typedef: this class:
		typedef Comb_enum_comb<UIntType> type_this;
		//typedef: the base class:
		typedef Comb_enum<UIntType> type_base;
		//typedef: the type of nature number:
		typedef typename type_base::type_nat type_nat;
		//typedef: the type of state:
		typedef typename type_base::type_stat type_stat;
		//typedef: the type of value of the combination:
		typedef typename type_base::type_val type_val;
		//typedef: its iterator:
		typedef type_this iterator;
		//typedefs for iterator:
		typedef std::forward_iterator_tag iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		Comb_enum_comb();
		Comb_enum_comb(type_nat const & num_total,type_nat const & num_select);
		Comb_enum_comb(type_this const & comb_enum_comb);
		Comb_enum_comb(type_this && comb_enum_comb);
		virtual ~Comb_enum_comb()=default;
		
		//operator and virtual: from type_base:
		type_this& operator =(type_this const & comb_enum_comb);
		type_this& operator =(type_this && comb_enum_comb);
		Type_Bool operator ==(type_this const & comb_enum_comb) const;
		virtual type_this& operator ++();
		virtual type_stat Set_begin();
		
		//virtual: deleted in type_base because of the return type:
		iterator begin() const;
		iterator end() const;
		
		//Set: re-initialise with new parameters:
		type_stat Set_begin(type_nat const & num_total,type_nat const & num_select);
		type_stat Set_begin(type_this const & comb_enum_comb);
		
		//Get: the parameters:
		type_nat Get_total() const;
		type_nat Get_select() const;
		
		//static:
		template<typename IndexType=type_nat,
			typename=typename std::enable_if<std::is_arithmetic<IndexType>::value>::type>
		static type_stat value_to_index
		(type_nat const & num_total,type_val const & val_comb,IndexType & index);
		static type_stat index_to_value
		(type_nat const & num_total,type_nat const & num_select,type_nat const & index,
		type_val & val_comb);
	protected:
		type_nat total,select;
};
//}(class Comb_enum_comb)
/****************************************************************************************************/
//class Comb_enum_perm{
//enumerate all permutations on {0,...,num_order-1}:
template<typename UIntType=Type_UInt>
class Comb_enum_perm: public Comb_enum<UIntType>{
	public:
		//typedef: this class:
		typedef Comb_enum_perm<UIntType> type_this;
		//typedef: the base class:
		typedef Comb_enum<UIntType> type_base;
		//typedef: the type of nature number:
		typedef typename type_base::type_nat type_nat;
		//typedef: the type of state:
		typedef typename type_base::type_stat type_stat;
		//typedef: the type of value of the combination:
		typedef typename type_base::type_val type_val;
		//typedef: its iterator:
		typedef type_this iterator;
		//typedefs for iterator:
		typedef std::forward_iterator_tag iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		Comb_enum_perm();
		Comb_enum_perm(type_nat const & num_order);
		Comb_enum_perm(type_this const & comb_enum_perm);
		Comb_enum_perm(type_this && comb_enum_perm);
		virtual ~Comb_enum_perm()=default;
		
		//operator and virtual: from type_base:
		type_this& operator =(type_this const & comb_enum_perm);
		type_this& operator =(type_this && comb_enum_perm);
		Type_Bool operator ==(type_this const & comb_enum_perm) const;
		virtual type_this& operator ++();
		virtual type_stat Set_begin();
		
		//virtual: deleted in type_base because of the return type:
		iterator begin() const;
		iterator end() const;
		
		//Set: re-initialise with new parameters:
		type_stat Set_begin(type_nat const & num_order);
		type_stat Set_begin(type_this const & comb_enum_perm);
		
		//Get: the parameters:
		type_nat Get_order() const;
		
		//static:
		template<typename IndexType=type_nat,
			typename=typename std::enable_if<std::is_arithmetic<IndexType>::value>::type>
		static type_stat value_to_index
		(type_val const & val_perm,IndexType & index);
		static type_stat index_to_value
		(type_nat const & num_order,type_nat const & index,type_val & val_perm);
	protected:
		type_nat order;
};
//}(class Comb_enum_perm)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/function.cpp"

#endif //#ifndef LZ_DEF_LZ_H_math_function