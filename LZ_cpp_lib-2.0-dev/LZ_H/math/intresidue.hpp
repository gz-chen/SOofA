#ifndef LZ_DEF_LZ_H_math_intresidue
#define LZ_DEF_LZ_H_math_intresidue 202105L

#include<iostream>
#include"LZ_H/math/function.hpp"

/****************************************************************************************************/
namespace Liuze{
namespace math{
namespace IntegerResidueClass{

/****************************************************************************************************/
//classes in this hpp{
template<typename IntType=Type_Int,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
class Ring_Base;

template<typename IntType>
using Ring=Ring_Base<IntType>;

class tag_type_space_default;

template<typename Tp>
class tag_type_space;

template<typename IntType=Type_Int,typename RingType=Ring<IntType>,
	typename=typename std::enable_if<std::is_integral<IntType>::value,
	typename tag_type_space<RingType>::type>::type>
class Element;
//}(classes in this hpp)
/****************************************************************************************************/

/****************************************************************************************************/
//the ring:

template<typename Tp>
class tag_type_space{
	public:
		static bool constexpr value=false;
};

class tag_type_space_default{
};

//`template<typename IntType=Type_Int,
//`	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
template<typename IntType,typename>
class Ring_Base{
	public:
		//Type: this class:
		typedef Ring_Base<IntType> type_this;
		//Type: integer:
		typedef IntType type_int;
		//Type: signed integer:
		typedef typename std::conditional<std::is_signed<type_int>::value,
			type_int,typename std::make_signed<type_int>::type>::type type_sint;
		//Type: unsigned integer:
		typedef typename std::conditional<std::is_unsigned<type_int>::value,
			type_int,typename std::make_unsigned<type_int>::type>::type type_uint;
		//Type: for storing the value:
		typedef type_sint type_val;
		typedef type_val value_type;
		//Type: for storing the modulo:
		typedef type_sint type_mod;
		//Type: for element:
		typedef Element<type_int,type_this> type_element;
		//Type: the actual ring:
		typedef type_this type_space_act;
		
		//static value:
		static type_mod constexpr val_modulo_universal=-1;
		static type_mod constexpr val_modulo_invalid=-10;
		
		//Constructor:
		Ring_Base();
		template<typename AuxModIntType,
			typename=typename std::enable_if<std::is_convertible<AuxModIntType,type_mod>::value>::type>
		explicit Ring_Base(AuxModIntType const & num_modulo);
		Ring_Base(type_this const & ring);
		Ring_Base(type_this && ring);
		
		//Destructor:
		~Ring_Base()=default;
		
		//operator:
		type_this & operator =(type_this const & ring);
		type_this & operator =(type_this && ring);
		Type_Bool operator ==(type_this const & ring) const;
		Type_Bool operator !=(type_this const & ring) const;
		//Converting `num_value' to the type of element of this ring:
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element operator ()(AuxIntType const & num_value) const;
		
		//Is:
		Type_Bool Is_valid() const;
		Type_Bool Is_invalid() const;
		Type_Bool Is_universal() const;
		Type_Bool Is_specific() const;
		Type_Bool Is_finite() const;
		
		//value:
		type_val size() const;
		type_mod modulo() const;
		type_space_act const & space_act() const;
		type_space_act & space_act();
		
		//other function:
		type_element element_zero() const;
		type_element element_one() const;
		//Converting `num_value' to the type of element of this ring:
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element element(AuxIntType const & num_value) const;
		
		//static:
		static type_this invalid();
		static type_this universal();
		static type_element element_invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_element element_invalid(AuxIntType const & num_value);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_element element_universal(AuxIntType const & num_value);
	protected:
		type_mod m_modulo;
}; //class Ring_Base;

template<typename IntType,typename TraitsCheck>
class tag_type_space<Ring_Base<IntType,TraitsCheck> >{
	public:
		typedef tag_type_space_default type;
		static bool constexpr value=true;
};

//(the ring)
/****************************************************************************************************/
/****************************************************************************************************/
//the element:

//`template<typename IntType=Type_Int,typename RingType=Ring<IntType>,
//`	typename=typename std::enable_if<std::is_integral<IntType>::value,
//`	typename tag_type_space<RingType>::type>::type>
template<typename IntType,typename RingType,typename>
class Element{
	LZ_DEF_func_check_traits
		((std::is_same<Element<IntType,RingType>,typename RingType::type_element>::value));
	public:
		//Type: this class:
		typedef Element<IntType,RingType> type_this;
		//Type: the ring:
		typedef RingType type_space;
		//Type: integer:
		typedef typename type_space::type_int type_int;
		//Type: signed integer:
		typedef typename type_space::type_sint type_sint;
		//Type: unsigned integer:
		typedef typename type_space::type_uint type_uint;
		//Type: for storing the value:
		typedef typename type_space::type_val type_val;
		typedef typename type_space::value_type value_type;
		//Type: for storing the modulo:
		typedef typename type_space::type_mod type_mod;
		//Type: the actual ring:
		typedef typename type_space::type_space_act type_space_act;
		
		//Constructor:
		Element();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		explicit Element(AuxIntType const & num_value);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_space const & ring);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_this const & element);
		Element(type_this const & element);
		Element(type_this && element);
		
		//Destructor:
		~Element()=default;
		
		//operator:
		type_this & operator =(type_this const & x);
		type_this & operator =(type_this && x);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_this & operator =(AuxIntType const & x);
		
		type_this & operator +=(type_this const & x);
		type_this & operator -=(type_this const & x);
		type_this & operator *=(type_this const & x);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this & operator *=(AuxIntType const & num_times);
		type_this & operator /=(type_this const & x);
		
		type_this & operator ++();
		type_this & operator --();
		type_this operator ++(int);
		type_this operator --(int);
		
		type_this operator +() const;
		type_this operator -() const;
		
		template<typename AuxIntType,typename AuxRingType>
		friend Element<AuxIntType,AuxRingType> operator +
		(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1);
		template<typename AuxIntType,typename AuxRingType>
		friend Element<AuxIntType,AuxRingType> operator -
		(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1);
		template<typename AuxIntType,typename AuxRingType>
		friend Element<AuxIntType,AuxRingType> operator *
		(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1);
		template<typename AuxIntType,typename AuxRingType>
		friend Element<AuxIntType,AuxRingType> operator /
		(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1);
		
		Type_Bool operator ==(type_this const & x) const;
		Type_Bool operator !=(type_this const & x) const;
		Type_Bool operator !() const;
		explicit operator Type_Bool() const;
		
		template<typename AuxType,
			typename=typename std::enable_if<std::is_arithmetic<AuxType>::value &&
			std::is_convertible<type_val,AuxType>::value>::type>
		explicit operator AuxType() const;
		
		//Is:
		//It is valid if it is universal or specific:
		Type_Bool Is_valid() const;
		Type_Bool Is_invalid() const;
		Type_Bool Is_universal() const;
		Type_Bool Is_specific() const;
		Type_Bool Is_zero() const;
		//True if it is invertible:
		Type_Bool Is_unit() const;
		
		//member value:
		type_val value() const;
		type_space const & space() const;
		type_space_act const & space_act() const;
		
		//other operator:
		type_this inv() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this pow(AuxIntType const & exp) const;
		
		//other function:
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<typename std::decay<typename InputIterType::value_type>::type,
				type_this>::value &&
				std::is_arithmetic<typename OutputIterType::value_type>::value &&
				std::is_convertible<type_val,typename OutputIterType::value_type>::value
			>::type>
		void convert_to_value
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,typename OutputIterType::value_type>::value &&
				std::is_convertible<typename InputIterType::value_type,type_val>::value
			>::type>
		void convert_from_value
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		
		//input and output:
		template<typename AuxIntType,typename AuxRingType>
		friend std::ostream & operator <<
		(std::ostream & out,Element<AuxIntType,AuxRingType> const & x);
		
		//static:
		static type_this invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_this invalid(AuxIntType const & num_value);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_this universal(AuxIntType const & num_value);
	protected:
		//member value:
		type_space m_ring;
		type_val m_value;
}; //class Element;

//operators of class Element{

template<typename IntType,typename RingType,typename AuxIntType,
	typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
Element<IntType,RingType> operator *
(Element<IntType,RingType> x,AuxIntType const & num_times);
template<typename IntType,typename RingType,typename AuxIntType,
	typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
Element<IntType,RingType> operator *
(AuxIntType const & num_times,Element<IntType,RingType> x);

//}(operators of class Element)

//(the element)
/****************************************************************************************************/

} //namespace IntegerResidueClass;
} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace Liuze{
namespace math{

//arithmetic functions for class `IntegerResidueClass::Element'{

//pseudo-specialisations for functions\
//`zero', `one', `Is_approx_zero', `inverse', `powint'\
//in namespace `Liuze::math';

//}(arithmetic functions for class `IntegerResidueClass::Element')

} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
//numeric_limits{
namespace std{
	template<typename IntType,typename RingType>
	class numeric_limits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >;
} //namespace std;
#ifdef LZ_DEF_extLIB_Eigen
namespace Eigen{
	template<typename IntType,typename RingType>
	class NumTraits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >;
} //namespace Eigen;
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(numeric_limits)
/****************************************************************************************************/

//implementation:
#include"LZ_H/math/src/intresidue.cpp"

#endif //#ifndef LZ_DEF_LZ_H_math_intresidue