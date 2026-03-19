#ifndef LZ_DEF_LZ_H_math_GaloisField
#define LZ_DEF_LZ_H_math_GaloisField 202105L

#include<LZ_H/math/polynomial.hpp>
#include<LZ_H/math/intresidue.hpp>

/****************************************************************************************************/
namespace Liuze{
namespace math{
namespace GaloisField{

/****************************************************************************************************/
//classes in this hpp{
template<typename IntType=Type_Int,
	typename CoefRingType=Liuze::math::IntegerResidueClass::Ring<IntType>,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
class Field_Base;

template<typename IntType=Type_Int,
	typename CoefRingType=Liuze::math::IntegerResidueClass::Ring<IntType> >
using Field=Field_Base<IntType,CoefRingType>;

template<typename IntType=Type_Int,
	typename CoefRingType=Liuze::math::IntegerResidueClass::Ring<IntType> >
class Field_tabulated;

template<typename Tp>
class is_field_type;

template<typename IntType=Type_Int,typename FieldType=Field<IntType>,
	typename=typename std::enable_if<std::is_integral<IntType>::value &&
	//`std::is_base_of<Field_Base<IntType,typename FieldType::type_ring_coef>,FieldType>::value>::type>
	is_field_type<FieldType>::value>::type>
class Element;

//`template<typename IntType,typename CoefRingType>
//`class Element<IntType,Field_tabulated<IntType,CoefRingType> >;
//}(classes in this hpp)
/****************************************************************************************************/

/****************************************************************************************************/
//the field:

template<typename Tp>
class is_field_type: public ::std::false_type{
};

//`template<typename IntType=Type_Int,
//`	typename CoefRingType=Liuze::math::IntegerResidueClass::Ring<IntType>,
//`	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
template<typename IntType,typename CoefRingType,typename>
class Field_Base{
	public:
		typedef Field_Base<IntType,CoefRingType> type_this;
		typedef IntType type_int;
		typedef CoefRingType type_ring_coef;
		typedef Liuze::math::IntegerResidueClass::Element<type_int,type_ring_coef> type_int_coef;
		typedef Polynomial<type_int_coef> type_poly;
		typedef Element<type_int,type_this> type_element;
		typedef typename type_int_coef::type_val type_val;
		typedef type_val value_type;
		
		//Constructor:
		Field_Base();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Field_Base
			(AuxIntType const & num_prime,AuxIntType const & num_power,type_poly const & irr_poly);
		Field_Base(type_this const & field);
		Field_Base(type_this && field);
		
		//Destructor:
		~Field_Base()=default;
		
		//operator:
		type_this & operator =(type_this const & field);
		type_this & operator =(type_this && field);
		Type_Bool operator ==(type_this const & field) const;
		Type_Bool operator !=(type_this const & field) const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element operator ()(AuxIntType const & num_value) const;
		
		//Is:
		Type_Bool Is_valid() const;
		Type_Bool Is_invalid() const;
		Type_Bool Is_universal() const;
		Type_Bool Is_specific() const;
		
		//value:
		type_val size() const;
		type_val characteristic() const;
		type_val power() const;
		type_poly irrPoly() const;
		type_poly const & irrPoly_cref() const;
		
		//other function:
		type_element element_zero() const;
		type_element element_one() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element element(AuxIntType const & num_value) const;
		type_val poly_to_num(type_poly const & poly_value) const;
		type_poly num_to_poly(type_val const & num_value) const;
		
		//static:
		static type_this invalid();
		static type_this universal();
		static type_element element_invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_element element_universal(AuxIntType const & num_value);
		//fun: Get_irrPoly_all:
		//it returns all irreducible monic polynomials
		//whose coefficients are on the ring Z/(`num_prime'*Z)
		//and degrees are from 1 to `num_degree_upper':
		static Type_ArrayTemp<type_poly> Get_irrPoly_all
			(type_val const & num_prime,type_val const & num_degree_upper);
	protected:
		type_val m_ch,m_pow,m_size; //GF(m_ch^m_pow), where m_ch is a prime;
		type_poly m_irr; //the irreducible polynomial of degree m_pow;
}; //class Field_Base;

template<typename IntType,typename CoefRingType,typename TraitsCheck>
class is_field_type<Field_Base<IntType,CoefRingType,TraitsCheck> >:
public ::std::true_type{
};

//`template<typename IntType=Type_Int,
//`	typename CoefRingType=Liuze::math::IntegerResidueClass::Ring<IntType> >
template<typename IntType,typename CoefRingType>
class Field_tabulated: public Field_Base<IntType,CoefRingType>{
	public:
		typedef Field_tabulated<IntType,CoefRingType> type_this;
		typedef Field_Base<IntType,CoefRingType> type_base;
		typedef typename type_base::type_int type_int;
		typedef typename type_base::type_ring_coef type_ring_coef;
		typedef typename type_base::type_int_coef type_int_coef;
		typedef typename type_base::type_poly type_poly;
		typedef Element<type_int,type_this> type_element;
		typedef typename type_base::type_val type_val;
		typedef typename type_base::value_type value_type;
		typedef Type_ArrayTemp<type_val> type_table;
		
		//Constructor:
		Field_tabulated();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Field_tabulated
			(AuxIntType const & num_prime,AuxIntType const & num_power,type_poly const & irr_poly);
		template<typename FieldType,
			typename=typename std::enable_if<std::is_base_of<type_base,FieldType>::value>::type>
		Field_tabulated(FieldType && field);
		Field_tabulated(type_this const & field);
		Field_tabulated(type_this && field);
		
		//Destructor:
		~Field_tabulated()=default;
		
		//operator:
		template<typename FieldType,
			typename=typename std::enable_if<std::is_base_of<type_base,FieldType>::value>::type>
		type_this & operator =(FieldType && field);
		type_this & operator =(type_this const & field);
		type_this & operator =(type_this && field);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element operator ()(AuxIntType const & num_value) const;
		
		//value:
		type_val const & table_add(type_val const & id) const;
		type_val const & table_neg(type_val const & id) const;
		type_val const & table_pow(type_val const & id) const;
		type_val const & table_log(type_val const & id) const;
		
		//other function:
		type_element primitive() const;
		type_element element_zero() const;
		type_element element_one() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_element element(AuxIntType const & num_value) const;
		
		//static:
		static type_this invalid();
		static type_this universal();
		static type_element element_invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_element element_universal(AuxIntType const & num_value);
	protected:
		//value:
		type_table m_tab_add; //addition table of length `m_size*(m_size+1)/2';
		type_table m_tab_neg; //stores the negatives;
		type_table m_tab_pow; //stores the powers of a primitive;
		type_table m_tab_log; //stores the logarithms base the primitive;
		
		//Cal:
		void Cal_table();
}; //class Field_tabulated;

template<typename IntType,typename CoefRingType>
class is_field_type<Field_tabulated<IntType,CoefRingType> >:
public ::std::true_type{
};

//(the field)
/****************************************************************************************************/
/****************************************************************************************************/
//the element:

//`template<typename IntType=Type_Int,typename FieldType=Field<IntType>,
//`	typename=typename std::enable_if<std::is_integral<IntType>::value &&
//`	std::is_base_of<Field_Base<IntType,typename FieldType::type_ring_coef>,FieldType>::value>::type>
template<typename IntType,typename FieldType,typename>
class Element{
	LZ_DEF_func_check_traits
		((std::is_same<Element<IntType,FieldType>,typename FieldType::type_element>::value));
	public:
		typedef Element<IntType,FieldType> type_this;
		typedef FieldType type_field;
		typedef type_field type_space;
		typedef typename type_field::type_int type_int;
		typedef typename type_field::type_int_coef type_int_coef;
		typedef typename type_field::type_poly type_poly;
		typedef typename type_field::type_val type_val;
		typedef typename type_field::value_type value_type;
		typedef std::shared_ptr<type_field const> type_ptr_field;
		typedef type_ptr_field type_ptr_space;
		
		//Constructor:
		Element();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		explicit Element(AuxIntType const & num_value);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_field const & field);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_ptr_field const & ptr_field);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_this const & element);
		Element(type_this const & element);
		Element(type_this && element);
		
		//Destructor:
		~Element()=default;
		
		//operator:
		type_this & operator =(type_this const & element);
		type_this & operator =(type_this && element);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_this & operator =(AuxIntType const & num_value);
		
		type_this & operator +=(type_this const & x);
		type_this & operator -=(type_this const & x);
		type_this & operator *=(type_this const & x);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this & operator *=(AuxIntType const num_times);
		type_this & operator /=(type_this const & x);
		
		type_this operator +() const;
		type_this operator -() const;
		
		Type_Bool operator ==(type_this const & element) const;
		Type_Bool operator !=(type_this const & element) const;
		Type_Bool operator !() const;
		explicit operator Type_Bool() const;
		
		template<typename AuxType,
			typename=typename std::enable_if<std::is_arithmetic<AuxType>::value &&
			std::is_convertible<type_val,AuxType>::value>::type>
		explicit operator AuxType() const;
		
		//Is:
		Type_Bool Is_invalid() const;
		Type_Bool Is_valid() const;
		Type_Bool Is_universal() const;
		Type_Bool Is_specific() const;
		Type_Bool Is_zero() const;
		Type_Bool Is_one() const;
		
		//value:
		type_val value_number() const;
		type_poly value_polynomial() const;
		type_field const & space() const;
		type_ptr_field const & ptr_space() const;
		
		//other operator:
		type_this zero() const;
		type_this one() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_this element(AuxIntType const & num_value) const;
		type_this inv() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this pow(AuxIntType const & num_exp) const;
		
		//Set:
		type_this & Set_invalid();
		
		//other function:
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,
				typename std::decay<typename InputIterType::value_type>::type
				>::value &&
				std::is_arithmetic<typename OutputIterType::value_type>::value &&
				std::is_convertible<type_val,typename OutputIterType::value_type>::value
			>::type>
		void convert_to_value_number
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,
				typename std::decay<typename InputIterType::value_type>::type
				>::value &&
				std::is_convertible<type_poly,typename OutputIterType::value_type>::value
			>::type>
		void convert_to_value_polynomial
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,typename OutputIterType::value_type>::value &&
				std::is_convertible<typename InputIterType::value_type,type_val>::value
			>::type>
		void convert_from_value_number
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		
		//static:
		static type_this invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_this universal(AuxIntType const & num_value);
	protected:
		//member value:
		type_ptr_field m_ptr_field;
		type_val m_val_num;
		type_poly m_val_poly;
		
		//static value:
		static type_ptr_field val_ptr_field_invalid;
		static type_ptr_field val_ptr_field_universal;
	private:
		//other function:
		void Set_init_val();
}; //class Element;

//operators of class Element{

template<typename IntType,typename FieldType>
Element<IntType,FieldType> operator +
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1);
template<typename IntType,typename FieldType>
Element<IntType,FieldType> operator -
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1);
template<typename IntType,typename FieldType>
Element<IntType,FieldType> operator *
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1);
template<typename IntType,typename FieldType,typename AuxIntType,
	typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
Element<IntType,FieldType> operator *
(Element<IntType,FieldType> x,AuxIntType const & num_times);
template<typename IntType,typename FieldType,typename AuxIntType,
	typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
Element<IntType,FieldType> operator *
(AuxIntType const & num_times,Element<IntType,FieldType> x);
template<typename IntType,typename FieldType>
Element<IntType,FieldType> operator /
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1);

template<typename IntType,typename FieldType>
std::ostream & operator <<
(std::ostream & out,Element<IntType,FieldType> const & element);

//}(operators of class Element)

template<typename IntType,typename CoefRingType>
class Element<IntType,Field_tabulated<IntType,CoefRingType> >{
	public:
		typedef Element<IntType,Field_tabulated<IntType,CoefRingType> > type_this;
		typedef Field_tabulated<IntType,CoefRingType> type_field;
		typedef type_field type_space;
		typedef typename type_field::type_int type_int;
		typedef typename type_field::type_int_coef type_int_coef;
		typedef typename type_field::type_poly type_poly;
		typedef typename type_field::type_val type_val;
		typedef typename type_field::value_type value_type;
		typedef std::shared_ptr<type_field const> type_ptr_field;
		typedef type_ptr_field type_ptr_space;
		
		//Constructor:
		Element();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		explicit Element(AuxIntType const & num_value);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_field const & field);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_ptr_field const & ptr_field);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		Element(AuxIntType const & num_value,type_this const & element);
		Element(type_this const & element);
		Element(type_this && element);
		
		//Destructor:
		~Element()=default;
		
		//operator:
		type_this & operator =(type_this const & element);
		type_this & operator =(type_this && element);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_this & operator =(AuxIntType const & num_value);
		
		type_this & operator +=(type_this const & x);
		type_this & operator -=(type_this const & x);
		type_this & operator *=(type_this const & x);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this & operator *=(AuxIntType const num_times);
		type_this & operator /=(type_this const & x);
		
		type_this operator +() const;
		type_this operator -() const;
		
		Type_Bool operator ==(type_this const & element) const;
		Type_Bool operator !=(type_this const & element) const;
		Type_Bool operator !() const;
		explicit operator Type_Bool() const;
		
		template<typename AuxType,
			typename=typename std::enable_if<std::is_arithmetic<AuxType>::value &&
			std::is_convertible<type_val,AuxType>::value>::type>
		explicit operator AuxType() const;
		
		//Is:
		Type_Bool Is_invalid() const;
		Type_Bool Is_valid() const;
		Type_Bool Is_universal() const;
		Type_Bool Is_specific() const;
		Type_Bool Is_zero() const;
		Type_Bool Is_one() const;
		
		//value:
		type_val value_number() const;
		type_poly value_polynomial() const;
		type_field const & space() const;
		type_ptr_field const & ptr_space() const;
		
		//other operator:
		type_this zero() const;
		type_this one() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		type_this element(AuxIntType const & num_value) const;
		type_this inv() const;
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this pow(AuxIntType const & num_exp) const;
		
		//Set:
		type_this & Set_invalid();
		
		//other function:
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,
				typename std::decay<typename InputIterType::value_type>::type
				>::value &&
				std::is_arithmetic<typename OutputIterType::value_type>::value &&
				std::is_convertible<type_val,typename OutputIterType::value_type>::value
			>::type>
		void convert_to_value_number
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,
				typename std::decay<typename InputIterType::value_type>::type
				>::value &&
				std::is_convertible<type_poly,typename OutputIterType::value_type>::value
			>::type>
		void convert_to_value_polynomial
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		template<typename InputIterType,typename OutputIterType,
			typename=typename std::enable_if<
				std::is_same<type_this,typename OutputIterType::value_type>::value &&
				std::is_convertible<typename InputIterType::value_type,type_val>::value
			>::type>
		void convert_from_value_number
		(InputIterType const & first,InputIterType const & last,OutputIterType result) const;
		
		//static:
		static type_this invalid();
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_convertible<AuxIntType,type_val>::value>::type>
		static type_this universal(AuxIntType const & num_value);
	protected:
		//member value:
		type_ptr_field m_ptr_field;
		type_val m_val_num;
		
		//static value:
		static type_ptr_field val_ptr_field_invalid;
		static type_ptr_field val_ptr_field_universal;
	private:
		//other function:
		void Set_init_val();
}; //class Element<IntType,Field_tabulated<IntType,CoefRingType> >;

//(the element)
/****************************************************************************************************/

} //namespace GaloisField;
} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace Liuze{
namespace math{

//arithmetic functions for class `GaloisField::Element'{

//pseudo-specialisations for functions\
//`zero', `one', `Is_approx_zero', `inverse', `powint'\
//in namespace `Liuze::math';

//}(arithmetic functions for class `GaloisField::Element')

} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
//numeric_limits{
namespace std{
	template<typename IntType,typename FieldType>
	class numeric_limits<Liuze::math::GaloisField::Element<IntType,FieldType> >;
} //namespace std;
#ifdef LZ_DEF_extLIB_Eigen
namespace Eigen{
	template<typename IntType,typename FieldType>
	class NumTraits<Liuze::math::GaloisField::Element<IntType,FieldType> >;
	template<typename IntType,typename CoefRingType>
	class NumTraits<
		Liuze::math::GaloisField::Element<
		IntType,Liuze::math::GaloisField::Field_tabulated<IntType,CoefRingType> > >;
} //namespace Eigen;
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(numeric_limits)
/****************************************************************************************************/

//implementation:
#include"LZ_H/math/src/GaloisField.cpp"

#endif //#ifndef LZ_DEF_LZ_H_math_GaloisField