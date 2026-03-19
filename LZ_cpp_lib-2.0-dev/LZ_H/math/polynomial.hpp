#ifndef LZ_DEF_LZ_H_math_polynomial
#define LZ_DEF_LZ_H_math_polynomial 202105L

#include<iostream>
#include"LZ_H/math/function.hpp"

namespace Liuze{
namespace math{
/****************************************************************************************************/
//class Polynomial:
template<typename CoefType>
class Polynomial{
	public:
		typedef Polynomial<CoefType> type_this;
		typedef CoefType type_coef;
		typedef Type_Size type_size;
		typedef Type_ArrayTemp<type_coef> type_store;
		
		//Constructor:
		Polynomial();
		template<typename InputIterType,
			typename=typename std::enable_if<
			std::is_convertible<typename InputIterType::value_type,type_coef>::value>::type>
		Polynomial(InputIterType iter_coef_0,InputIterType const & iter_coef_end);
		template<typename InputIterType,
			typename=typename std::enable_if<
			std::is_convertible<typename InputIterType::value_type,type_coef>::value>::type>
		Polynomial(InputIterType iter_coef_0,type_size const & length);
		template<typename CoefVecType,
			typename=typename std::enable_if<
			std::is_convertible<typename CoefVecType::value_type,type_coef>::value>::type>
		explicit Polynomial(CoefVecType const & coef_vector);
		Polynomial(type_this const & poly);
		Polynomial(type_this && poly);
		
		//Destructor:
		~Polynomial()=default;
		
		//operator:
		type_this & operator =(type_this const & poly);
		type_this & operator =(type_this && poly);
		template<typename CoefVecType,
			typename=typename std::enable_if<
			std::is_convertible<typename CoefVecType::value_type,type_coef>::value>::type>
		type_this & operator =(CoefVecType const & coef_vector);
		
		type_this & operator +=(type_this const & poly);
		type_this & operator -=(type_this const & poly);
		type_this & operator *=(type_this const & poly);
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this & operator *=(AuxIntType const & num_times);
		type_this & operator *=(type_coef const & val_coef);
		type_this & operator /=(type_this const & poly);
		type_this & operator %=(type_this const & poly);
		type_this operator +() const;
		type_this operator -() const;
		
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator +
		(Polynomial<AuxCoefType> poly0,Polynomial<AuxCoefType> const & poly1);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator -
		(Polynomial<AuxCoefType> poly0,Polynomial<AuxCoefType> const & poly1);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator *
		(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator *
		(typename Polynomial<AuxCoefType>::type_coef const & val_coef,Polynomial<AuxCoefType> poly);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator *
		(Polynomial<AuxCoefType> poly,typename Polynomial<AuxCoefType>::type_coef const & val_coef);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator /
		(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1);
		template<typename AuxCoefType>
		friend Polynomial<AuxCoefType> operator %
		(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1);
		
		Type_Bool operator ==(type_this const & poly) const;
		Type_Bool operator !=(type_this const & poly) const;
		Type_Bool operator !() const;
		explicit operator Type_Bool() const;
		
		//Evaluating at `x':
		template<typename ValType>
		ValType operator ()(ValType const & x) const;
		
		//Get coefficient (not safe):
		type_coef const & operator [](type_size const & i_deg) const;
		
		//Is:
		Type_Bool Is_zero() const;
		
		//value:
		//The order of a polynomial is degree plus 1:
		type_size order() const;
		type_size degree() const;
		//Get the zero of the coefficient type:
		type_coef coef_zero() const;
		//Get coefficient (safe):
		type_coef coef(type_size const & i_deg) const;
		
		//other operator:
		//Get the zero polynomial:
		type_this zero() const;
		//Get the polynomial (`coeff_deg0'):
		type_this constant(type_coef const & coeff_deg0) const;
		template<typename InputIterType,
			typename=typename std::enable_if<
			std::is_convertible<typename InputIterType::value_type,type_coef>::value>::type>
		type_this & assign(InputIterType iter_coef_0,InputIterType const & iter_coef_end);
		template<typename InputIterType,
			typename=typename std::enable_if<
			std::is_convertible<typename InputIterType::value_type,type_coef>::value>::type>
		type_this & assign(InputIterType iter_coef_0,type_size const & length);
		template<typename CoefVecType,
			typename=typename std::enable_if<
			std::is_convertible<typename CoefVecType::value_type,type_coef>::value>::type>
		type_this & assign(CoefVecType const & coef_vector);
		//Get the `num_exp'-th power of this polynoimal:
		template<typename AuxIntType,
			typename=typename std::enable_if<std::is_integral<AuxIntType>::value>::type>
		type_this pow(AuxIntType const & num_exp) const;
		
		//Set:
		type_this & Set_coef(type_size const & i_deg,type_coef const & val_coef);
		
		//input and output:
		template<typename AuxCoefType>
		friend std::ostream & operator <<
		(std::ostream & out,Polynomial<AuxCoefType> const & poly);
		
		//static:
		static Type_Stat div
		(type_this const & numer,type_this const & denom,
		type_this * ptr_quot=NULL,type_this * ptr_rem=NULL);
	protected:
		//member value:
		type_store m_arr_coef;
		type_size m_order;
		
		//operator:
		type_coef & operator [](type_size const & i_deg);
}; //class Polynomial;

//operators of class Polynomial{

template<typename IntType,typename CoefType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Polynomial<CoefType> operator *
(IntType const & num_times,Polynomial<CoefType> poly);
template<typename IntType,typename CoefType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Polynomial<CoefType> operator *
(Polynomial<CoefType> poly,IntType const & num_times);

//}(operators of class Polynomial)

//arithmetic functions{

//pseudo-specialisations for functions\
//`zero', `one', `Is_approx_zero', `powint'\
//in namespace `Liuze::math';

template<typename CoefType>
Type_Bool div
(Polynomial<CoefType> const & numer,Polynomial<CoefType> const & denom,
Polynomial<CoefType> * ptr_quot=NULL,Polynomial<CoefType> * ptr_rem=NULL);

//}(arithmetic functions)

//(class Polynomial)
/****************************************************************************************************/
} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/polynomial.cpp"

#endif //#ifndef LZ_DEF_LZ_H_math_polynomial