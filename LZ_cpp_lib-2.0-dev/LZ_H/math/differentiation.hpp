#ifndef LZ_DEF_LZ_H_math_differentiation
#define LZ_DEF_LZ_H_math_differentiation 202105L

#include"LZ_H/math/function.hpp"

#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace math{

/****************************************************************************************************/
template<typename RealType=Type_Real>
class Derivative{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		typedef Derivative<RealType> type_this;
		typedef RealType type_real;
		typedef Type_Size type_size;
		typedef Type_MatTemp<type_real> type_matrix;
		typedef std::function<type_matrix(type_matrix const &)> type_func;
		typedef Type_Stat type_stat;
		//type: for storing the precision information:
		class type_pre{
			public:
				//values:
				type_real err,length_init; //the maximum error and the initial step length;
				type_size n_iter_max; //the maximum iteration number;
				type_stat err_type; //the error types;
				
				//static values:
				static type_stat constexpr val_err_invalid=0; //as if the error is always very large;
				static type_stat constexpr val_err_in_abs=1; //absolute error of input;
				static type_stat constexpr val_err_in_rel=2; //relative error of input;
				static type_stat constexpr val_err_out_abs=3; //absolute error of output;
				static type_stat constexpr val_err_out_rel=4; //relative error of output;
				static type_size constexpr val_size_infinity=0;
				
				//Constructor:
				type_pre
					(type_real length_initial=1,type_real error_max=LZ_DEF_const_default_precision,
					type_size iteration_max=100,type_stat error_type=type_pre::val_err_out_abs);
				type_pre(type_pre const & precision)=default;
				
				//operator:
				type_pre & operator =(type_pre const & precision)=default;
		}; //class type_pre;
	public:
		//Construction:
		Derivative();
		template<typename FuncType,
			typename=typename std::enable_if<std::is_convertible<FuncType&&,type_func>::value>::type>
		Derivative
			(FuncType && func,
			type_size nrow_out=0,type_size ncol_out=0,type_size nrow_in=0,type_size ncol_in=0);
		Derivative(type_this const & obj)=delete;
		Derivative(type_this && obj)=delete;
		
		//Destructor:
		~Derivative();
		
		//operator:
		type_this & operator =(type_this const & obj)=delete;
		type_this & operator =(type_this && obj)=delete;
		//The following operator () calculates the derivative numerically.
		//If the primitive function is from R^(m cross n) to R^(p cross q),
		//then the derivative is from R^(m cross n) to R^(mp cross nq);
		//i.e., if the primitive is x |-> y, then the derivative is x |-> z,
		//with z(im+k,jn+l) = (partial y(i,j))/(partial x(k,l)).
		type_matrix operator ()(type_matrix const & x);
		
		//Get:
		type_matrix Get_function(type_matrix const & x);
		
		//Set:
		void Set_function();
		template<typename FuncType,
			typename=typename std::enable_if<std::is_convertible<FuncType&&,type_func>::value>::type>
		void Set_function
			(FuncType && func,
			type_size nrow_out=0,type_size ncol_out=0,type_size nrow_in=0,type_size ncol_in=0);
		void Set_precision
			(type_real length_initial=1,type_real error_max=LZ_DEF_const_default_precision,
			type_size iteration_max=100,type_stat error_type=type_pre::val_err_out_abs);
	protected:
		//values:
		type_size in_nrow,in_ncol,out_nrow,out_ncol,der_nrow,der_ncol; //dim.'s of input, output, derivative;
		type_func * pt_func; //the primitive function;
		type_stat stat; //the state;
		type_pre pre; //the precision information;
		
		//static values:
		static type_stat constexpr val_stat_null=0; //without primitive function;
		static type_stat constexpr val_stat_complete=1; //with primitive function, with dimension;
		static type_stat constexpr val_stat_incomplete=2; //with primitive function, without dimension;
		
		//Cal:
		Type_Bool Cal_dim(type_matrix const & x);
		Type_Bool Cal_err
			(type_real const & x0,type_real const & x1,type_real const & y0,type_real const & y1) const;
}; //class Derivative;
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/differentiation.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_math_differentiation