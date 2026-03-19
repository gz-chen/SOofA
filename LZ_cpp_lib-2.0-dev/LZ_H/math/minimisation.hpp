#ifndef LZ_DEF_LZ_H_math_minimisation
#define LZ_DEF_LZ_H_math_minimisation 202105L

#include"LZ_H/math/differentiation.hpp"

#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace math{

/****************************************************************************************************/
/****************************************************************************************************/
//iterative descending algorithm for minimisation without or with constraints{
/****************************************************************************************************/
//minimisation algorithm{
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Descent{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Descent<RealType> type_this; //this class;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Size type_size; //non-negative integer;
		typedef Type_Stat type_stat; //state;
		//type: for storing the precision information:
		class type_pre{
			public:
				//values:
				type_real err; //the maximum error;
				type_size n_iter_max; //the maximum iteration number;
				type_stat err_type; //the error types;
				
				//static values:
				static type_stat constexpr val_err_invalid=0; //as if the error is always very large;
				static type_stat constexpr val_err_in_abs=1; //absolute error of input;
				static type_stat constexpr val_err_in_rel=2; //relative error of input;
				static type_stat constexpr val_err_out_abs=3; //absolute error of output;
				static type_stat constexpr val_err_out_rel=4; //relative error of output;
				static type_stat constexpr val_err_slope=5; //(absolute error of output)/(that of input);
				static type_size constexpr val_size_infinity=0;
				
				//Constructor:
				type_pre
					(type_real error_max=LZ_DEF_const_default_precision,
					type_size iteration_max=100,type_stat error_type=type_pre::val_err_out_abs);
				type_pre(type_pre const & precision)=default;
				
				//operator:
				type_pre & operator =(type_pre const & precision)=default;
		}; //class type_pre;
		
		//static values:
		static type_stat constexpr val_stat_null=0; //no result;
		static type_stat constexpr val_stat_good=1; //a successful result;
		static type_stat constexpr val_stat_bad=2; //a doubtful result;
		
		//Constructor:
		Minimisation_Descent();
		Minimisation_Descent(type_this const & solver)=delete;
		Minimisation_Descent(type_this && solver)=delete;
		
		//Destructor:
		~Minimisation_Descent()=default;
		
		//operator:
		type_this & operator =(type_this const & solver)=delete;
		type_this & operator =(type_this && solver)=delete;
		//the following operator () does the same as this->solve(...):
		template<typename FuncType,typename FuncDerType,typename DirType,typename StepType>
		type_this & operator ()
			(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
			DirType & get_direct,StepType & get_step);
		//the following operator () uses numerical derivative:
		template<typename FuncType,typename DirType,typename StepType>
		type_this & operator ()
			(FuncType & func,type_matrix const & sol_init,DirType & get_direct,StepType & get_step);
		
		//Is:
		Type_Bool Is_null() const;
		Type_Bool Is_good() const;
		Type_Bool Is_bad() const;
		
		//Get:
		type_real const & Get_value() const;
		type_matrix const & Get_solution() const;
		
		//Set:
		void Set_precision
			(type_real error_max=LZ_DEF_const_default_precision,
			type_size iteration_max=100,type_stat error_type=type_pre::val_err_out_abs);
		
		//Cal:
		//solve: it solves the minimisation problem:
		//its parameters are:{
		//func: the target function whose parameter list is (type_matrix const &)
		//and whose return type is type_matrix of 1 by 1 or type_real;
		//sol_init: the initial solution;
		//func_der: the derivative of func, whose type is like type_matrix(type_matrix const &);
		//get_direct: get_direct(func,sol_now,gradient) gives
		//the (maybe un-normalised) direction for descending,
		//and it should have a member function Set_begin()
		//if it has to be initialised at the first iteration of a minimisation procedure;
		//get_step: get_step(func,sol_now,gradient,direction,value_as_real(func(sol_now)))
		//gives the step length on this direction,
		//and it should have a member function Set_begin()
		//if it has to be initialised at the first iteration of a minimisation procedure;
		//}:
		template<typename FuncType,typename FuncDerType,typename DirType,typename StepType>
		type_this & solve
			(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
			DirType & get_direct,StepType & get_step);
		
		//static function:
		static type_real value_as_real(type_real const & x);
		static type_real value_as_real(type_matrix const & x);
		static type_matrix value_as_matrix(type_real const & x);
		static type_matrix value_as_matrix(type_matrix const & x);
	protected:
		type_stat stat; //the state;
		type_real min_val; //the function value of the solution;
		type_matrix min_sol; //the solution;
		type_pre pre; //the precision information;
		
		//functions:
		Type_Bool Is_valid_direction(type_matrix const & direct) const;
		Type_Bool Cal_err
			(type_matrix const & x0,type_matrix const & x1,
			type_real const & y0,type_real const & y1) const;
		
		//static functions:
		template<typename Tp>
		static void invoke_Set_begin(Tp & obj);
	private:
		//static functions:
		//invoke_Set_begin(Tp &) calls the following one when Tp has the member function Set_begin():
		template<typename Tp>
		static
		typename std::enable_if<std::is_convertible<
		typename std::add_pointer<decltype(std::declval<Tp>().Set_begin())>::type,void*>::value>::type
		invoke_Set_begin(Tp & obj,int);
		//invoke_Set_begin(Tp &) calls the following one when Tp has no member function Set_begin():
		template<typename Tp>
		static void invoke_Set_begin(Tp & obj,unsigned char);
}; //class Minimisation_Descent;
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Constrained_Descent{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Constrained_Descent<RealType> type_this; //this class;
		typedef Minimisation_Descent<RealType> type_solver; //actual solver for unconstrained minimisation;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Size type_size; //non-negative integer;
		typedef Type_Stat type_stat; //state;
		//type: for storing the precision information:
		class type_pre{
			public:
				//values:
				type_real coef_con_init; //the initial penalty coefficient (>0);
				type_real coef_con_rate; //the increasing rate of the penalty coefficient (>1);
				type_real err_con; //the maximum error for the constraints (>0);
				type_size n_iter_max; //the maximum iteration number for the penalty coefficient;
				
				//static values:
				static type_size constexpr val_size_infinity=0;
				static inline type_real Get_val_coef_penalty_init_auto(){return type_real(-1);}
				
				//Constructor:
				type_pre
					(type_real coef_penalty_init=type_pre::Get_val_coef_penalty_init_auto(),
					type_real coef_penalty_rate=2,
					type_real err_constraint=LZ_DEF_const_default_precision,
					type_size iteration_max=20);
				type_pre(type_pre const & precision)=default;
				
				//operator:
				type_pre & operator =(type_pre const & precision)=default;
		}; //class type_pre;
	protected:
		//type: the augmented function:
		template<typename FuncType,typename FuncDerType,
			typename FuncConEqType,typename FuncConLeqType,
			typename FuncDerConEqType,typename FuncDerConLeqType>
		class type_func_aug;
	public:
		//static values:
		static type_stat constexpr val_stat_null=0; //no result;
		static type_stat constexpr val_stat_good=1; //a successful result;
		static type_stat constexpr val_stat_bad=2; //a doubtful result;
		
		//Constructor:
		Minimisation_Constrained_Descent();
		Minimisation_Constrained_Descent(type_this const & solver)=delete;
		Minimisation_Constrained_Descent(type_this && solver)=delete;
		
		//Destructor:
		~Minimisation_Constrained_Descent()=default;
		
		//operator:
		type_this & operator =(type_this const & solver)=delete;
		type_this & operator =(type_this && solver)=delete;
		
		//Is:
		Type_Bool Is_null() const;
		Type_Bool Is_good() const;
		Type_Bool Is_bad() const;
		
		//Get:
		type_real const & Get_value() const;
		type_matrix const & Get_value_constraint_eq() const;
		type_matrix const & Get_value_constraint_leq() const;
		type_matrix const & Get_solution() const;
		
		//non-const Get:
		type_solver & solver();
		
		//Set:
		void Set_precision
			(type_real coef_penalty_init=type_pre::Get_val_coef_penalty_init_auto(),
			type_real coef_penalty_rate=2,
			type_real err_constraint=LZ_DEF_const_default_precision,
			type_size iteration_max=20);
		
		//Cal:
		//solve: it solves the minimisation problem:
		//its parameters are:{
		//func: the target function whose parameter list is (type_matrix const &)
		//and whose return type is type_matrix of 1 by 1 or type_real;
		//sol_init: the initial solution;
		//func_der: the derivative of func, whose type is like type_matrix(type_matrix const &);
		//get_direct: get_direct(func,sol_now,gradient) gives
		//the (maybe un-normalised) direction for descending,
		//and it should have a member function Set_begin()
		//if it has to be initialised at the first iteration of a minimisation procedure;
		//get_step: get_step(func,sol_now,gradient,direction,value_as_real(func(sol_now)))
		//gives the step length on this direction,
		//and it should have a member function Set_begin()
		//if it has to be initialised at the first iteration of a minimisation procedure;
		//func_con_eq: the constraint function that is component-wise equal to 0,
		//whose parameter list is (type_matrix const &)
		//and whose return type is type_matrix or type_real;
		//func_con_leq: the constraint function that is component-wise less than or equal to 0,
		//whose parameter list is (type_matrix const &)
		//and whose return type is type_matrix or type_real;
		//func_der_con_eq: the derivative of func_con_eq,
		//whose type is like type_matrix(type_matrix const &);
		//func_der_con_leq: the derivative of func_con_leq,
		//whose type is like type_matrix(type_matrix const &);
		//}:
		template<typename FuncType,typename FuncDerType,typename DirType,typename StepType,
			typename FuncConEqType,typename FuncConLeqType,
			typename FuncDerConEqType,typename FuncDerConLeqType>
		type_this & solve
			(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
			DirType & get_direct,StepType & get_step,
			FuncConEqType & func_con_eq,FuncConLeqType & func_con_leq,
			FuncDerConEqType & func_der_con_eq,FuncDerConLeqType & func_der_con_leq);
	protected:
		type_stat stat; //the state;
		type_real min_val; //the function value at the solution;
		type_matrix min_val_con_eq,min_val_con_leq; //the eq and leq-constraint values at the solution;
		type_matrix min_sol; //the solution;
		type_pre pre; //the precision information;
		type_solver min_solver; //the actual solver;
		
		//static function:
		static type_matrix value_as_matrix(type_real const & x);
		static type_matrix value_as_matrix(type_matrix const & x);
}; //class Minimisation_Constrained_Descent;
/****************************************************************************************************/
//}(minimisation algorithm)
/****************************************************************************************************/
//getting descending direction{
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Descent_GetDirect_NegGradient{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Descent_GetDirect_NegGradient<RealType> type_this; //this class;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Stat type_stat; //state;
		
		//static values:
		static type_stat constexpr val_stat_normalised=0; //the direction is normalised;
		static type_stat constexpr val_stat_unnormalised=1; //the direction is not normalised;
		
		//Constructor:
		Minimisation_Descent_GetDirect_NegGradient
			(type_stat tag_normalised=type_this::val_stat_normalised);
		Minimisation_Descent_GetDirect_NegGradient(type_this const &)=delete;
		Minimisation_Descent_GetDirect_NegGradient(type_this &&)=delete;
		
		//operator:
		type_this & operator =(type_this const &)=delete;
		type_this & operator =(type_this &&)=delete;
		template<typename FuncType>
		type_matrix operator ()
			(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient) const;
		
		//Set:
		void Set_normalised(type_stat tag_normalised=type_this::val_stat_normalised);
	protected:
		type_stat tag_norm; //whether the direction is normalised;
}; //class Minimisation_Descent_GetDirect_NegGradient;
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Descent_GetDirect_BFGS{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Descent_GetDirect_BFGS<RealType> type_this; //this class;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Stat type_stat; //state;
		
		//static values:
		static type_stat constexpr val_stat_normalised=0; //the direction is normalised;
		static type_stat constexpr val_stat_unnormalised=1; //the direction is not normalised;
		
		//Constructor:
		Minimisation_Descent_GetDirect_BFGS
			(type_matrix quasiHessianInverse_init=type_matrix(0,0),
			type_stat tag_normalised=type_this::val_stat_normalised);
		Minimisation_Descent_GetDirect_BFGS(type_this const &)=delete;
		Minimisation_Descent_GetDirect_BFGS(type_this &&)=delete;
		
		//operator:
		type_this & operator =(type_this const &)=delete;
		type_this & operator =(type_this &&)=delete;
		template<typename FuncType>
		type_matrix operator ()
			(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient) const;
		
		//Set:
		void Set_begin();
		void Set_begin
			(type_matrix quasiHessianInverse_init,
			type_stat tag_normalised=type_this::val_stat_normalised);
		void Set_normalised(type_stat tag_normalised=type_this::val_stat_normalised);
	protected:
		mutable type_matrix Hinv; //the inverse of the quasi-Hessian matrix;
		mutable type_matrix sol_last,grad_last; //the vectorised solution and gradient at the last iteration;
		mutable Type_Bool tag_init; //tag_init is true if it is at the first iteration;
		type_matrix Hinv_init; //the initial setting of Hinv;
		type_stat tag_norm; //whether the direction is normalised;
}; //class Minimisation_Descent_GetDirect_BFGS;
/****************************************************************************************************/
//}(getting descending direction)
/****************************************************************************************************/
//getting step length{
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Descent_GetStep_GoldenSection{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Descent_GetStep_GoldenSection<RealType> type_this; //this class;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Size type_size; //non-negative integer;
		typedef Type_Stat type_stat; //state;
		//type: for storing the precision information:
		class type_pre{
			public:
				//values:
				type_real step_init; //the initial step length;
				type_size n_iter_max; //the maximum iteration number;
				type_real *pt_step_last; //the last step length;
				
				//Constructor:
				type_pre
					(type_real step_initial=1,type_size iteration_max=10,Type_Bool adaptive=true);
				type_pre(type_pre const & precision);
				
				//Destructor:
				~type_pre();
				
				//operator:
				type_pre & operator =(type_pre const & precision);
		}; //class type_pre;
	public:
		//Constructor:
		Minimisation_Descent_GetStep_GoldenSection();
		Minimisation_Descent_GetStep_GoldenSection(type_this const &)=delete;
		Minimisation_Descent_GetStep_GoldenSection(type_this &&)=delete;
		
		//Destructor:
		~Minimisation_Descent_GetStep_GoldenSection()=default;
		
		//operator:
		type_this & operator =(type_this const &)=delete;
		type_this & operator =(type_this &&)=delete;
		//getting the step length:
		template<typename FuncType>
		type_real operator ()
			(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient,
			type_matrix const & direction,type_real const & val_current) const;
		
		//Set:
		void Set_begin();
		void Set_precision
			(type_real step_initial=1,type_size iteration_max=10,Type_Bool adaptive=true);
	protected:
		type_pre pre; //precision information;
}; //class Minimisation_Descent_GetStep_GoldenSection;
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Minimisation_Descent_GetStep_ArmijoBin{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//types:
		typedef Minimisation_Descent_GetStep_ArmijoBin<RealType> type_this; //this class;
		typedef RealType type_real; //real number;
		typedef Type_MatTemp<type_real> type_matrix; //matrix;
		typedef Type_Size type_size; //non-negative integer;
		typedef Type_Stat type_stat; //state;
		//type: for storing the precision information:
		class type_pre{
			public:
				//values:
				type_real rate_relax; //the rate in (0,0.5);
				type_real step_init; //the initial step length;
				type_size n_iter_max; //the maximum iteration number;
				type_real *pt_step_last; //the last step length;
				
				//Constructor:
				type_pre
					(type_real rate_relaxing=0.25,type_real step_initial=1,type_size iteration_max=10,
					Type_Bool adaptive=true);
				type_pre(type_pre const & precision);
				
				//Destructor:
				~type_pre();
				
				//operator:
				type_pre & operator =(type_pre const & precision);
		}; //class type_pre;
	public:
		//Constructor:
		Minimisation_Descent_GetStep_ArmijoBin();
		Minimisation_Descent_GetStep_ArmijoBin(type_this const &)=delete;
		Minimisation_Descent_GetStep_ArmijoBin(type_this &&)=delete;
		
		//Destructor:
		~Minimisation_Descent_GetStep_ArmijoBin()=default;
		
		//operator:
		type_this & operator =(type_this const &)=delete;
		type_this & operator =(type_this &&)=delete;
		//getting the step length:
		template<typename FuncType>
		type_real operator ()
			(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient,
			type_matrix const & direction,type_real const & val_current) const;
		
		//Set:
		void Set_begin();
		void Set_precision
			(type_real rate_relaxing=0.25,type_real step_initial=1,type_size iteration_max=10,
			Type_Bool adaptive=true);
	protected:
		type_pre pre; //precision information;
}; //class Minimisation_Descent_GetStep_ArmijoBin;
/****************************************************************************************************/
//}(getting step length)
/****************************************************************************************************/
//}(iterative descending algorithm for minimisation without or with constraints)
/****************************************************************************************************/
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/minimisation.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_math_minimisation