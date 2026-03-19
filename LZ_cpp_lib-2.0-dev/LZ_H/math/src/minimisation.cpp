#if LZ_DEF_LZ_H_math_minimisation!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_minimisation_CPP
#define LZ_DEF_LZ_H_math_src_minimisation_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//class Minimisation_Descent{
//public:
	
	//type: type_pre{
	//public:
		
		template<typename RealType>
		Minimisation_Descent<RealType>::type_pre::type_pre
		(type_real error_max,type_size iteration_max,type_stat error_type):
		err(error_max>type_real(0) ? error_max : type_real(LZ_DEF_const_default_precision)),
		n_iter_max((iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
		iteration_max : type_size(100)),
		err_type(error_type){
		}
		
	//}(type_pre)
	
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Descent<RealType>::Minimisation_Descent
	():
	stat(type_this::val_stat_null),min_val(0),min_sol(0,0),pre(){
	}
	
	//operator:
	
	template<typename RealType>
	template<typename FuncType,typename FuncDerType,typename DirType,typename StepType>
	inline
	typename Minimisation_Descent<RealType>::type_this &
	Minimisation_Descent<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
	DirType & get_direct,StepType & get_step){
		return this->solve(func,sol_init,func_der,get_direct,get_step);
	}
	
	template<typename RealType>
	template<typename FuncType,typename DirType,typename StepType>
	typename Minimisation_Descent<RealType>::type_this &
	Minimisation_Descent<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_init,DirType & get_direct,StepType & get_step){
		auto func_tar=[&func](type_matrix const & x)->type_matrix{
				return type_this::value_as_matrix(func(x));
			};
		Derivative<type_real> func_der(func_tar,(type_size)1,(type_size)1,sol_init.rows(),sol_init.cols());
		return this->solve(func_tar,sol_init,func_der,get_direct,get_step);
	}
	
	//Is:
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Descent<RealType>::Is_null
	() const {
		return this->stat==type_this::val_stat_null;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Descent<RealType>::Is_good
	() const {
		return this->stat==type_this::val_stat_good;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Descent<RealType>::Is_bad
	() const {
		return this->stat==type_this::val_stat_bad;
	}
	
	//Get:
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_real const &
	Minimisation_Descent<RealType>::Get_value
	() const {
		return this->min_val;
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_matrix const &
	Minimisation_Descent<RealType>::Get_solution
	() const {
		return this->min_sol;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent<RealType>::Set_precision
	(type_real error_max,type_size iteration_max,type_stat error_type){
		pre.err= error_max>type_real(0) ? error_max : type_real(LZ_DEF_const_default_precision);
		pre.n_iter_max= (iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
			iteration_max : type_size(100);
		pre.err_type=error_type;
	}
	
	//Cal:
	
	template<typename RealType>
	template<typename FuncType,typename FuncDerType,typename DirType,typename StepType>
	typename Minimisation_Descent<RealType>::type_this &
	Minimisation_Descent<RealType>::solve
	(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
	DirType & get_direct,StepType & get_step){
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<DirType
			(FuncType &,type_matrix const &,type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_real,
			typename std::result_of<StepType
			(FuncType &,type_matrix const &,type_matrix const &,type_matrix const &,type_real const &)
			>::type>::value));
		type_size i_iter;
		type_matrix direct,grad,sol0;
		type_real val0;
		min_sol=sol_init;
		min_val=type_this::value_as_real(func(sol_init));
		type_this::invoke_Set_begin(get_direct);
		type_this::invoke_Set_begin(get_step);
		for(i_iter=0;i_iter<pre.n_iter_max || pre.n_iter_max==type_pre::val_size_infinity;++i_iter){
			sol0=min_sol;
			val0=min_val;
			grad=func_der(sol0);
			if(grad.template lpNorm<Eigen::Infinity>()<=pre.err){
				stat=type_this::val_stat_good;
				return *this;
			}
			direct=get_direct(func,sol0,grad);
			if(!(this->Is_valid_direction(direct))){
				stat=type_this::val_stat_bad;
				return *this;
			}
			min_sol+=(get_step(func,sol0,grad,direct,val0)*direct);
			min_val=type_this::value_as_real(func(min_sol));
			if(this->Cal_err(sol0,min_sol,val0,min_val)){
				stat=type_this::val_stat_good;
				return *this;
			}
		}
		stat=type_this::val_stat_bad;
		return *this;
	}
	
	//static function:
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_real
	Minimisation_Descent<RealType>::value_as_real
	(type_real const & x){
		return x;
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_real
	Minimisation_Descent<RealType>::value_as_real
	(type_matrix const & x){
		return x(0,0);
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_matrix
	Minimisation_Descent<RealType>::value_as_matrix
	(type_real const & x){
		return type_matrix::Constant(1,1,x);
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Descent<RealType>::type_matrix
	Minimisation_Descent<RealType>::value_as_matrix
	(type_matrix const & x){
		return type_matrix::Constant(1,1,x(0,0));
	}
	
//protected:
	
	//functions:
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Descent<RealType>::Is_valid_direction
	(type_matrix const & direct) const {
		return (direct.rows()>(type_size)0 && direct.cols()>(type_size)0 &&
			direct!=type_matrix::Zero(direct.rows(),direct.cols()));
	}
	
	template<typename RealType>
	Type_Bool
	Minimisation_Descent<RealType>::Cal_err
	(type_matrix const & x0,type_matrix const & x1,type_real const & y0,type_real const & y1) const {
		type_size i_dimr,i_dimc;
		type_real diff_y,diff_x;
		switch(pre.err_type){
			case type_pre::val_err_in_abs:
				return (x0-x1).template lpNorm<Eigen::Infinity>()<=pre.err;
			case type_pre::val_err_in_rel:
				for(i_dimr=0;i_dimr<x1.rows();++i_dimr){
					for(i_dimc=0;i_dimc<x1.cols();++i_dimc){
						if(x1(i_dimr,i_dimc)!=type_real(0)){
							if(Liuze::math::abs((x0(i_dimr,i_dimc)-x1(i_dimr,i_dimc))/x1(i_dimr,i_dimc))
							>pre.err) return false;
						} else if(x0(i_dimr,i_dimc)!=type_real(0) && type_real(1)>pre.err){
							return false;
						}
					}
				}
				return true;
			case type_pre::val_err_out_abs:
				return Liuze::math::abs(y0-y1)<=pre.err;
			case type_pre::val_err_out_rel:
				if(y1!=(type_real)0) return Liuze::math::abs((y0-y1)/y1)<=pre.err;
				return y0==(type_real)0 ? true : (type_real)1<=pre.err;
			case type_pre::val_err_slope:
				diff_y=y0-y1;
				for(i_dimr=0;i_dimr<x1.rows();++i_dimr){
					for(i_dimc=0;i_dimc<x1.cols();++i_dimc){
						diff_x=x0(i_dimr,i_dimc)-x1(i_dimr,i_dimc);
						if(diff_x!=type_real(0) && Liuze::math::abs(diff_y/diff_x)>pre.err){
							return false;
						}
					}
				}
				return true;
			default:
				return false;
		}
	}
	
	//static functions:
	
	template<typename RealType>
	template<typename Tp>
	inline
	void
	Minimisation_Descent<RealType>::invoke_Set_begin
	(Tp & obj){
		type_this::invoke_Set_begin(obj,int(0));
	}
	
//private:

	//static functions:
	
	template<typename RealType>
	template<typename Tp>
	inline
	typename std::enable_if<std::is_convertible<
	typename std::add_pointer<decltype(std::declval<Tp>().Set_begin())>::type,void*>::value>::type
	Minimisation_Descent<RealType>::invoke_Set_begin
	(Tp & obj,int){
		obj.Set_begin();
	}
	
	template<typename RealType>
	template<typename Tp>
	inline
	void
	Minimisation_Descent<RealType>::invoke_Set_begin
	(Tp & obj,unsigned char){
	}
	
//}(class Minimisation_Descent)
/****************************************************************************************************/
//class Minimisation_Constrained_Descent{
//public:
	
	//type: type_pre{
	//public:
		
		template<typename RealType>
		Minimisation_Constrained_Descent<RealType>::type_pre::type_pre
		(type_real coef_penalty_init,type_real coef_penalty_rate,type_real err_constraint,
		type_size iteration_max):
		coef_con_init((coef_penalty_init>type_real(0) ||
		coef_penalty_init==type_pre::Get_val_coef_penalty_init_auto()) ?
		coef_penalty_init : type_pre::Get_val_coef_penalty_init_auto()),
		coef_con_rate(coef_penalty_rate>type_real(1) ? coef_penalty_rate : type_real(2)),
		err_con(err_constraint>type_real(0) ? err_constraint : type_real(LZ_DEF_const_default_precision)),
		n_iter_max((iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
		iteration_max : type_size(20)){
		}
		
		
	//}(type_pre)
	
//protected:
	
	//type: type_func_aug{
	
	template<typename RealType>
	template<typename FuncType,typename FuncDerType,
		typename FuncConEqType,typename FuncConLeqType,
		typename FuncDerConEqType,typename FuncDerConLeqType>
	class Minimisation_Constrained_Descent<RealType>::type_func_aug{
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncConEqType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncConEqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerConEqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncConLeqType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncConLeqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerConLeqType(type_matrix const &)>::type>::value));
		protected:
			typedef type_size & type_size_ref;
		public:
			//Constructor:
			type_func_aug
				(FuncType & ref_func,FuncDerType & ref_func_der,
				FuncConEqType & ref_func_con_eq,FuncConLeqType & ref_func_con_leq,
				FuncDerConEqType & ref_func_der_con_eq,FuncDerConLeqType & ref_func_der_con_leq,
				type_real & ref_coef_con,
				type_size & ref_nr_x,type_size & ref_nc_x,
				type_size & ref_nr_con_eq,type_size & ref_nc_con_eq,
				type_size & ref_nr_con_leq,type_size & ref_nc_con_leq,
				type_size & ref_n_x,type_size & ref_n_con_leq,type_size & ref_n_x_aug):
			func(ref_func),func_der(ref_func_der),
			func_con_eq(ref_func_con_eq),func_con_leq(ref_func_con_leq),
			func_der_con_eq(ref_func_der_con_eq),func_der_con_leq(ref_func_der_con_leq),
			coef_con(ref_coef_con),
			nr_x(ref_nr_x),nc_x(ref_nc_x),
			nr_con_eq(ref_nr_con_eq),nc_con_eq(ref_nc_con_eq),
			nr_con_leq(ref_nr_con_leq),nc_con_leq(ref_nc_con_leq),
			n_x(ref_n_x),n_con_leq(ref_n_con_leq),n_x_aug(ref_n_x_aug){
			}
			
			//augmented function (from R^(n_x_aug cross 1) to R):
			type_real func_aug(type_matrix const & x_aug){
				x.resize(nr_x,nc_x); //x is the non-augmented solution;
				tran_from_aug(x_aug,x);
				return type_solver::value_as_real(func(x))+coef_con*
					(func_penalty_con_eq_sum(type_this::value_as_matrix(func_con_eq(x)))+
					func_penalty_con_leq_sum(type_this::value_as_matrix(func_con_leq(x)),x_aug));
			}
			
			//the derivative of the augmented function (from R^(n_x_aug cross 1) to R^(n_x_aug cross 1)):
			type_matrix func_der_aug(type_matrix const & x_aug){
				type_real val_penalty;
				x.resize(nr_x,nc_x); //x is the non-augmented solution;
				tran_from_aug(x_aug,x);
				type_matrix res_aug(n_x_aug,1);
				type_matrix res=type_matrix::Zero(nr_x,nc_x);
				type_matrix val_con=type_this::value_as_matrix(func_con_leq(x)); //the leq-constraint;
				type_matrix val_der_con=func_der_con_leq(x); //the derivative of the leq-constraint;
				for(ic=0;ic<nc_con_leq;++ic){
					for(ir=0;ir<nr_con_leq;++ir){
						val_penalty=func_penalty_der_con_leq
							(val_con(ir,ic),func_penalty_con_leq_slack(x_aug,ir,ic));
						res+=(val_penalty*val_der_con.block(ir*nr_x,ic*nc_x,nr_x,nc_x));
						Get_aux_con_leq(res_aug,ir,ic)=
							(val_penalty*func_penalty_der_con_leq_slack(x_aug,ir,ic))*coef_con;
					}
				}
				val_con=type_this::value_as_matrix(func_con_eq(x)); //the eq-constraint;
				val_der_con=func_der_con_eq(x); //the derivative of the eq-constraint;
				for(ic=0;ic<nc_con_eq;++ic){
					for(ir=0;ir<nr_con_eq;++ir){
						res+=(func_penalty_der_con_eq(val_con(ir,ic))*
							val_der_con.block(ir*nr_x,ic*nc_x,nr_x,nc_x));
					}
				}
				res*=coef_con;
				res+=func_der(x);
				tran_to_aug(res,res_aug);
				return res_aug;
			}
			
			//transforming solution to augmented type:
			void tran_to_aug(type_matrix const & x,type_matrix & x_aug){
				for(ic=0;ic<nc_x;++ic) x_aug.block(ic*nr_x,0,nr_x,1)=x.col(ic);
			}
			//transforming solution to original type:
			void tran_from_aug(type_matrix const & x_aug,type_matrix & x){
				for(ic=0;ic<nc_x;++ic) x.col(ic)=x_aug.block(ic*nr_x,0,nr_x,1);
			}
			//*
			//initialising:
			void Set_sol_init_aug(type_matrix const & sol_init,type_matrix & sol_init_aug){
				tran_to_aug(sol_init,sol_init_aug);
				sol_init_aug.block(n_x,0,n_con_leq,1).setZero();
			}
			
			type_real & Get_aux_con_leq
			(type_matrix & x_aug,type_size const & i_row,type_size const & i_col){
				return x_aug(n_x+i_col*nr_con_leq+i_row,0);
			}
			type_real const & Get_aux_con_leq
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return x_aug(n_x+i_col*nr_con_leq+i_row,0);
			}
			
			//the penalty function for eq-constraint:
			type_real func_penalty_con_eq(type_real const & val_con){
				return val_con*val_con;
			}
			type_real func_penalty_der_con_eq(type_real const & val_con){
				return type_real(2)*val_con;
			}
			type_real func_penalty_con_eq_sum(type_matrix const & val_con){
				type_real res=0;
				for(ic=0;ic<nc_con_eq;++ic){
					for(ir=0;ir<nr_con_eq;++ir){
						res+=func_penalty_con_eq(val_con(ir,ic));
					}
				}
				return res;
			}
			//the penalty function for leq-constraint:
			type_real func_penalty_con_leq_slack
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return type_real(exp(Get_aux_con_leq(x_aug,i_row,i_col)));
			}
			type_real func_penalty_der_con_leq_slack
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return type_real(exp(Get_aux_con_leq(x_aug,i_row,i_col)));
			}
			type_real func_penalty_con_leq
			(type_real const & val_con,type_real const & val_slack){
				return func_penalty_con_eq(val_con+val_slack);
			}
			type_real func_penalty_der_con_leq
			(type_real const & val_con,type_real const & val_slack){
				return func_penalty_der_con_eq(val_con+val_slack);
			}
			type_real func_penalty_con_leq_sum
			(type_matrix const & val_con,type_matrix const & x_aug){
				type_real res=0;
				for(ic=0;ic<nc_con_eq;++ic){
					for(ir=0;ir<nr_con_eq;++ir){
						res+=func_penalty_con_leq(val_con(ir,ic),func_penalty_con_leq_slack(x_aug,ir,ic));
					}
				}
				return res;
			}
			//*/
			/*
			//initialising:
			void Set_sol_init_aug(type_matrix const & sol_init,type_matrix & sol_init_aug){
				tran_to_aug(sol_init,sol_init_aug);
			}
			
			type_real & Get_aux_con_leq
			(type_matrix & x_aug,type_size const & i_row,type_size const & i_col){
				return x_aug(0,0);
			}
			type_real const & Get_aux_con_leq
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return x_aug(0,0);
			}
			
			//the penalty function for eq-constraint:
			type_real func_penalty_con_eq(type_real const & val_con){
				return val_con*val_con;
			}
			type_real func_penalty_der_con_eq(type_real const & val_con){
				return type_real(2)*val_con;
			}
			type_real func_penalty_con_eq_sum(type_matrix const & val_con){
				type_real res=0;
				for(ic=0;ic<nc_con_eq;++ic){
					for(ir=0;ir<nr_con_eq;++ir){
						res+=func_penalty_con_eq(val_con(ir,ic));
					}
				}
				return res;
			}
			//the penalty function for leq-constraint:
			type_real func_penalty_con_leq_slack
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return type_real(0);
			}
			type_real func_penalty_der_con_leq_slack
			(type_matrix const & x_aug,type_size const & i_row,type_size const & i_col){
				return type_real(0);
			}
			type_real func_penalty_con_leq
			(type_real const & val_con,type_real const & val_slack){
				return val_con<=type_real(0) ? type_real(0) : val_con*val_con;
			}
			type_real func_penalty_der_con_leq
			(type_real const & val_con,type_real const & val_slack){
				return val_con<=type_real(0) ? type_real(0) : type_real(2)*val_con;
			}
			type_real func_penalty_con_leq_sum
			(type_matrix const & val_con,type_matrix const & x_aug){
				type_real res=0;
				for(ic=0;ic<nc_con_eq;++ic){
					for(ir=0;ir<nr_con_eq;++ir){
						res+=func_penalty_con_leq(val_con(ir,ic),func_penalty_con_leq_slack(x_aug,ir,ic));
					}
				}
				return res;
			}
			//*/
		protected:
			FuncType & func;
			FuncDerType & func_der;
			FuncConEqType & func_con_eq;
			FuncConLeqType & func_con_leq;
			FuncDerConEqType & func_der_con_eq;
			FuncDerConLeqType & func_der_con_leq;
			type_real & coef_con;
			type_size_ref nr_x,nc_x,nr_con_eq,nc_con_eq,nr_con_leq,nc_con_leq,n_x,n_con_leq,n_x_aug;
			
			type_matrix x;
			type_size ir,ic,i;
	};
	
	//}(type_func_aug)
	
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Constrained_Descent<RealType>::Minimisation_Constrained_Descent
	():
	stat(type_this::val_stat_null),min_val(0),min_val_con_eq(),min_val_con_leq(),min_sol(),pre(),
	min_solver(){
	}
	
	//Is:
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Constrained_Descent<RealType>::Is_null
	() const {
		return this->stat==type_this::val_stat_null;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Constrained_Descent<RealType>::Is_good
	() const {
		return this->stat==type_this::val_stat_good;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Minimisation_Constrained_Descent<RealType>::Is_bad
	() const {
		return this->stat==type_this::val_stat_bad;
	}
	
	//Get:
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_real const &
	Minimisation_Constrained_Descent<RealType>::Get_value
	() const {
		return this->min_val;
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_matrix const &
	Minimisation_Constrained_Descent<RealType>::Get_value_constraint_eq
	() const {
		return this->min_val_con_eq;
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_matrix const &
	Minimisation_Constrained_Descent<RealType>::Get_value_constraint_leq
	() const {
		return this->min_val_con_leq;
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_matrix const &
	Minimisation_Constrained_Descent<RealType>::Get_solution
	() const {
		return this->min_sol;
	}
	
	//non-const Get:
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_solver &
	Minimisation_Constrained_Descent<RealType>::solver
	(){
		return this->min_solver;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Constrained_Descent<RealType>::Set_precision
	(type_real coef_penalty_init,type_real coef_penalty_rate,type_real err_constraint,
	type_size iteration_max){
		pre.coef_con_init= (coef_penalty_init>type_real(0) ||
			coef_penalty_init==type_pre::Get_val_coef_penalty_init_auto()) ?
			coef_penalty_init : type_pre::Get_val_coef_penalty_init_auto();
		pre.coef_con_rate= coef_penalty_rate>type_real(1) ? coef_penalty_rate : type_real(2);
		pre.err_con= err_constraint>type_real(0) ? err_constraint :
			type_real(LZ_DEF_const_default_precision);
		pre.n_iter_max= (iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
			iteration_max : type_size(20);
	}
	
	//Cal:
	
	template<typename RealType>
	template<typename FuncType,typename FuncDerType,typename DirType,typename StepType,
		typename FuncConEqType,typename FuncConLeqType,
		typename FuncDerConEqType,typename FuncDerConLeqType>
	typename Minimisation_Constrained_Descent<RealType>::type_this &
	Minimisation_Constrained_Descent<RealType>::solve
	(FuncType & func,type_matrix const & sol_init,FuncDerType & func_der,
	DirType & get_direct,StepType & get_step,
	FuncConEqType & func_con_eq,FuncConLeqType & func_con_leq,
	FuncDerConEqType & func_der_con_eq,FuncDerConLeqType & func_der_con_leq){
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<DirType
			(FuncType &,type_matrix const &,type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_real,
			typename std::result_of<StepType
			(FuncType &,type_matrix const &,type_matrix const &,type_matrix const &,type_real const &)
			>::type>::value));
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncConEqType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncConEqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerConEqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncConLeqType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncConLeqType(type_matrix const &)>::type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::result_of<FuncDerConLeqType(type_matrix const &)>::type>::value));
		typedef
			type_func_aug<FuncType,FuncDerType,FuncConEqType,FuncConLeqType,
			FuncDerConEqType,FuncDerConLeqType>
			type_func_augment;
		min_sol=sol_init;
		min_val_con_eq=type_this::value_as_matrix(func_con_eq(sol_init));
		min_val_con_leq=type_this::value_as_matrix(func_con_leq(sol_init));
		type_size nr_x=sol_init.rows(),nc_x=sol_init.cols(); //dim. of solution;
		type_size nr_con_eq=min_val_con_eq.rows(),nc_con_eq=min_val_con_eq.cols(); //dim. of eq-constraint;
		type_size nr_con_leq=min_val_con_leq.rows(),nc_con_leq=min_val_con_leq.cols(); //dim. of leq-constraint;
		type_size n_x=nr_x*nc_x,n_con_leq=nr_con_leq*nc_con_leq;
		type_size n_x_aug=n_x+n_con_leq; //the dim. of augmented solution;
		type_real coef_con=pre.coef_con_init; //the coefficient of the penalty function;
		//the augmented function:
		type_func_augment
			func_augment(func,func_der,func_con_eq,func_con_leq,func_der_con_eq,func_der_con_leq,
			coef_con,nr_x,nc_x,nr_con_eq,nc_con_eq,nr_con_leq,nc_con_leq,n_x,n_con_leq,n_x_aug);
		type_matrix sol_init_aug(n_x_aug,1);
		func_augment.Set_sol_init_aug(sol_init,sol_init_aug); //set the initial augmented solution;
		//deal with auto initial penalty coefficient:
		if(pre.coef_con_init==type_pre::Get_val_coef_penalty_init_auto()){
			type_real scale_con=
				func_augment.func_penalty_con_eq_sum(min_val_con_eq)+
				func_augment.func_penalty_con_leq_sum(min_val_con_leq,sol_init_aug);
			if(scale_con<=pre.err_con){
				coef_con=type_real(0.1);
			} else {
				type_real scale_func=Liuze::math::abs(type_solver::value_as_real(func(sol_init)));
				coef_con= scale_func>type_real(0) ? scale_func/scale_con*type_real(0.1) : type_real(0.1);
			}
		}
		//the augmented target function:
		std::function<type_real(type_matrix const &)> func_aug=
			std::bind(&type_func_augment::func_aug,std::ref(func_augment),std::placeholders::_1);
		//the derivative of augmented target function:
		std::function<type_matrix(type_matrix const &)> func_der_aug=
			std::bind(&type_func_augment::func_der_aug,std::ref(func_augment),std::placeholders::_1);
		type_size i_iter;
		for(i_iter=0;i_iter<pre.n_iter_max || pre.n_iter_max==type_pre::val_size_infinity;++i_iter){
			min_solver.solve(func_aug,sol_init_aug,func_der_aug,get_direct,get_step);
			func_augment.tran_from_aug(min_solver.Get_solution(),min_sol);
			min_val_con_eq=type_this::value_as_matrix(func_con_eq(min_sol));
			min_val_con_leq=type_this::value_as_matrix(func_con_leq(min_sol));
			if((min_val_con_eq.array().abs()<=pre.err_con).all() &&
			(min_val_con_leq.array()<=pre.err_con).all() && min_solver.Is_good()){
				stat=type_this::val_stat_good;
				min_val=type_solver::value_as_real(func(min_sol));
				return *this;
			}
			coef_con*=pre.coef_con_rate;
		}
		stat=type_this::val_stat_bad;
		min_val=type_solver::value_as_real(func(min_sol));
		return *this;
	}
	
//protected:
	
	//static function:
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_matrix
	Minimisation_Constrained_Descent<RealType>::value_as_matrix
	(type_real const & x){
		return type_matrix::Constant(1,1,x);
	}
	
	template<typename RealType>
	inline
	typename Minimisation_Constrained_Descent<RealType>::type_matrix
	Minimisation_Constrained_Descent<RealType>::value_as_matrix
	(type_matrix const & x){
		return x;
	}
	
//}(class Minimisation_Constrained_Descent)
/****************************************************************************************************/
//class Minimisation_Descent_GetDirect_NegGradient{
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Descent_GetDirect_NegGradient<RealType>::Minimisation_Descent_GetDirect_NegGradient
	(type_stat tag_normalised):
	tag_norm(tag_normalised==type_this::val_stat_normalised ?
	type_this::val_stat_normalised : type_this::val_stat_unnormalised){
	}
	
	//operator:
	
	template<typename RealType>
	template<typename FuncType>
	inline
	typename Minimisation_Descent_GetDirect_NegGradient<RealType>::type_matrix
	Minimisation_Descent_GetDirect_NegGradient<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient) const {
		if(tag_norm==type_this::val_stat_unnormalised) return -gradient;
		if(gradient.size()==(Type_Size)0) return type_matrix(0,0);
		type_real len=gradient.template lpNorm<2>();
		return len>(type_real)0 ? -gradient/len : type_matrix(0,0);
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetDirect_NegGradient<RealType>::Set_normalised
	(type_stat tag_normalised){
		tag_norm= tag_normalised==type_this::val_stat_normalised ?
			type_this::val_stat_normalised : type_this::val_stat_unnormalised;
	}
	
//}(class Minimisation_Descent_GetDirect_NegGradient)
/****************************************************************************************************/
//class Minimisation_Descent_GetDirect_BFGS{
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Descent_GetDirect_BFGS<RealType>::Minimisation_Descent_GetDirect_BFGS
	(type_matrix quasiHessianInverse_init,type_stat tag_normalised):
	Hinv(),sol_last(),grad_last(),tag_init(true),Hinv_init(quasiHessianInverse_init),
	tag_norm(tag_normalised==type_this::val_stat_normalised ?
	type_this::val_stat_normalised : type_this::val_stat_unnormalised){
	}
	
	//operator:
	
	template<typename RealType>
	template<typename FuncType>
	typename Minimisation_Descent_GetDirect_BFGS<RealType>::type_matrix
	Minimisation_Descent_GetDirect_BFGS<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient) const {
		Type_Size dim=gradient.size();
		if(dim==Type_Size(0)) return type_matrix(0,0);
		type_matrix sol(dim,1),grad(dim,1); //the vectorised sol_current and gradient;
		if(gradient.cols()==Type_Size(1)){
			sol=sol_current;
			grad=gradient;
		} else {
			for(Type_Size i=0;i<gradient.cols();++i){
				sol.block(i*gradient.rows(),0,gradient.rows(),1)=sol_current.col(i);
				grad.block(i*gradient.rows(),0,gradient.rows(),1)=gradient.col(i);
			}
		}
		//renew Hinv:
		if(tag_init){
			if(Hinv_init.rows()!=dim || Hinv_init.cols()!=dim) Hinv=type_matrix::Identity(dim,dim);
			else Hinv=Hinv_init;
			tag_init=false;
		} else {
			type_matrix dsol=sol-sol_last,dgrad=grad-grad_last;
			type_real prod=(dsol.transpose()*dgrad)(0,0);
			if(prod>type_real(LZ_DEF_const_default_precision)){
				dgrad=type_matrix::Identity(dim,dim)-dgrad*dsol.transpose()/prod;
				Hinv=dgrad.transpose()*Hinv*dgrad+dsol*dsol.transpose()/prod;
			}
		}
		sol_last=sol;
		grad_last=grad;
		grad=-Hinv*grad;
		if(tag_norm==type_this::val_stat_normalised){
			type_real len=grad.template lpNorm<2>();
			if(len<=type_real(0)) return type_matrix(0,0);
			grad/=len;
		}
		if(gradient.cols()==Type_Size(1)) return grad;
		sol.resize(gradient.rows(),gradient.cols());
		for(Type_Size i=0;i<gradient.cols();++i){
			sol.col(i)=grad.block(i*gradient.rows(),0,gradient.rows(),1);
		}
		return sol;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetDirect_BFGS<RealType>::Set_begin
	(){
		tag_init=true;
	}
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetDirect_BFGS<RealType>::Set_begin
	(type_matrix quasiHessianInverse_init,type_stat tag_normalised){
		Hinv_init=quasiHessianInverse_init;
		tag_init=true;
		tag_norm= tag_normalised==type_this::val_stat_normalised ?
			type_this::val_stat_normalised : type_this::val_stat_unnormalised;
	}
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetDirect_BFGS<RealType>::Set_normalised
	(type_stat tag_normalised){
		tag_norm= tag_normalised==type_this::val_stat_normalised ?
			type_this::val_stat_normalised : type_this::val_stat_unnormalised;
	}
	
//}(class Minimisation_Descent_GetDirect_BFGS)
/****************************************************************************************************/
//class Minimisation_Descent_GetStep_GoldenSection{
//public:
	
	//type: type_pre{
	//public:
		
		//Constructor:
		
		template<typename RealType>
		Minimisation_Descent_GetStep_GoldenSection<RealType>::type_pre::type_pre
		(type_real step_initial,type_size iteration_max,Type_Bool adaptive):
		step_init(step_initial>type_real(0) ? step_initial : type_real(1)),
		n_iter_max(iteration_max>type_size(0) ? iteration_max : type_size(10)),
		pt_step_last(adaptive ? new type_real : NULL){
			if(pt_step_last) *pt_step_last=step_init;
		}
		
		template<typename RealType>
		Minimisation_Descent_GetStep_GoldenSection<RealType>::type_pre::type_pre
		(type_pre const & precision):
		step_init(precision.step_init),n_iter_max(precision.n_iter_max),
		pt_step_last(precision.pt_step_last ? new type_real(*(precision.pt_step_last)) : NULL){
		}
		
		//Destructor:
		
		template<typename RealType>
		Minimisation_Descent_GetStep_GoldenSection<RealType>::type_pre::~type_pre
		(){
			if(pt_step_last) delete pt_step_last;
		}
		
		//operator:
		
		template<typename RealType>
		typename Minimisation_Descent_GetStep_GoldenSection<RealType>::type_pre &
		Minimisation_Descent_GetStep_GoldenSection<RealType>::type_pre::operator =
		(type_pre const & precision){
			step_init=precision.step_init;
			n_iter_max=precision.n_iter_max;
			if(pt_step_last){
				if(precision.pt_step_last){
					*pt_step_last=*(precision.pt_step_last);
				} else {
					delete pt_step_last;
					pt_step_last=NULL;
				}
			} else if(precision.pt_step_last){
				pt_step_last= new type_real(*(precision.pt_step_last));
			}
			return *this;
		}
		
	//}(type_pre)
	
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Descent_GetStep_GoldenSection<RealType>::Minimisation_Descent_GetStep_GoldenSection
	():
	pre(){
	}
	
	//operator:
	
	template<typename RealType>
	template<typename FuncType>
	typename Minimisation_Descent_GetStep_GoldenSection<RealType>::type_real
	Minimisation_Descent_GetStep_GoldenSection<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient,
	type_matrix const & direction,type_real const & val_current) const {
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value));
		type_size i_iter=0;
		type_real step_best,val_best;
		type_real step[4],val[4];
		step[0]=type_real(0);
		val[0]=val_current;
		step[1]=(pre.pt_step_last ? *(pre.pt_step_last) : pre.step_init);
		val[1]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[1]*direction));
		step[3]=type_real(-1); //step[3] is invalid;
		//such shape for function values from step[0] to step[1]: "/":
		while(val[0]<val[1]){
			step[3]=step[1]; //step[3] is made valid;
			val[3]=val[1];
			step[1]=type_real(0.5)*step[3];
			if(step[1]<=type_real(LZ_DEF_const_default_precision)) return type_real(0);
			val[1]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[1]*direction));
			++i_iter;
		}
		//now such shape for function values from step[0] to step[1]: "\".
		if(i_iter>=pre.n_iter_max){
			if(pre.pt_step_last) *(pre.pt_step_last)=step[1];
			return step[1];
		}
		step_best=step[1];
		val_best=val[1];
		//such shape for function values at step[0], step[1] and step[3]: "\/":
		if(step[0]<step[3]) goto GTS_unimode;
		step[2]=step[1]+step[1];
		val[2]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[2]*direction));
		//such shape for function values from step[0] to step[2]: "\/":
		if(val[1]<val[2]){
			step[3]=step[2];
			val[3]=val[2];
			goto GTS_unimode;
		}
		//now such shape for function values from step[0] to step[2]: "\\".
		while(i_iter<pre.n_iter_max){
			step[3]=step[2]+step[2];
			val[3]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[3]*direction));
			//such shape for function values from step[0] to step[3]: "\\/":
			if(val[2]<val[3]){
				step_best=step[2];
				val_best=val[2];
				step[0]=step[1];
				val[0]=val[1];
				goto GTS_unimode;
			}
			//now such shape for function values from step[0] to step[3]: "\\\".
			step[1]=step[2];
			val[1]=val[2];
			step[2]=step[3];
			val[2]=val[3];
			++i_iter;
		}
		//iter==pre.n_iter_max and now such shape for function values from step[0] to step[3]: "\\\".
		if(pre.pt_step_last) *(pre.pt_step_last)=step[3];
		return step[3];
		//now the unimode interval has been found as [step[0],step[3]]:
		GTS_unimode:
		step[1]=type_real(0.618)*step[0]+type_real(0.382)*step[3];
		step[2]=step[3]-(step[1]-step[0]);
		val[1]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[1]*direction));
		val[2]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[2]*direction));
		//now using the golden section method to find the best step length in [step[0],step[3]]:
		while(i_iter<pre.n_iter_max){
			if(val[1]<val[2]){
				step[3]=step[2];
				val[3]=val[2];
				step[2]=step[1];
				val[2]=val[1];
				step[1]=step[0]+(step[3]-step[2]);
				val[1]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[1]*direction));
			} else {
				step[0]=step[1];
				val[0]=val[1];
				step[1]=step[2];
				val[1]=val[2];
				step[2]=step[3]-(step[1]-step[0]);
				val[2]=Minimisation_Descent<type_real>::value_as_real(func(sol_current+step[2]*direction));
			}
			++i_iter;
		}
		i_iter= val[1]<val[2] ? 1 : 2; //step[i_iter] is the best step length obtained by the above search;
		if(val[i_iter]<=val_best) step_best=step[i_iter];
		if(pre.pt_step_last) *(pre.pt_step_last)=step_best;
		return step_best;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetStep_GoldenSection<RealType>::Set_begin
	(){
		if(pre.pt_step_last) *(pre.pt_step_last)=pre.step_init;
	}
	
	template<typename RealType>
	void
	Minimisation_Descent_GetStep_GoldenSection<RealType>::Set_precision
	(type_real step_initial,type_size iteration_max,Type_Bool adaptive){
		pre.step_init= step_initial>type_real(0) ? step_initial : type_real(1);
		pre.n_iter_max= iteration_max>type_size(0) ? iteration_max : type_size(10);
		if(adaptive){
			if(pre.pt_step_last==NULL) pre.pt_step_last=new type_real(pre.step_init);
		} else if(pre.pt_step_last){
			delete pre.pt_step_last;
			pre.pt_step_last=NULL;
		}
	}
	
//}(class Minimisation_Descent_GetStep_GoldenSection)
/****************************************************************************************************/
//class Minimisation_Descent_GetStep_ArmijoBin{
//public:
	
	//type: type_pre{
	//public:
		
		//Constructor:
		
		template<typename RealType>
		Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_pre::type_pre
		(type_real rate_relaxing,type_real step_initial,type_size iteration_max,Type_Bool adaptive):
		rate_relax((rate_relaxing>type_real(0) && rate_relaxing<type_real(0.5)) ?
		rate_relaxing : type_real(0.25)),
		step_init(step_initial>type_real(0) ? step_initial : type_real(1)),
		n_iter_max(iteration_max>type_size(0) ? iteration_max : type_size(10)),
		pt_step_last(adaptive ? new type_real : NULL){
			if(pt_step_last) *pt_step_last=step_init;
		}
		
		template<typename RealType>
		Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_pre::type_pre
		(type_pre const & precision):
		rate_relax(precision.rate_relax),step_init(precision.step_init),
		n_iter_max(precision.n_iter_max),
		pt_step_last(precision.pt_step_last ? new type_real(*(precision.pt_step_last)) : NULL){
		}
		
		//Destructor:
		
		template<typename RealType>
		Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_pre::~type_pre
		(){
			if(pt_step_last) delete pt_step_last;
		}
		
		//operator:
		
		template<typename RealType>
		typename Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_pre &
		Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_pre::operator =
		(type_pre const & precision){
			rate_relax=precision.rate_relax;
			step_init=precision.step_init;
			n_iter_max=precision.n_iter_max;
			if(pt_step_last){
				if(precision.pt_step_last){
					*pt_step_last=*(precision.pt_step_last);
				} else {
					delete pt_step_last;
					pt_step_last=NULL;
				}
			} else if(precision.pt_step_last){
				pt_step_last= new type_real(*(precision.pt_step_last));
			}
			return *this;
		}
		
	//}(type_pre)
	
//public:
	
	//Constructor:
	
	template<typename RealType>
	Minimisation_Descent_GetStep_ArmijoBin<RealType>::Minimisation_Descent_GetStep_ArmijoBin
	():
	pre(){
	}
	
	//operator:
	
	template<typename RealType>
	template<typename FuncType>
	typename Minimisation_Descent_GetStep_ArmijoBin<RealType>::type_real
	Minimisation_Descent_GetStep_ArmijoBin<RealType>::operator ()
	(FuncType & func,type_matrix const & sol_current,type_matrix const & gradient,
	type_matrix const & direction,type_real const & val_current) const {
		LZ_DEF_func_check_traits
			((std::is_same<type_real,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value ||
			std::is_same<type_matrix,
			typename std::result_of<FuncType(type_matrix const &)>::type>::value));
		type_real const rate_more=10,rate_less=0.5;
		type_real val_unit=(gradient.array()*direction.array()).sum()*pre.rate_relax;
		type_real step=(pre.pt_step_last ? *(pre.pt_step_last) : pre.step_init),step0=0;
		type_real step_diff=step;
		type_size i_iter=0;
		Type_Bool tag_step0=false; //tag_step0 is true if step0 satisfies the Armijo condition;
		//making the needed step length in [step0,step]:
		while(i_iter<pre.n_iter_max &&
		Minimisation_Descent<type_real>::value_as_real(func(sol_current+step*direction))
		<=val_current+step*val_unit){
			tag_step0=true;
			step0=step;
			step_diff*=rate_more;
			step+=step_diff;
			++i_iter;
		}
		step*=rate_less;
		step_diff=step;
		while(i_iter<pre.n_iter_max ||
		(!tag_step0 && step_diff>type_real(LZ_DEF_const_default_precision))){
			step_diff*=rate_less;
			if(Minimisation_Descent<type_real>::value_as_real(func(sol_current+step*direction))
			<=val_current+step*val_unit){
				tag_step0=true;
				if(step0<step) step0=step;
				step+=step_diff;
			} else {
				step-=step_diff;
			}
			++i_iter;
		}
		step=(Minimisation_Descent<type_real>::value_as_real(func(sol_current+step*direction))
			<=val_current+step*val_unit) ? step : step0;
		if(pre.pt_step_last && step>type_real(0)) *(pre.pt_step_last)=step;
		return step;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Minimisation_Descent_GetStep_ArmijoBin<RealType>::Set_begin
	(){
		if(pre.pt_step_last) *(pre.pt_step_last)=pre.step_init;
	}
	
	template<typename RealType>
	void
	Minimisation_Descent_GetStep_ArmijoBin<RealType>::Set_precision
	(type_real rate_relaxing,type_real step_initial,type_size iteration_max,Type_Bool adaptive){
		pre.rate_relax= (rate_relaxing>type_real(0) && rate_relaxing<type_real(0.5)) ?
			rate_relaxing : type_real(0.25);
		pre.step_init= step_initial>type_real(0) ? step_initial : type_real(1);
		pre.n_iter_max= iteration_max>type_size(0) ? iteration_max : type_size(10);
		if(adaptive){
			if(pre.pt_step_last==NULL) pre.pt_step_last=new type_real(pre.step_init);
		} else if(pre.pt_step_last){
			delete pre.pt_step_last;
			pre.pt_step_last=NULL;
		}
	}
	
//}(class Minimisation_Descent_GetStep_ArmijoBin)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_minimisation_CPP
#endif //#if LZ_DEF_LZ_H_math_minimisation!=202105L