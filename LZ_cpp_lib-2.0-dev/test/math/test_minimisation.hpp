#ifndef LZ_DEF_LZ_H_test_math_minimisation
#define LZ_DEF_LZ_H_test_math_minimisation

#include<LZ_H/math/minimisation.hpp>
using namespace Liuze;
using namespace Liuze::math;

#if LZ_DEF_compile_part(0)
int fun_test_Minimisation_Descent(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	class Func_Rosenbrock{
		public:
			Func_Rosenbrock(type_size const & dimension=2):
			dim(dimension>(type_size)0 ? dimension : (type_size)2){}
			type_matrix operator ()(type_matrix const & x) const {
				if(x.rows()!=dim || x.cols()!=(type_size)1) return type_matrix(0,0);
				type_matrix res=type_matrix::Zero(1,1);
				type_real tr0;
				for(type_size i_dim=1;i_dim<dim;++i_dim){
					tr0=x(i_dim-1,0);
					tr0=x(i_dim,0)-tr0*tr0;
					res(0,0)+=(type_real(100)*tr0*tr0);
					tr0=type_real(1)-x(i_dim-1,0);
					res(0,0)+=tr0*tr0;
				}
				return res;
			}
		protected:
			type_size dim;
	};
	type_size dim=5;
	Func_Rosenbrock func_Rosenbrock(dim);
	type_matrix sol_init=type_matrix::Zero(dim,1);
	Minimisation_Descent<type_real> solver;
	Minimisation_Descent_GetDirect_NegGradient<type_real> get_direct_neggrad;
	Minimisation_Descent_GetStep_ArmijoBin<type_real> get_step_armijo;
	solver.Set_precision(1e-5,Minimisation_Descent<type_real>::type_pre::val_size_infinity);
	if(solver(func_Rosenbrock,sol_init,get_direct_neggrad,get_step_armijo).Is_bad()){
		cout<<"bad"<<endl;
	}
	cout<<solver.Get_value()<<endl<<solver.Get_solution().transpose()<<endl;
	return 0;
} //fun: fun_test_Minimisation_Descent;
LZ_DEF_compile_part_setfun(0,fun_test_Minimisation_Descent);
#endif

#if LZ_DEF_compile_part(1)
int fun_test_Minimisation_Descent_constraint_test(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	type_size dim=2;
	type_real coef_con=type_real(1);
	type_real coef_con_rate=2;
	type_real err=1e-5;
	auto func_tar=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(1,1,x.array().exp().sum());
		};
	auto func_con_eq=[](type_matrix const & x)->type_matrix{
			type_real diff=x.array().exp().prod()-type_real(1);
			return type_matrix::Constant(1,1,diff);
		};
	auto func_der_tar=[](type_matrix const & x)->type_matrix{
			return x.array().exp().matrix();
		};
	auto func_der_con_eq=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(x.rows(),x.cols(),x.array().exp().prod());
		};
	Minimisation_Descent<type_real> solver;
	Minimisation_Descent_GetDirect_NegGradient<type_real> get_direct_neggrad;
	Minimisation_Descent_GetStep_ArmijoBin<type_real> get_step_armijo;
	solver.Set_precision(1e-10,Minimisation_Descent<type_real>::type_pre::val_size_infinity,
		Minimisation_Descent<type_real>::type_pre::val_err_in_abs);
	type_matrix sol_init=type_matrix::Constant(dim+1,1,-2.0);
	type_real y_tar=solver.value_as_real(func_tar(sol_init.block(0,0,dim,1)));
	type_real y_con=solver.value_as_real(func_con_eq(sol_init.block(0,0,dim,1)));
	type_size i_iter=0;
	auto func_obj_penalty=[&func_tar,&func_con_eq,&coef_con,&dim,&sol_init]
		(type_matrix const & x)->type_matrix{
			type_matrix val_con_eq=func_con_eq(x);
			return func_tar(x)+coef_con*type_matrix::Constant(1,1,val_con_eq(0,0)*val_con_eq(0,0));
		};
	auto func_der_obj_penalty=[&func_der_tar,&func_con_eq,&func_der_con_eq,&coef_con,&dim,&sol_init]
		(type_matrix const & x)->type_matrix{
			return func_der_tar(x)+type_real(2)*coef_con*func_con_eq(x)(0,0)*func_der_con_eq(x);
		};
	auto func_obj_ALM=[&func_tar,&func_con_eq,&coef_con,&dim,&sol_init]
		(type_matrix const & x)->type_matrix{
			type_matrix val_con_eq=func_con_eq(x);
			return func_tar(x)+sol_init(dim,0)*val_con_eq+
				type_real(0.5)*coef_con*type_matrix::Constant(1,1,val_con_eq(0,0)*val_con_eq(0,0));
		};
	auto func_der_obj_ALM=[&func_der_tar,&func_con_eq,&func_der_con_eq,&coef_con,&dim,&sol_init]
		(type_matrix const & x)->type_matrix{
			return func_der_tar(x)+(sol_init(dim,0)+coef_con*func_con_eq(x)(0,0))*func_der_con_eq(x);
		};
	if(false){
		while(i_iter<200){
			solver(func_obj_ALM,sol_init.block(0,0,dim,1),func_der_obj_ALM,
				get_direct_neggrad,get_step_armijo);
			sol_init.block(0,0,dim,1)=solver.Get_solution();
			y_tar=solver.value_as_real(func_tar(sol_init.block(0,0,dim,1)));
			y_con=solver.value_as_real(func_con_eq(sol_init.block(0,0,dim,1)));
			sol_init(dim,0)+=coef_con*y_con;
			cout<<i_iter<<": "<<solver.Get_value()<<" "<<y_tar<<" "<<y_con<<endl;
			cout<<sol_init.transpose()<<endl;
			if(solver.Is_bad()) cout<<"bad: "<<i_iter<<endl<<endl;
			else if(Liuze::math::abs(y_con)<=1e-5) break;
			++i_iter;
		}
	} else {
		while(i_iter<200){
			solver(func_obj_penalty,sol_init.block(0,0,dim,1),func_der_obj_penalty,
				get_direct_neggrad,get_step_armijo);
			sol_init.block(0,0,dim,1)=solver.Get_solution();
			y_tar=solver.value_as_real(func_tar(sol_init.block(0,0,dim,1)));
			y_con=solver.value_as_real(func_con_eq(sol_init.block(0,0,dim,1)));
			cout<<i_iter<<": "<<solver.Get_value()<<" "<<y_tar<<" "<<y_con<<endl;
			cout<<sol_init.transpose()<<endl;
			if(solver.Is_bad()) cout<<"bad: "<<i_iter<<endl<<endl;
			else if(Liuze::math::abs(y_con)<=1e-5) break;
			++i_iter;
			coef_con*=type_real(2);
		}
	}
	return 0;
} //fun: fun_test_Minimisation_Descent_constraint_test;
LZ_DEF_compile_part_setfun(1,fun_test_Minimisation_Descent_constraint_test);
#endif

#if LZ_DEF_compile_part(2)
int fun_test_Minimisation_Constrained_Descent(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	type_size dim=2;
	/*
	auto func_tar=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(1,1,(x.array()*x.array()).sum());
		};
	auto func_con_eq=[](type_matrix const & x)->type_matrix{
			type_real diff=(x.array()*x.array()).prod()-type_real(1);
			return type_matrix::Constant(1,1,diff*diff);
		};
	auto func_der_tar=[](type_matrix const & x)->type_matrix{
			return type_real(2)*x;
		};
	auto func_der_con_eq=[](type_matrix const & x)->type_matrix{
			type_real p=(x.array()*x.array()).prod();
			if(p==type_real(1)) return type_matrix::Zero(x.rows(),x.cols());
			type_size ir,ic;
			type_matrix res(x.rows(),x.cols());
			for(ir=0;ir<x.rows();++ir){
				for(ic=0;ic<x.cols();++ic){
					res(ir,ic)= x(ir,ic)==type_real(0) ? type_real(0) : p/x(ir,ic);
				}
			}
			return type_real(4)*(p-type_real(1))*res;
		};
	*/
	auto func_tar=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(1,1,x.sum());
		};
	auto func_con_eq=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(1,1,x.prod()-type_real(1));
		};
	auto func_con_leq=[](type_matrix const & x)->type_matrix{
			type_matrix res=-x;
			res(0,0)+=type_real(2);
			return res;
		};
	auto func_der_tar=[](type_matrix const & x)->type_matrix{
			return type_matrix::Constant(x.rows(),x.cols(),type_real(1));
		};
	auto func_der_con_eq=[](type_matrix const & x)->type_matrix{
			type_matrix res=type_matrix::Constant(x.rows(),x.cols(),x.prod());
			type_matrix x0=x;
			type_size ir,ic;
			for(ir=0;ir<x.rows();++ir){
				for(ic=0;ic<x.cols();++ic){
					if(x(ir,ic)!=type_real(0)){
						res(ir,ic)/=x(ir,ic);
					} else {
						x0(ir,ic)=type_real(1);
						res(ir,ic)=x0.prod();
						x0(ir,ic)=x(ir,ic);
					}
				}
			}
			return res;
		};
	auto func_der_con_leq=[](type_matrix const & x)->type_matrix{
			type_size nr_x=x.rows(),nc_x=x.cols();
			type_matrix res=type_matrix::Zero(nr_x*nr_x,nc_x*nc_x);
			type_size ir,ic;
			for(ir=0;ir<x.rows();++ir){
				for(ic=0;ic<x.cols();++ic){
					res(ir*nr_x+ir,ic*nc_x+ic)=type_real(-1);
				}
			}
			return res;
		};
	Minimisation_Constrained_Descent<type_real> solver;
	Minimisation_Descent_GetDirect_NegGradient<type_real> get_direct_neggrad;
	Minimisation_Descent_GetDirect_BFGS<type_real> get_direct_bfgs;
	Minimisation_Descent_GetStep_GoldenSection<type_real> get_step_goldsect;
	Minimisation_Descent_GetStep_ArmijoBin<type_real> get_step_armijo;
	solver.Set_precision(
		Minimisation_Constrained_Descent<type_real>::type_pre::Get_val_coef_penalty_init_auto(),
		2,1e-5,20);
	solver.solver().Set_precision(
		1e-5,
		Minimisation_Descent<type_real>::type_pre::val_size_infinity,
		Minimisation_Descent<type_real>::type_pre::val_err_in_rel);
	type_matrix sol_init;
	auto fun_solve_neggrad_armijo=
		[&solver,&func_tar,&sol_init,&func_der_tar,&get_direct_neggrad,&get_step_armijo,
		&func_con_eq,&func_con_leq,&func_der_con_eq,&func_der_con_leq]
		(void)->void{
			solver.solve(func_tar,sol_init,func_der_tar,get_direct_neggrad,get_step_armijo,
				func_con_eq,func_con_leq,func_der_con_eq,func_der_con_leq);
		};
	auto fun_solve_neggrad_goldsect=
		[&solver,&func_tar,&sol_init,&func_der_tar,&get_direct_neggrad,&get_step_goldsect,
		&func_con_eq,&func_con_leq,&func_der_con_eq,&func_der_con_leq]
		(void)->void{
			solver.solve(func_tar,sol_init,func_der_tar,get_direct_neggrad,get_step_goldsect,
				func_con_eq,func_con_leq,func_der_con_eq,func_der_con_leq);
		};
	auto fun_solve_bfgs_armijo=
		[&solver,&func_tar,&sol_init,&func_der_tar,&get_direct_bfgs,&get_step_armijo,
		&func_con_eq,&func_con_leq,&func_der_con_eq,&func_der_con_leq]
		(void)->void{
			solver.solve(func_tar,sol_init,func_der_tar,get_direct_bfgs,get_step_armijo,
				func_con_eq,func_con_leq,func_der_con_eq,func_der_con_leq);
		};
	auto fun_solve_bfgs_goldsect=
		[&solver,&func_tar,&sol_init,&func_der_tar,&get_direct_bfgs,&get_step_goldsect,
		&func_con_eq,&func_con_leq,&func_der_con_eq,&func_der_con_leq]
		(void)->void{
			solver.solve(func_tar,sol_init,func_der_tar,get_direct_bfgs,get_step_goldsect,
				func_con_eq,func_con_leq,func_der_con_eq,func_der_con_leq);
		};
	Type_ArrayTemp<std::function<void(void)> > fun_solve
		{fun_solve_neggrad_armijo,fun_solve_neggrad_goldsect,
		fun_solve_bfgs_armijo,fun_solve_bfgs_goldsect};
	type_size i_solve;
	for(i_solve=0;i_solve<fun_solve.size();++i_solve){
		cout<<"Solve round "<<i_solve<<":"<<endl;
		sol_init=type_matrix::Constant(1,dim,-1);
		fun_solve[i_solve]();
		if(solver.Is_bad()) cout<<"bad_con"<<endl;
		if(solver.solver().Is_bad()) cout<<"bad_uncon"<<endl;
		cout<<solver.Get_value()<<endl<<endl;
		cout<<solver.Get_solution()<<endl<<endl;
		cout<<solver.Get_value_constraint_eq()<<endl<<endl;
		cout<<solver.Get_value_constraint_leq()<<endl<<endl;
		cout<<solver.solver().Get_value()<<endl<<endl;
		cout<<solver.solver().Get_solution().transpose()<<endl<<endl;
	}
	return 0;
} //fun: fun_test_Minimisation_Constrained_Descent;
LZ_DEF_compile_part_setfun(2,fun_test_Minimisation_Constrained_Descent);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_minimisation