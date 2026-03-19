#ifndef LZ_DEF_LZ_H_test_math_differentiation
#define LZ_DEF_LZ_H_test_math_differentiation

#include<LZ_H/math/differentiation.hpp>
using namespace Liuze;
using namespace Liuze::math;

#if LZ_DEF_compile_part(0)
int fun_test_Derivative(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	type_size dim;
	cout<<"dim=";
	cin>>dim;
	if(dim<=(type_size)0) dim=1;
	Type_ArrayTemp<type_matrix> coef(3);
	type_matrix x(dim,1);
	Derivative<type_real> der1,der2;
	while(true){
		coef[2]=type_matrix::Random(dim,dim);
		coef[1]=type_matrix::Random(dim,1);
		coef[0]=type_matrix::Random(1,1);
		auto f0=[&coef](type_matrix const & x)->type_matrix{
				return x.transpose()*coef[2]*x+coef[1].transpose()*x+coef[0];
			};
		auto f1=[&coef](type_matrix const & x)->type_matrix{
				return (coef[2]+coef[2].transpose())*x+coef[1];
			};
		auto f2=[&coef](type_matrix const & x)->type_matrix{
				return coef[2]+coef[2].transpose();
			};
		der1.Set_function(f0);
		auto der1tran=[&der1](type_matrix const & x)->type_matrix{
				return der1(x).transpose();
			};
		der2.Set_function(der1tran);
		while(true){
			x=type_matrix::Random(dim,1);
			cout<<(f1(x)-der1(x)).transpose()<<endl<<endl;
			cout<<(f2(x)-der2(x))<<endl<<endl;
			cout<<"continue=";
			cin>>x(0,0);
			if(x(0,0)<=type_real(1e-3)) break;
		}
		der2.Set_function();
		der1.Set_function();
		cout<<"continue=";
		cin>>x(0,0);
		if(x(0,0)<=type_real(1e-3)) break;
	}
	return 0;
} //fun: fun_test_Derivative;
LZ_DEF_compile_part_setfun(0,fun_test_Derivative);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_differentiation