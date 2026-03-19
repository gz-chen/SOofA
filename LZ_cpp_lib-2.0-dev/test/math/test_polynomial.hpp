#ifndef LZ_DEF_LZ_H_test_math_polynomial
#define LZ_DEF_LZ_H_test_math_polynomial

#include<LZ_H/math/polynomial.hpp>
using namespace Liuze;
using namespace Liuze::math;

#if LZ_DEF_compile_part(0)
int fun_test_Polynomial_1(){
	typedef Type_Int type_int;
	typedef type_int type_coef;
	typedef Polynomial<type_coef> type_poly;
	auto fun_output_poly=[]
		(ostream & out,type_poly const & poly,string const & var_nom)->ostream &{
			type_int i_deg=poly.order();
			if(i_deg==type_int(0)){
				out<<"0";
				return out;
			}
			while(i_deg>type_int(0)){
				--i_deg;
				if(Is_approx_zero(poly.coef(i_deg))) continue;
				out<<" + ("<<poly.coef(i_deg)<<")*"<<var_nom<<"^"<<i_deg;
			}
			return out;
		};
	Type_ArrayTemp<type_coef> arr_coef_poly_0={0,1,0,3};
	Type_ArrayTemp<type_coef> arr_coef_poly_1={5,9,-5,4,0};
	type_poly poly_0(arr_coef_poly_0),poly_1(arr_coef_poly_1);
	Type_ArrayTemp<type_poly> list_poly={poly_0,poly_1};
	list_poly.push_back(poly_0+poly_1);
	list_poly.push_back(poly_0*poly_1);
	cout<<poly_0(type_int(2))<<endl;
	type_int i_poly;
	for(i_poly=0;i_poly<list_poly.size();++i_poly){
		fun_output_poly(cout,list_poly[i_poly],"x")<<endl;
	}
	return 0;
} //fun: fun_test_Polynomial_1;
LZ_DEF_compile_part_setfun(0,fun_test_Polynomial_1);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_polynomial