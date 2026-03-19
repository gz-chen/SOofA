#ifndef LZ_DEF_LZ_H_test_math_function
#define LZ_DEF_LZ_H_test_math_function

#include<LZ_H/math/function.hpp>
using namespace Liuze;
using namespace Liuze::math;

#include<iostream>
template<typename Tp>
std::ostream & operator <<(std::ostream & out,Type_ArrayTemp<Tp> const & x){
	if(x.size()>Type_Size(0)) out<<x[0];
	for(Type_Size i=1;i<x.size();++i) out<<" "<<x[i];
	return out;
}

#if LZ_DEF_compile_part(0)
int fun_test_Comb_enum_pow_1(){
	typedef Type_Int type_int;
	typedef typename Comb_enum_pow<type_int>::type_val type_val;
	type_int i,i1,n_err=0;
	type_int n_base=3,n_pow=3;
	type_int N=powint(n_base,n_pow);
	Comb_enum_pow<type_int> iter(n_base,n_pow);
	type_val val1;
	for(i=0;i<N;++i){
		Comb_enum_pow<type_int>::index_to_value(n_base,n_pow,i,val1);
		Comb_enum_pow<type_int>::value_to_index(n_base,*iter,i1);
		if(val1!=*iter || i1!=i){
			++n_err;
			cout<<i<<":("<<(*iter)<<"); "<<i1<<"; ("<<val1<<")"<<endl;
		}
		++iter;
	}
	cout<<"Errors: "<<n_err<<endl;
	return 0;
} //fun: fun_test_Comb_enum_pow_1;

int fun_test_Comb_enum_comb_1(){
	typedef Type_Int type_int;
	typedef typename Comb_enum_comb<type_int>::type_val type_val;
	type_int i,i1,n_err=0;
	type_int n_total=5,n_select=3;
	type_int N=Comb_num::comb<type_int>(n_total,n_select);
	Comb_enum_comb<type_int> iter(n_total,n_select);
	type_val val1;
	for(i=0;i<N;++i){
		Comb_enum_comb<type_int>::index_to_value(n_total,n_select,i,val1);
		Comb_enum_comb<type_int>::value_to_index(n_total,*iter,i1);
		if(val1!=*iter || i1!=i){
			++n_err;
			cout<<i<<":("<<(*iter)<<"); "<<i1<<"; ("<<val1<<")"<<endl;
		}
		++iter;
	}
	cout<<"Errors: "<<n_err<<endl;
	return 0;
} //fun: fun_test_Comb_enum_comb_1;

int fun_test_Comb_enum_perm_1(){
	typedef Type_Int type_int;
	typedef typename Comb_enum_perm<type_int>::type_val type_val;
	type_int i,i1,n_err=0;
	type_int n_order=5;
	type_int N=Comb_num::perm<type_int>(n_order,n_order);
	Comb_enum_perm<type_int> iter(n_order);
	type_val val1;
	for(i=0;i<N;++i){
		Comb_enum_perm<type_int>::index_to_value(n_order,i,val1);
		Comb_enum_perm<type_int>::value_to_index(*iter,i1);
		if(val1!=*iter || i1!=i){
			++n_err;
			cout<<i<<":("<<(*iter)<<"); "<<i1<<"; ("<<val1<<")"<<endl;
		}
		++iter;
	}
	cout<<"Errors: "<<n_err<<endl;
	return 0;
} //fun: fun_test_Comb_enum_perm_1;

int fun_test_Comb_enum_powcombperm_1(){
	fun_test_Comb_enum_pow_1();
	cout<<endl;
	fun_test_Comb_enum_comb_1();
	cout<<endl;
	fun_test_Comb_enum_perm_1();
	cout<<endl;
	return 0;
} //fun: fun_test_Comb_enum_powcombperm_1;
LZ_DEF_compile_part_setfun(0,fun_test_Comb_enum_powcombperm_1);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_function