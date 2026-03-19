#ifndef LZ_DEF_LZ_H_test_math_prime
#define LZ_DEF_LZ_H_test_math_prime

#include<LZ_H/math/prime.hpp>
using namespace Liuze;
using namespace Liuze::math;

#if LZ_DEF_compile_part(0)
#include<LZ_H/math/function.hpp>
int fun_test_integer_prime_decomposition(){
	typedef Type_Int type_int;
	function<void(Type_ArrayTemp<Type_ArrayTemp<type_int> > const &)> infun_shower_decomp=
		[](Type_ArrayTemp<Type_ArrayTemp<type_int> > const & decomp)->void{
			Type_Size i;
			type_int prod=1;
			for(i=0;i<decomp.size();++i){
				cout<<"("<<decomp[i][0]<<","<<decomp[i][1]<<") ";
				prod*=powint(decomp[i][0],decomp[i][1]);
			}
			cout<<endl<<endl;
			cout<<"product: "<<prod<<endl;
		};
	type_int num;
	while(true){
		cout<<"number: ";
		cin>>num;
		if(num==type_int(0)) break;
		infun_shower_decomp(integer_prime_decomposition(num));
		cout<<endl;
	}
	return 0;
} //fun: fun_test_integer_prime_decomposition;
LZ_DEF_compile_part_setfun_auto(fun_test_integer_prime_decomposition);
#endif

#if LZ_DEF_compile_part(1)
int fun_test_PrimePowerQuerier(){
	typedef Type_UInt type_int;
	function<void(PrimePowerQuerier<type_int> const &)> infun_PrimePowerQuerier_shower=
		[](PrimePowerQuerier<type_int> const & querier)->void{
			cout<<"upper: "<<querier.get_upper()<<endl<<endl;
			Type_Size i,j;
			for(i=0;i<querier.get_table_base().size();++i){
				for(j=0;j<querier.get_table_base()[i].size();++j){
					cout<<querier.get_table_base()[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
			for(i=0;i<querier.get_table_ascending().size();++i){
				cout<<querier.get_table_ascending()[i][0]
					<<"("<<querier.get_table_ascending()[i][1]<<","
					<<querier.get_table_ascending()[i][2]<<") ";
			}
			cout<<endl<<endl;
		};
	PrimePowerQuerier<type_int> PPQ;
	infun_PrimePowerQuerier_shower(PPQ);
	type_int upper;
	while(true){
		cout<<"upper bound for prime powers: ";
		cin>>upper;
		if(upper==type_int(0)) break;
		PPQ.set_upper(upper);
		infun_PrimePowerQuerier_shower(PPQ);
	}
	return 0;
} //fun: fun_test_PrimePowerQuerier;
LZ_DEF_compile_part_setfun_auto(fun_test_PrimePowerQuerier);
#endif

#if LZ_DEF_compile_part(2)
#include<iomanip>
int fun_test_Get_prime_all(){
	typedef Type_UInt type_uint;
	typedef Type_ArrayTemp<type_uint> type_arr_uint;
	type_arr_uint arr_prime=Get_prime_all(type_uint(14219673));
	Type_Size i_bg=
		std::lower_bound(arr_prime.begin(),arr_prime.end(),type_uint(14196736))-arr_prime.begin();
	Type_Size i_ed=arr_prime.size();
	Type_Size i;
	cout<<i_ed-i_bg<<endl<<endl;
	for(i=i_bg;i<i_ed;++i){
		cout<<setbase(10)<<i<<"\t"<<arr_prime[i]<<"\t"<<setbase(16)<<uppercase<<arr_prime[i]<<endl;
	}
	return 0;
} //fun: fun_test_Get_prime_all;
LZ_DEF_compile_part_setfun_auto(fun_test_Get_prime_all);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_prime