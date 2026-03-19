#ifndef LZ_DEF_LZ_H_test_math_GaloisField
#define LZ_DEF_LZ_H_test_math_GaloisField

#include<LZ_H/math/GaloisField.hpp>
using namespace Liuze;
using namespace Liuze::math;

#if LZ_DEF_compile_part(0)
int fun_test_GaloisField_0(){
	typedef Type_Size type_int;
	typedef IntegerResidueClass::Element<type_int> type_intr;
	typedef Polynomial<type_intr> type_poly;
	typedef GaloisField::Field<type_int> type_GF_base;
	typedef GaloisField::Field_tabulated<type_int> type_GF_tab;
	typedef GaloisField::Field_const_shared<type_GF_base> type_GF_base_csh;
	typedef GaloisField::Field_const_shared<type_GF_tab> type_GF_tab_csh;
	typedef type_GF_base_csh type_GF;
	typedef typename type_GF::type_space_act type_GF_act;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int i_poly,i_ele,j_ele;
	type_int GF_prime=2,GF_power=2;
	cin>>GF_prime>>GF_power;
	Type_ArrayTemp<type_poly> list_irrPoly=type_GF_act::Get_irrPoly_all(GF_prime,GF_power);
	for(i_poly=0;i_poly<list_irrPoly.size();++i_poly){
		cout<<list_irrPoly[i_poly]<<endl;
	}
	cout<<endl;
	--i_poly;
	if(list_irrPoly[i_poly].degree()!=GF_power) return 1;
	type_poly GF_irr=list_irrPoly[i_poly];
	type_GF GF(type_GF_act(GF_prime,GF_power,GF_irr));
	type_int GF_size=GF.size();
	Type_ArrayTemp<type_eleGF> ele_GF(GF_size);
	for(i_ele=0;i_ele<GF_size;++i_ele){
		ele_GF[i_ele]=type_eleGF(i_ele,GF);
		cout<<i_ele<<": "<<ele_GF[i_ele].value_polynomial()<<endl;
	}
	cout<<endl;
	for(i_ele=0;i_ele<GF_size;++i_ele){
		for(j_ele=0;j_ele<GF_size;++j_ele){
			cout<<(ele_GF[i_ele]+ele_GF[j_ele]).value_number()<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	for(i_ele=0;i_ele<GF_size;++i_ele){
		for(j_ele=0;j_ele<GF_size;++j_ele){
			cout<<(ele_GF[i_ele]*ele_GF[j_ele]).value_number()<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
	for(i_ele=0;i_ele<GF_size;++i_ele){
		cout<<ele_GF[i_ele].inv()<<" ";
	}
	cout<<endl<<endl;
	Type_MatTemp<type_eleGF> vec_ele_GF(GF_size,1);
	for(i_ele=0;i_ele<GF_size;++i_ele){
		vec_ele_GF(i_ele,0)=ele_GF[i_ele];
	}
	cout<<(-vec_ele_GF).transpose()<<endl<<endl;
	cout<<(vec_ele_GF*vec_ele_GF.transpose())<<endl;
	cout<<endl;
	return 0;
} //fun: fun_test_GaloisField_0;
LZ_DEF_compile_part_setfun(0,fun_test_GaloisField_0);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_math_GaloisField