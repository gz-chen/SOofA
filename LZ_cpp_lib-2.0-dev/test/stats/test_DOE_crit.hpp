#ifndef LZ_DEF_LZ_H_test_stats_DOE_crit
#define LZ_DEF_LZ_H_test_stats_DOE_crit

//$#include<LZ_H/math/distribution.hpp>
#include<LZ_H/stats/DOE_crit.hpp>
#include<LZ_H/stats/DOE_design_uniform.hpp>
using namespace Liuze;
using namespace Liuze::stats;

#if LZ_DEF_compile_part(0)
int fun_test_Discrepancy_Leb_infty(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	typedef Type_ArrayTemp<type_matrix> type_des;
	math::Type_RandEngine rand(100);
	uniform_real_distribution<type_real> distrU(0,1);
	type_size size=200,dim=2;
	type_des sample(size);
	type_size i_size,i_dim;
	for(i_size=0;i_size<sample.size();++i_size){
		sample[i_size].resize(dim,1);
		for(i_dim=0;i_dim<dim;++i_dim) sample[i_size](i_dim,0)=distrU(rand);
	}
	Iterator_constval<type_real> weight_des((type_real)1/type_real(sample.size()),sample.size());
	//\crit::Discrepancy_Leb_infty<type_real> disc_Linf(dim);
	//\crit::Discrepancy_canon_Leb_infty_approxLambda<type_real> disc_Linf_approx;
	crit::Discrepancy_canon_Leb_2_kernel<type_real> disc_L2_wrap
		(dim,crit::Discrepancy_canon_Leb_2_kernel<type_real>::val_stat_ker_wrap,type_real(1));
	//\cout<<disc_Linf(sample.begin(),sample.end(),weight_des.begin(),weight_des.end())<<endl;
	//\cout<<disc_Linf_approx(sample.begin(),sample.end(),weight_des.begin(),weight_des.end())<<endl;
	cout<<disc_L2_wrap(sample.begin(),sample.end(),weight_des.begin(),weight_des.end())<<endl;
	return 0;
} //fun: fun_test_Discrepancy_Leb_infty;
LZ_DEF_compile_part_setfun(0,fun_test_Discrepancy_Leb_infty);
#endif

#if LZ_DEF_compile_part(1)
int fun_test_Iterator_design_canon(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	typedef Type_ArrayTemp<type_matrix> type_des;
	math::Type_RandEngine rand(100),rand1(200);
	type_size dim=2,size=10;
	Iterator_design_canon_RandomUniform<type_real> iter_pt(dim,rand),iter_pt1(dim,rand1);
	type_size i_size;
	for(i_size=0;i_size<size;++i_size){
		cout<<(*iter_pt).transpose()<<endl;
		++iter_pt;
	}
	cout<<endl;
	for(i_size=0;i_size<size;++i_size){
		cout<<(*iter_pt1).transpose()<<endl;
		++iter_pt1;
	}
	cout<<endl;
	cout<<type_size(rand==rand1)<<" ";
	cout<<type_size(iter_pt==iter_pt1)<<" ";
	iter_pt1=iter_pt;
	cout<<type_size(iter_pt==iter_pt1)<<" ";
	iter_pt1.Set_rand_engine(rand1);
	//$cout<<type_size(iter_pt.decltype(iter_pt)::type_base::operator ==(iter_pt1))<<" ";
	cout<<type_size(iter_pt==iter_pt1)<<" ";
	cout<<endl<<endl;
	Iterator_design_canon_GoodPoint<type_real> iter_gp_cyc
		(Iterator_design_canon_GoodPoint<type_real>::Get_GoodPoint_cyclotomic(dim,7));
	Iterator_design_canon_GoodPoint<type_real> iter_gp_pow
		(Iterator_design_canon_GoodPoint<type_real>::Get_GoodPoint_power(dim,2));
	Get_design_GoodLatticePoint<type_real,type_size> CGLP(size);
	crit::Discrepancy_canon_Leb_2_kernel<type_real> disc_L2_wrap
		(dim,crit::Discrepancy_canon_Leb_2_kernel<type_real>::val_stat_ker_wrap,type_real(1));
	CGLP.Cal_best_power(dim,disc_L2_wrap);
	type_des des_glp(size);
	if(CGLP.Is_valid()) CGLP.Get_design_canon(des_glp.begin());
	Type_ArrayTemp<Iterator_wrapper<type_matrix const> > iter_wrapper{iter_gp_cyc,iter_gp_pow,des_glp.begin()};
	for(i_size=0;i_size<size;++i_size){
		cout<<(*iter_gp_cyc).transpose()<<endl;
		++iter_gp_cyc;
	}
	cout<<endl;
	for(i_size=0;i_size<size;++i_size){
		cout<<(*iter_gp_pow).transpose()<<endl;
		++iter_gp_pow;
	}
	cout<<endl;
	for(type_size i_iter=0;i_iter<iter_wrapper.size();++i_iter){
		for(i_size=0;i_size<size;++i_size){
			cout<<(*(iter_wrapper[i_iter])).transpose()<<endl;
			++(iter_wrapper[i_iter]);
		}
		cout<<endl;
	}
	iter_wrapper={iter_gp_cyc,iter_gp_pow,des_glp.begin()};
	Iterator_design_canon_WrapAroundShifter<type_real> iter_shift({0,0,des_glp.size()},iter_wrapper);
	for(i_size=0;i_size<3*size;++i_size){
		cout<<(*iter_shift).transpose()<<endl;
		++iter_shift;
	}
	cout<<endl;
	return 0;
} //fun: fun_test_Iterator_design_canon;
LZ_DEF_compile_part_setfun(1,fun_test_Iterator_design_canon);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_stats_DOE_crit