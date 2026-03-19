#ifndef LZ_DEF_LZ_H_test_stats_DOE_crit_aberration
#define LZ_DEF_LZ_H_test_stats_DOE_crit_aberration

#include<LZ_H/stats/DOE_crit_aberration.hpp>
using namespace Liuze;
using namespace Liuze::stats;

#if LZ_DEF_compile_part(0)
int fun_test_WordlengthPattern(){
	typedef Type_Size type_size;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_size> type_des;
	type_size n_level=3,n_fac=2;
	type_size i_des,i_fac;
	Liuze::math::Comb_enum_pow<type_size> iter_pow(n_level,n_fac);
	Type_ArrayTemp<type_des> des(Liuze::math::powint(n_level,n_fac));
	i_des=0;
	while(iter_pow.Is_stat_deref()){
		des[i_des].resize(n_fac,1);
		for(i_fac=0;i_fac<n_fac;++i_fac) des[i_des](i_fac,0)=(*iter_pow)[i_fac];
		++iter_pow;
		++i_des;
	}
	Iterator_constval<type_real> wt_fulldes((type_real)1/type_real(des.size()),des.size());
	Type_ArrayTemp<type_real> wt_Doptdes={0.1458,0.0802,0.1458,0.0802,0.0960,0.0802,0.1458,0.0802,0.1458};
	crit::WordlengthPattern<type_real,type_size> const WLP;
	Type_ArrayTemp<type_real> beta_fulldes=
		WLP.Alphabeta(des.begin(),des.end(),wt_fulldes.begin(),wt_fulldes.end(),
		Type_ArrayTemp<type_size>(0),Type_ArrayTemp<type_size>(n_fac,n_level));
	for(i_fac=0;i_fac<beta_fulldes.size();++i_fac) cout<<beta_fulldes[i_fac]<<" ";
	cout<<endl;
	Type_ArrayTemp<type_real> beta_Doptdes=
		WLP.Alphabeta(des.begin(),des.end(),wt_Doptdes.begin(),wt_Doptdes.end(),
		Type_ArrayTemp<type_size>(0),Type_ArrayTemp<type_size>(n_fac,n_level));
	for(i_fac=0;i_fac<beta_fulldes.size();++i_fac) cout<<beta_Doptdes[i_fac]<<" ";
	cout<<endl;
	return 0;
} //fun: fun_test_WordlengthPattern;
LZ_DEF_compile_part_setfun(0,fun_test_WordlengthPattern);
#endif

#if LZ_DEF_compile_part(1)
template<typename Tp>
std::ostream & operator <<(std::ostream & out,Type_ArrayTemp<Tp> const & x){
	if(x.size()>Type_Size(0)) out<<x[0];
	for(Type_Size i=1;i<x.size();++i) out<<" "<<x[i];
	return out;
}
int fun_test_contrast_Haar(){
	typedef Type_UInt type_int;
	typedef Type_Real type_real;
	typedef Type_MatTemp<type_real> type_mat_real;
	typedef crit::WordlengthPattern<type_real,type_int> type_WLP;
	typedef crit::StratificationPattern<type_real,type_int> type_SP;
	type_int n_lev,base,i_lev,j_lev;
	type_mat_real haarval;
	Type_ArrayTemp<type_int> haarord;
	cout.precision(2);
	cout<<fixed;
	while(true){
		cout<<"n_lev, base: ";
		cin>>n_lev>>base;
		if(n_lev<1 || base<2) break;
		haarval=type_WLP::get_contrast_Haar_value(n_lev,base);
		haarord=type_WLP::get_contrast_Haar_order(n_lev,base);
		cout<<endl<<haarord<<endl<<endl;
		cout<<haarval<<endl<<endl;
		cout<<(haarval*haarval.transpose())<<endl<<endl;
		cout<<
			(haarval*haarval.transpose()-type_real(n_lev)*type_mat_real::Identity(n_lev,n_lev))
			.squaredNorm()
			<<endl<<endl;
		cout<<"***"<<endl;
		haarval=type_SP::get_contrast_Haar_value(n_lev,base);
		haarord=type_SP::get_contrast_Haar_order(n_lev,base);
		cout<<endl<<haarord<<endl<<endl;
		cout<<haarval<<endl<<endl;
		cout<<(haarval*haarval.transpose())<<endl<<endl;
		cout<<
			(haarval*haarval.transpose()-type_real(n_lev)*type_mat_real::Identity(n_lev,n_lev))
			.squaredNorm()
			<<endl<<endl;
		cout<<"***"<<endl;
	}
	return 0;
} //fun: fun_test_contrast_Haar;
LZ_DEF_compile_part_setfun(1,fun_test_contrast_Haar);
#endif

#endif //#ifndef LZ_DEF_LZ_H_test_stats_DOE_crit_aberration