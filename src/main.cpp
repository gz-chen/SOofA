template<int fun_id=-1>
class Fun_main{
	public:
	int operator ()(int,char*) const {return 0;}
};
#ifndef LZ_DEF_compile_part_val_makefile
	#define LZ_DEF_compile_part_val 7
#else
	#define LZ_DEF_compile_part_val LZ_DEF_compile_part_val_makefile
#endif
#define LZ_DEF_compile_part(x) LZ_DEF_compile_part_val==x
#define LZ_DEF_compile_part_setfun(id,fun) \
	template<>class Fun_main<id>\
	{public:int operator()(int argc,char *argv[])const{return fun(argc,argv);}}
#define LZ_DEF_compile_part_setfun_auto(fun) LZ_DEF_compile_part_setfun(LZ_DEF_compile_part_val,fun)
#define LZ_DEF_run_mainfun(x) return Fun_main<x>()(argc,argv)

#define DEF_debug_pos std::cout<<__FILE__<<"("<<__LINE__<<")"<<std::endl;
#define DEF_debug_out_head std::cout<<__FILE__<<"("<<__LINE__<<"): "
#define DEF_debug_out(x) DEF_debug_out_head<<(x)<<std::endl

#include<string>
#include<iostream>
#include<fstream>
using namespace std;
//#include<LZ_H/math/intresidue.hpp>
#include<LZ_H/math/prime.hpp>
#include<LZ_H/math/GaloisField.hpp>
using namespace Liuze;
using namespace Liuze::math;

#include"OofA_des_strat.hpp"

typedef Type_Size type_int;
typedef Type_Real type_real;
typedef Type_MatTemp<type_real> type_mat_real;
typedef Type_MatTemp<type_int> type_mat_int;
typedef Type_ArrayTemp<type_real> type_arr_real;
typedef Type_ArrayTemp<type_int> type_arr_int;

#if LZ_DEF_compile_part(0)
int main_COAfrac(int argc,char *argv[]){
	//$typedef GaloisField::Field_const_shared<GaloisField::Field<type_int> > type_GF;
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int n_comp,strth,n_run_frac,stat_tran_level;
	type_int i_argv=0;
	if(argc==1){
		cout<<"n_component strength: ";
		cin>>n_comp>>strth;
	} else {
		n_comp=type_int(atoi(argv[++i_argv]));
		strth=type_int(atoi(argv[++i_argv]));
	}
	type_int n_run;
	type_mat_int des;
	while(true){
		if(n_comp==type_int(0) && strth==type_int(0)){
			if(argc==1){
				cout<<"n_run n_component: ";
				cin>>n_run>>n_comp;
				if(n_run<=type_int(0) || n_comp<=type_int(0)){
					cout<<"The parameters are wrong!"<<endl;
					return 11;
				}
				des.resize(n_run,n_comp);
				cout<<"design:"<<endl;
				type_int ir,ic;
				for(ir=0;ir<n_run;++ir){
					for(ic=0;ic<n_comp;++ic) cin>>des(ir,ic);
				}
			} else {
				return 11;
			}
		} else if(n_comp<=type_int(0) || strth<=type_int(0)){
			return 11;
		} else if(strth==type_int(1)){
			des=Get_OofADesign_cycle
				(Liuze::math::Comb_num::perm<type_int,type_int>(n_comp,strth),n_comp);
		} else if(strth==type_int(2) && n_comp==type_int(10)){
			Type_ArrayTemp<type_poly> list_irrPoly=
				type_GF::type_space_act::Get_irrPoly_all(type_int(3),type_int(2));
			type_int i_poly=list_irrPoly.size();
			--i_poly;
			if(argc==1){
				cout<<"reduction_of_run_by_fraction: ";
				cin>>n_run_frac;
				cout<<"tran_level: "
					"0: origin; "
					"1: ordered_post0; 2: collapsed_post0; 3: ordered_pre0; 4: collapsed_pre0:";
				cin>>stat_tran_level;
			} else {
				n_run_frac=type_int(atoi(argv[++i_argv]));
				stat_tran_level=type_int(atoi(argv[++i_argv]));
			}
			type_GF GF(typename type_GF::type_space_act
				(type_int(3),type_int(2),list_irrPoly[i_poly]));
			des=Get_COA2frac(GF,n_run_frac).block(0,0,type_int(9)/n_run_frac,type_int(9));
			applyOnMat_entrywise
				(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(stat_tran_level)),des);
			type_mat_int mat_aug(20,10);
			mat_aug<<
				4,3,7,8,5,6,9,1,2,0,
				5,8,1,2,6,9,0,3,7,4,
				6,2,3,7,9,0,4,8,1,5,
				9,7,8,1,0,4,5,2,3,6,
				0,1,2,3,4,5,6,7,8,9,
				7,5,9,6,1,4,2,3,8,0,
				1,6,3,8,4,2,0,5,9,7,
				4,8,5,9,2,0,7,6,3,1,
				2,9,6,3,0,7,1,8,5,4,
				0,3,8,5,7,1,4,9,6,2,
				8,6,5,4,7,3,9,2,1,0,
				7,4,2,1,3,9,0,6,5,8,
				3,1,6,5,9,0,8,4,2,7,
				9,5,4,2,0,8,7,1,6,3,
				0,2,1,6,8,7,3,5,4,9,
				5,7,9,3,2,8,1,6,4,0,
				2,3,6,4,8,1,0,7,9,5,
				8,4,7,9,1,0,5,3,6,2,
				1,9,3,6,0,5,2,4,7,8,
				0,6,4,7,5,2,8,9,3,1;
			des=Get_COA_addComponent_custom(des,mat_aug);
		} else if(strth==type_int(2) || strth==type_int(3)){
			if(n_comp==type_int(1)){
				des.resize(1,1);
				des(0,0)=type_int(0);
			} else {
				Type_ArrayTemp<type_arr_int> decomp=
					integer_prime_decomposition(strth==type_int(2) ? n_comp : n_comp-type_int(1));
				if(decomp.size()!=Type_Size(1)){
					cout<<"The parameters are wrong!"<<endl;
					return 12;
				}
				Type_ArrayTemp<type_poly> list_irrPoly=
					type_GF::type_space_act::Get_irrPoly_all(decomp[0][0],decomp[0][1]);
				type_int i_poly=list_irrPoly.size();
				if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=decomp[0][1]){
					cout<<"No irr. poly."<<endl;
					return 13;
				}
				--i_poly;
				if(argc==1){
					cout<<"reduction_of_run_by_fraction: ";
					cin>>n_run_frac;
					cout<<"tran_level: "
						"0: origin; "
						"1: ordered_post0; 2: collapsed_post0; 3: ordered_pre0; 4: collapsed_pre0:";
					cin>>stat_tran_level;
				} else {
					n_run_frac=type_int(atoi(argv[++i_argv]));
					stat_tran_level=type_int(atoi(argv[++i_argv]));
				}
				type_GF GF(typename type_GF::type_space_act
					(decomp[0][0],decomp[0][1],list_irrPoly[i_poly]));
				des=Get_COA2frac(GF,n_run_frac);
				applyOnMat_entrywise
					(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(stat_tran_level)),des);
				if(strth==type_int(3)){
					des=Get_COA3frac(GF,n_run_frac,&des);
				}
			}
		} else if(strth==type_int(4) && n_comp==type_int(7)){
			Type_ArrayTemp<type_poly> list_irrPoly=
				type_GF::type_space_act::Get_irrPoly_all(type_int(5),type_int(1));
			type_int i_poly=list_irrPoly.size();
			--i_poly;
			if(argc==1){
				cout<<"reduction_of_run_by_fraction: ";
				cin>>n_run_frac;
				cout<<"tran_level: "
					"0: origin; "
					"1: ordered_post0; 2: collapsed_post0; 3: ordered_pre0; 4: collapsed_pre0:";
				cin>>stat_tran_level;
			} else {
				n_run_frac=type_int(atoi(argv[++i_argv]));
				stat_tran_level=type_int(atoi(argv[++i_argv]));
			}
			type_GF GF(typename type_GF::type_space_act
				(type_int(5),type_int(1),list_irrPoly[i_poly]));
			des=Get_COA2frac(GF,n_run_frac);
			applyOnMat_entrywise
				(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(stat_tran_level)),des);
			des=Get_COA3frac(GF,n_run_frac,&des);
			type_mat_int mat_aug(14,7);
			mat_aug<<
				6,0,1,2,3,4,5,
				0,6,1,2,3,5,4,
				0,1,6,2,3,4,5,
				0,1,2,6,3,5,4,
				0,1,2,5,6,3,4,
				0,1,2,4,3,6,5,
				0,1,2,3,4,5,6,
				6,0,1,2,5,3,4,
				0,6,1,2,5,4,3,
				0,1,6,2,4,5,3,
				0,1,2,6,4,3,5,
				0,1,2,5,6,4,3,
				0,1,2,4,5,6,3,
				0,1,2,3,5,4,6;
			des=Get_COA_addComponent_custom(des,mat_aug);
		} else {
			return 11;
		}
		break;
	}
	n_run=des.rows();
	if(n_run==type_int(0)) return 15;
	type_int flag_file;
	if(argc==1){
		cout<<"Write to file? (0/1): ";
		cin>>flag_file;
	} else {
		flag_file=type_int(atoi(argv[++i_argv]));
	}
	string file_dir="./des/";
	string filename=
		string("OofAdes_")+to_string(n_run)+"_"+to_string(n_comp)+"_"+to_string(strth);
	string file_comment=string("#")
		+"("+to_string(n_run)+","+to_string(n_comp)+","+to_string(strth);
	if(strth==type_int(2) || strth==type_int(3)){
		filename=filename+"_"+to_string(n_run_frac)+"_"+to_string(stat_tran_level);
		file_comment=file_comment+","+to_string(n_run_frac)+","+to_string(stat_tran_level);
	}
	file_comment=file_comment+")";
	fstream fout;
	multi_ostream<0> mout;
	if(flag_file!=type_int(0)){
		fout.open(file_dir+filename+".txt",ios::out);
		mout={&cout,&fout};
	} else {
		mout={&cout};
	}
	mout<<file_comment<<endl<<endl;
	mout<<des<<endl;
	if(argc==1){
		type_int i_comp;
		char cstr_sep[26];
		for(i_comp=0;i_comp<25;++i_comp) cstr_sep[i_comp]='*';
		cstr_sep[25]=0;
		Type_ArrayTemp<type_int> arr_coll(n_comp);
		type_mat_int des_coll;
		type_int strth_coll;
		while(true){
			cout<<endl<<cstr_sep<<endl;
			cout<<"to apply another collapsing? (1/0): ";
			cin>>stat_tran_level;
			if(stat_tran_level==type_int(0)) break;
			cout<<"input the collapsing: ";
			for(i_comp=0;i_comp<n_comp;++i_comp) cin>>arr_coll[i_comp];
			if(flag_file!=type_int(0)){
				fout<<endl<<cstr_sep<<endl<<"another collapsing:";
				for(i_comp=0;i_comp<n_comp;++i_comp) fout<<" "<<arr_coll[i_comp];
				fout<<endl;
			}
			des_coll=evalOnMat_entrywise(LevelCollapsing<type_int>(arr_coll),des);
			mout<<endl<<des_coll<<endl;
			if(stat_tran_level==type_int(2)){
				cout<<endl<<"strength_of_stratum_orthogonality: ";
				cin>>strth_coll;
				mout<<endl<<"Strength of stratum orthogonality: "<<
					(Check_PofC_SO(des,strth_coll,LevelCollapsing<type_int>(arr_coll)) ?
					strth_coll : type_int(0))
					<<endl;
			}
		}
	}
	if(flag_file!=type_int(0)) fout.close();
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_COAfrac);
#endif

#if LZ_DEF_compile_part(1)

#include<random>
#include<LZ_H/stats/DOE_crit_aberration.hpp>
using namespace Liuze::stats;

template<typename FieldType,typename IntType=Liuze::Type_UInt,
	typename RandEngineType=std::default_random_engine,
	typename=typename std::enable_if<
		std::is_integral<typename FieldType::value_type>::value &&
		std::is_integral<IntType>::value &&
		std::is_same<typename Liuze::math::GaloisField::tag_type_space<FieldType>::type,
		Liuze::math::GaloisField::tag_type_space_tabulated>::value
	>::type>
Liuze::Type_MatTemp<IntType> Get_COA2frac_rand
(FieldType const & field,IntType n_run_frac=0,Liuze::Type_Stat tran_level=0,
RandEngineType * rand_eng=NULL){
	typedef typename FieldType::value_type type_val;
	typedef typename FieldType::type_element type_eleGF;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_comp=field.size();
	if(n_comp<=type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return type_res::Zero(1,1);
	if(n_run_frac<=type_int(0)) n_run_frac=1;
	if(n_comp%n_run_frac==type_int(0)){
		uniform_int_distribution<type_int> Distr_unif(type_int(0),n_run_frac-type_int(1));
		RandEngineType * rand_eng_act= rand_eng ? rand_eng : new RandEngineType;
		type_int n_add=n_comp/n_run_frac;
		type_int n_run=n_add*(n_comp-type_int(1));
		type_res res(n_run,n_comp);
		type_int ic,ir0,ir1,im,ir0r;
		for(ir1=0,im=1;ir1<n_run;ir1+=n_add,++im){
			for(ir0=0;ir0<n_add;++ir0){
				ir0r=ir0*n_run_frac+Distr_unif(*rand_eng_act);
				for(ic=0;ic<n_comp;++ic){
					res(ir1+ir0,ic)=type_int(field(im)*field(ic)+field(ir0r));
				}
			}
		}
		if(tran_level==Type_Stat(2)){
			res/=n_run_frac;
		}
		if(!rand_eng) delete rand_eng_act;
		return res;
	} else if((n_comp-type_int(1))%n_run_frac==type_int(0)){
		type_int n_mul=(n_comp-type_int(1))/n_run_frac;
		type_int n_run=n_comp*n_mul;
		type_res res(n_run,n_comp);
		type_int ic,ir0,ir1;
		type_eleGF eleGFprim(type_int(field.space_act().primitive()),field);
		type_eleGF t_sum;
		for(ir0=0;ir0<n_comp;++ir0){
			for(ic=0;ic<n_comp;++ic){
				t_sum=field(ir0)+field(ic);
				res(ir0,ic)=type_int(t_sum);
				for(ir1=1;ir1<n_mul;++ir1){
					res(ir1*n_comp+ir0,ic)=type_int(eleGFprim.pow(ir1)*t_sum);
				}
			}
		}
		if(tran_level!=Type_Stat(0)){
			Liuze::Type_ArrayTemp<type_int> arr_tran(n_comp);
			if(tran_level==Type_Stat(1)){
				for(ic=1;ic<n_comp;++ic){
					ir0=field.space_act().table_log(ic-type_int(1));
					arr_tran[ic]=(ir0%n_mul)*n_run_frac+(ir0/n_mul);
				}
				arr_tran[0]=n_comp-type_int(1);
			} else if(tran_level==Type_Stat(2)){
				for(ic=1;ic<n_comp;++ic){
					arr_tran[ic]=Liuze::math::mod
						(type_int(field.space_act().table_log(ic-type_int(1))),n_mul);
				}
				arr_tran[0]=n_mul;
			}
			for(ir0=0;ir0<n_run;++ir0){
				for(ic=0;ic<n_comp;++ic){
					res(ir0,ic)=arr_tran[res(ir0,ic)];
				}
			}
		}
		return res;
	} else {
		return type_res(0,0);
	}
} //fun: Get_COA2frac_rand;

int main_test_StratificationPattern(int argc,char *argv[]){
	typedef crit::WordlengthPattern<type_real,type_int> type_crit_WLP;
	typedef crit::StratificationPattern<type_real,type_int> type_crit_SP;
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	std::default_random_engine rand_eng;
	type_crit_WLP WLP;
	type_crit_SP SP;
	Iterator_MatRowTranspose<type_mat_int const> iter_des_bg,iter_des_ed;
	Iterator_constval<type_real> iter_wt_bg,iter_wt_ed;
	type_arr_int arr_n_level,arr_base_strat;
	type_int i_ans;
	type_int n_comp,strth,n_run_frac,stat_tran_level,ord_max;
	type_int i_argv=0;
	if(argc==1){
		cout<<"n_component strength: ";
		cin>>n_comp>>strth;
	} else {
		n_comp=type_int(atoi(argv[++i_argv]));
		strth=type_int(atoi(argv[++i_argv]));
	}
	type_mat_int des;
	while(true){
		if(strth==type_int(2)){
			type_int n_comp_prime,n_comp_power;
			if(argc==1){
				cout<<"prime power: ";
				cin>>n_comp_prime>>n_comp_power;
			} else {
				n_comp_prime=type_int(atoi(argv[++i_argv]));
				n_comp_power=type_int(atoi(argv[++i_argv]));
			}
			if(Liuze::math::powint(n_comp_prime,n_comp_power)!=n_comp){
				cout<<"The parameters are wrong!"<<endl;
				return 12;
			}
			Type_ArrayTemp<type_poly> list_irrPoly=
				type_GF::type_space_act::Get_irrPoly_all(n_comp_prime,n_comp_power);
			type_int i_poly=list_irrPoly.size();
			if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=n_comp_power){
				cout<<"No irr. poly."<<endl;
				return 13;
			}
			--i_poly;
			if(argc==1){
				cout<<"reduction_of_run_by_fraction: ";
				cin>>n_run_frac;
				cout<<"tran_level: "
					"0: origin; "
					"1: ordered_post0; 2: collapsed_post0; 3: ordered_pre0; 4: collapsed_pre0:";
				cin>>stat_tran_level;
				cout<<"maximum_order: ";
				cin>>ord_max;
			} else {
				n_run_frac=type_int(atoi(argv[++i_argv]));
				stat_tran_level=type_int(atoi(argv[++i_argv]));
				ord_max=type_int(atoi(argv[++i_argv]));
			}
			iter_wt_bg.assign(type_real(1));
			iter_wt_ed=iter_wt_bg.end();
			arr_n_level.assign(n_comp,n_comp);
			arr_base_strat.assign(n_comp,n_comp/n_run_frac);
			cout<<"base_for_stratification: "<<n_comp/n_run_frac<<endl<<endl;
			i_ans=1;
			while(i_ans){
				des=Get_COA2frac_rand(
					type_GF(typename type_GF::type_space_act
					(n_comp_prime,n_comp_power,list_irrPoly[i_poly])),
					n_run_frac,Type_Stat(stat_tran_level),&rand_eng);
				//`cout<<des<<endl<<endl;
				iter_des_bg=des;
				iter_des_ed=iter_des_bg.end();
				/*
				cout<<"GWLP: "
					<<WLP.Alphabeta
					(iter_des_bg,iter_des_ed,iter_wt_bg,iter_wt_ed,arr_n_level,type_arr_int(0))
					<<endl<<endl;
				cout<<"BWLP: "
					<<WLP.Alphabeta
					(iter_des_bg,iter_des_ed,iter_wt_bg,iter_wt_ed,type_arr_int(0),arr_n_level)
					<<endl<<endl;
				*/
				cout<<"SP: "
					<<SP.Stratification
					(iter_des_bg,iter_des_ed,iter_wt_bg,iter_wt_ed,arr_n_level,arr_base_strat)
					<<endl<<endl;
				cout<<"SP: "
					<<SP.Stratification
					(iter_des_bg,iter_des_ed,iter_wt_bg,iter_wt_ed,arr_n_level,arr_base_strat,ord_max)
					<<endl<<endl;
				cout<<"continue ? ";
				cin>>i_ans;
				cout<<Type_ArrayTemp<char>(25,'*')<<endl<<endl;
			}
		} else {
			return 11;
		}
		break;
	}
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_test_StratificationPattern);
#endif

#if LZ_DEF_compile_part(2)
int main_Find_COA2_normSO3strat3(int argc,char *argv[]){
	//$typedef GaloisField::Field_const_shared<GaloisField::Field<type_int> > type_GF;
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int n_comp;
	type_int i_argv=0;
	if(argc==1){
		cout<<"n_component: ";
		cin>>n_comp;
	} else {
		n_comp=type_int(atoi(argv[++i_argv]));
	}
	type_mat_int des;
	type_int n_comp_prime,n_comp_power;
	if(argc==1){
		cout<<"prime power: ";
		cin>>n_comp_prime>>n_comp_power;
	} else {
		n_comp_prime=type_int(atoi(argv[++i_argv]));
		n_comp_power=type_int(atoi(argv[++i_argv]));
	}
	if(Liuze::math::powint(n_comp_prime,n_comp_power)!=n_comp){
		cout<<"The parameters are wrong!"<<endl;
		return 12;
	}
	Type_ArrayTemp<type_poly> list_irrPoly=
		type_GF::type_space_act::Get_irrPoly_all(n_comp_prime,n_comp_power);
	type_int i_poly=list_irrPoly.size();
	if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=n_comp_power){
		cout<<"No irr. poly."<<endl;
		return 13;
	}
	--i_poly;
	type_GF GF(typename type_GF::type_space_act
		(n_comp_prime,n_comp_power,list_irrPoly[i_poly]));
	type_int res=Find_COA2_normSO3strat3(GF);
	cout<<endl<<res<<endl<<endl;
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_Find_COA2_normSO3strat3);
#endif

#if LZ_DEF_compile_part(3)
int main_test_GYD(int argc,char *argv[]){
	#if false
	type_int n_comp=13;
	type_mat_int BBD(4,n_comp);
	BBD<<
		0,0,0,0,1,1,1,2,2,3,3,4,5,
		1,2,4,6,2,5,7,3,6,4,7,8,9,
		3,8,5,10,4,6,11,5,7,6,8,9,10,
		9,12,7,11,10,8,12,11,9,12,10,11,12;
	BBD.transposeInPlace();
	type_mat_int des(n_comp,n_comp);
	des.block(0,0,n_comp,BBD.cols())=Get_GYD_rowregular(BBD,n_comp);
	des.block(0,BBD.cols(),n_comp,n_comp-BBD.cols())=
		Get_GYD_rowregular(Get_rowBIBD_complement(BBD,n_comp),n_comp);
	cout<<des<<endl<<endl;
	cout<<Get_PofC(des)<<endl<<endl;
	#endif
	#if false
	type_mat_int BIBD(31,10);
	BIBD<<
		28,1,6,9,12,17,18,22,23,25,
		28,2,0,10,13,18,19,23,24,26,
		28,3,1,11,7,19,20,24,25,27,
		28,4,2,12,8,20,14,25,26,21,
		28,5,3,13,9,14,15,26,27,22,
		28,6,4,7,10,15,16,27,21,23,
		28,0,5,8,11,16,17,21,22,24,
		29,8,13,16,19,3,4,22,23,25,
		29,9,7,17,20,4,5,23,24,26,
		29,10,8,18,14,5,6,24,25,27,
		29,11,9,19,15,6,0,25,26,21,
		29,12,10,20,16,0,1,26,27,22,
		29,13,11,14,17,1,2,27,21,23,
		29,7,12,15,18,2,3,21,22,24,
		30,15,20,2,5,10,11,22,23,25,
		30,16,14,3,6,11,12,23,24,26,
		30,17,15,4,0,12,13,24,25,27,
		30,18,16,5,1,13,7,25,26,21,
		30,19,17,6,2,7,8,26,27,22,
		30,20,18,0,3,8,9,27,21,23,
		30,14,19,1,4,9,10,21,22,24,
		3,6,5,10,13,12,17,20,19,21,
		4,0,6,11,7,13,18,14,20,22,
		5,1,0,12,8,7,19,15,14,23,
		6,2,1,13,9,8,20,16,15,24,
		0,3,2,7,10,9,14,17,16,25,
		1,4,3,8,11,10,15,18,17,26,
		2,5,4,9,12,11,16,19,18,27,
		28,29,30,0,1,2,3,4,5,6,
		29,30,28,7,8,9,10,11,12,13,
		30,29,28,14,15,16,17,18,19,20;
	BIBD=Get_rowBIBD_derived(BIBD);
	cout<<BIBD<<endl<<endl;
	type_int n_comp=10;
	type_mat_int des(BIBD.rows(),n_comp);
	des.block(0,0,BIBD.rows(),BIBD.cols())=Get_GYD_rowregular(BIBD,n_comp);
	des.block(0,BIBD.cols(),BIBD.rows(),n_comp-BIBD.cols())=
		Get_GYD_rowregular(Get_rowBIBD_complement(BIBD,n_comp),n_comp);
	cout<<des.transpose()<<endl<<endl;
	cout<<Get_PofC(des).transpose()<<endl<<endl;
	#endif
	#if true
	{
		type_int n_comp=10;
		type_mat_int BIBD0[3],BIBD(30,4);
		type_int i;
		for(i=0;i<3;++i) BIBD0[i].resize(4,15);
		BIBD0[0].block(0,0,4,15)<<
			0,0,0,0,0,0,1,1,1,1,2,2,2,3,3,
			1,1,2,3,5,6,2,3,4,5,3,4,5,4,4,
			2,4,4,7,7,8,7,6,7,6,5,8,6,5,6,
			3,5,6,8,9,9,8,9,9,8,9,9,7,8,7;
		BIBD0[1].block(0,0,4,15)<<
			0,0,0,0,0,0,1,1,1,1,2,2,2,3,3,
			1,1,2,3,5,6,2,3,4,5,3,4,5,4,4,
			2,4,4,7,7,8,7,6,7,6,5,8,6,5,6,
			3,5,6,8,9,9,9,9,8,8,8,9,7,9,7;
		BIBD0[2].block(0,0,4,15)<<
			0,0,0,0,0,0,1,1,1,1,2,2,2,3,3,
			1,1,2,3,4,5,2,3,4,5,3,4,5,4,5,
			2,4,6,8,6,7,8,6,7,6,4,7,6,6,7,
			3,5,7,9,8,9,9,7,9,8,5,8,9,9,8;
		for(i=0;i<3;++i) BIBD0[i].transposeInPlace();
		type_mat_int des(BIBD.rows(),n_comp);
		for(i=0;i<3;++i){
			BIBD.block(15,0,15,4)=BIBD.block(0,0,15,4)=BIBD0[i];
			des.block(0,0,BIBD.rows(),BIBD.cols())=Get_GYD_rowregular(BIBD,n_comp);
			des.block(0,BIBD.cols(),BIBD.rows(),n_comp-BIBD.cols())=
				Get_GYD_rowregular(Get_rowBIBD_complement(BIBD,n_comp),n_comp);
			cout<<n_comp<<" "<<BIBD.rows()<<" "<<BIBD.cols()<<" "<<i<<": "<<endl<<endl;
			cout<<des.transpose()<<endl<<endl;
			cout<<Get_PofC(des).transpose()<<endl<<endl;
		}
	}
	{
		type_int n_comp=6;
		type_mat_int BIBD(3,30);
		BIBD.block(0,0,3,10)<<
			0,0,0,0,0,1,1,1,2,2,
			1,1,2,3,4,2,3,4,3,3,
			2,3,4,5,5,5,4,5,4,5;
		BIBD.transposeInPlace();
		BIBD.block(20,0,10,3)=BIBD.block(10,0,10,3)=BIBD.block(0,0,10,3);
		type_mat_int des(BIBD.rows(),n_comp);
		des.block(0,0,BIBD.rows(),BIBD.cols())=Get_GYD_rowregular(BIBD,n_comp);
		des.block(0,BIBD.cols(),BIBD.rows(),n_comp-BIBD.cols())=
			Get_GYD_rowregular(Get_rowBIBD_complement(BIBD,n_comp),n_comp);
		cout<<n_comp<<" "<<BIBD.rows()<<" "<<BIBD.cols()<<": "<<endl<<endl;
		cout<<des.transpose()<<endl<<endl;
		cout<<Get_PofC(des).transpose()<<endl<<endl;
	}
	#endif
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_test_GYD);
#endif

#if LZ_DEF_compile_part(4)
#include"DOE22S20_simu.hpp"
int main_DOE22S20_simu(int argc,char *argv[]){
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int n_comp,strth,n_run_frac;
	type_int i_argv=0;
	if(argc==1){
		cout<<"n_component strength: ";
		cin>>n_comp>>strth;
	} else {
		n_comp=type_int(atoi(argv[++i_argv]));
		strth=type_int(atoi(argv[++i_argv]));
	}
	type_int n_run;
	type_mat_int mat_PofC_train,mat_PofC_test;
	while(true){
		if(strth==n_comp+type_int(1)){
			if(argc==1){
				cout<<"n_run n_component: ";
				cin>>n_run>>n_comp;
				if(n_run<=type_int(0) || n_comp<=type_int(0)){
					cout<<"The parameters are wrong!"<<endl;
					return 11;
				}
				mat_PofC_train.resize(n_run,n_comp);
				cout<<"design:"<<endl;
				type_int ir,ic;
				for(ir=0;ir<n_run;++ir){
					for(ic=0;ic<n_comp;++ic) cin>>mat_PofC_train(ir,ic);
				}
			} else {
				return 11;
			}
		} else if(n_comp<=type_int(0) || strth<=type_int(0)){
			return 11;
		} else if(strth==type_int(1)){
			mat_PofC_train=Get_OofADesign_cycle
				(Liuze::math::Comb_num::perm<type_int,type_int>(n_comp,strth),n_comp);
		} else if(strth==type_int(2) || strth==type_int(3)){
			if(n_comp==type_int(1)){
				mat_PofC_train.resize(1,1);
				mat_PofC_train(0,0)=type_int(0);
			} else {
				Type_ArrayTemp<type_arr_int> decomp=
					integer_prime_decomposition(strth==type_int(2) ? n_comp : n_comp-type_int(1));
				if(decomp.size()!=Type_Size(1)){
					cout<<"The parameters are wrong!"<<endl;
					return 12;
				}
				Type_ArrayTemp<type_poly> list_irrPoly=
					type_GF::type_space_act::Get_irrPoly_all(decomp[0][0],decomp[0][1]);
				type_int i_poly=list_irrPoly.size();
				if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=decomp[0][1]){
					cout<<"No irr. poly."<<endl;
					return 13;
				}
				--i_poly;
				if(argc==1){
					cout<<"reduction_of_run_by_fraction: ";
					cin>>n_run_frac;
				} else {
					n_run_frac=type_int(atoi(argv[++i_argv]));
				}
				type_GF GF(typename type_GF::type_space_act
					(decomp[0][0],decomp[0][1],list_irrPoly[i_poly]));
				mat_PofC_train=Get_COA2frac(GF,n_run_frac);
				applyOnMat_entrywise
					(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(1)),mat_PofC_train);
				if(strth==type_int(3)){
					mat_PofC_train=Get_COA3frac(GF,n_run_frac,&mat_PofC_train);
				}
			}
		} else {
			return 11;
		}
		break;
	}
	if(argc==1){
		cout<<"strength: ";
		cin>>strth;
	} else {
		strth=type_int(atoi(argv[++i_argv]));
	}
	while(true){
		if(strth==n_comp+type_int(1)){
			if(argc==1){
				cout<<"n_run: ";
				cin>>n_run;
				if(n_run<=type_int(0) || n_comp<=type_int(0)){
					cout<<"The parameters are wrong!"<<endl;
					return 11;
				}
				mat_PofC_test.resize(n_run,n_comp);
				cout<<"design:"<<endl;
				type_int ir,ic;
				for(ir=0;ir<n_run;++ir){
					for(ic=0;ic<n_comp;++ic) cin>>mat_PofC_test(ir,ic);
				}
			} else {
				return 11;
			}
		} else if(n_comp<=type_int(0) || strth<=type_int(0)){
			return 11;
		} else if(strth==type_int(1)){
			mat_PofC_test=Get_OofADesign_cycle
				(Liuze::math::Comb_num::perm<type_int,type_int>(n_comp,strth),n_comp);
		} else if(strth==type_int(2) || strth==type_int(3)){
			if(n_comp==type_int(1)){
				mat_PofC_test.resize(1,1);
				mat_PofC_test(0,0)=type_int(0);
			} else {
				Type_ArrayTemp<type_arr_int> decomp=
					integer_prime_decomposition(strth==type_int(2) ? n_comp : n_comp-type_int(1));
				if(decomp.size()!=Type_Size(1)){
					cout<<"The parameters are wrong!"<<endl;
					return 12;
				}
				Type_ArrayTemp<type_poly> list_irrPoly=
					type_GF::type_space_act::Get_irrPoly_all(decomp[0][0],decomp[0][1]);
				type_int i_poly=list_irrPoly.size();
				if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=decomp[0][1]){
					cout<<"No irr. poly."<<endl;
					return 13;
				}
				--i_poly;
				if(argc==1){
					cout<<"reduction_of_run_by_fraction: ";
					cin>>n_run_frac;
				} else {
					n_run_frac=type_int(atoi(argv[++i_argv]));
				}
				type_GF GF(typename type_GF::type_space_act
					(decomp[0][0],decomp[0][1],list_irrPoly[i_poly]));
				mat_PofC_test=Get_COA2frac(GF,n_run_frac);
				applyOnMat_entrywise
					(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(1)),mat_PofC_test);
				if(strth==type_int(3)){
					mat_PofC_test=Get_COA3frac(GF,n_run_frac,&mat_PofC_test);
				}
			}
		} else {
			return 11;
		}
		break;
	}
	type_int i_comp;
	type_int seed_rand;
	type_real pow_cost;
	type_arr_int vec_strat_id(n_comp);
	for(i_comp=0;i_comp<n_comp;++i_comp) vec_strat_id[i_comp]=i_comp;
	type_int i_vec_strat,n_vec_strat;
	if(argc==1){
		cout<<"number_of_stratification_mappings: ";
		cin>>n_vec_strat;
	} else {
		n_vec_strat=type_int(atoi(argv[++i_argv]));
	}
	Type_ArrayTemp<type_arr_int> list_vec_strat(n_vec_strat,type_arr_int(n_comp));
	if(argc==1){
		cout<<"seed_for_random_engine: ";
		cin>>seed_rand;
		cout<<"power_for_cost: ";
		cin>>pow_cost;
		cout<<"stratification_mapping: ";
		for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
			for(i_comp=0;i_comp<n_comp;++i_comp) cin>>list_vec_strat[i_vec_strat][i_comp];
		}
	} else {
		seed_rand=type_int(atoi(argv[++i_argv]));
		pow_cost=type_real(atof(argv[++i_argv]));
		for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
			for(i_comp=0;i_comp<n_comp;++i_comp){
				list_vec_strat[i_vec_strat][i_comp]=type_int(atoi(argv[++i_argv]));
			}
		}
	}
	type_int i_model,n_model_0=5,n_model=n_model_0+n_vec_strat;
	Type_ArrayTemp<string> list_model_name={"PWO","poly1","poly2","poly2var2","CP"};
	list_model_name.resize(n_model);
	Type_ArrayTemp<function<type_mat_real(type_mat_int const &)> >
		list_model_Get_mat_model(n_model);
	Type_ArrayTemp<function<type_int(type_int const &)> > list_model_Get_rank_model(n_model);
	list_model_Get_mat_model[0]=Get_mat_model_OofA_PWO<type_real,type_mat_int>;
	list_model_Get_mat_model[1]=Get_mat_model_OofA_poly1<type_real,type_mat_int>;
	list_model_Get_mat_model[2]=Get_mat_model_OofA_poly2<type_real,type_mat_int>;
	list_model_Get_mat_model[3]=Get_mat_model_OofA_poly2var2<type_real,type_mat_int>;
	list_model_Get_mat_model[4]=bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
		std::placeholders::_1,vec_strat_id);
	list_model_Get_rank_model[0]=Get_rank_model_OofA_PWO<type_int>;
	list_model_Get_rank_model[1]=Get_rank_model_OofA_poly1<type_int>;
	list_model_Get_rank_model[2]=Get_rank_model_OofA_poly2<type_int>;
	list_model_Get_rank_model[3]=Get_rank_model_OofA_poly2var2<type_int>;
	list_model_Get_rank_model[4]=[&vec_strat_id](type_int const &)->type_int{
			return Get_rank_model_OofA_SCP1<type_int>(vec_strat_id);
		};
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		list_model_name[n_model_0+i_vec_strat]="SCP";
		list_model_Get_mat_model[n_model_0+i_vec_strat]=
			bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
			std::placeholders::_1,list_vec_strat[i_vec_strat]);
		list_model_Get_rank_model[n_model_0+i_vec_strat]=
			[&list_vec_strat,i_vec_strat](type_int const &)->type_int{
				return Get_rank_model_OofA_SCP1<type_int>(list_vec_strat[i_vec_strat]);
			};
	}
	DOE22S20_simu_LinReg<type_real> data_linreg;
	DOE22S20_data_generator<type_real>::rand.seed(seed_rand);
	DOE22S20_data_generator<type_real> data_gen(n_comp,pow_cost);
	mat_PofC_train=act_perm_OofADesign_left(
		Get_OofARandDesign(type_int(1),n_comp,DOE22S20_data_generator<type_real>::rand),
		mat_PofC_train);
	mat_PofC_test=act_perm_OofADesign_left(
		Get_OofARandDesign(type_int(1),n_comp,DOE22S20_data_generator<type_real>::rand),
		mat_PofC_test);
	type_mat_real mat_resp_train=data_gen.response(Get_PofC(mat_PofC_train));
	type_mat_real mat_resp_test=data_gen.response(Get_PofC(mat_PofC_test));
	type_real SSerr;
	cout<<"time:"<<endl;
	cout<<data_gen.time()<<endl<<endl;
	cout<<"weight:"<<endl;
	cout<<data_gen.weight()<<endl<<endl;
	cout<<"stratification:"<<endl;
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		cout<<list_vec_strat[i_vec_strat]<<endl;
	}
	cout<<endl;
	cout<<"model rank_model rank_mat_model est_var SS_err_train SS_err_test"<<endl;
	for(i_model=0;i_model<n_model;++i_model){
		cout<<list_model_name[i_model]<<" "<<list_model_Get_rank_model[i_model](n_comp)<<" ";
		data_linreg.compute(list_model_Get_mat_model[i_model](mat_PofC_train),mat_resp_train);
		cout<<data_linreg.rank()<<" ";
		SSerr=(mat_resp_train-data_linreg.mat_model()*data_linreg.mat_coef()).squaredNorm();
		cout<<SSerr/(type_real(mat_PofC_train.rows())-type_real(data_linreg.rank()))<<" "
			<<SSerr/type_real(mat_PofC_train.rows())<<" ";
		SSerr=(mat_resp_test-
			list_model_Get_mat_model[i_model](mat_PofC_test)*data_linreg.mat_coef()).squaredNorm();
		cout<<SSerr/type_real(mat_PofC_test.rows())<<endl;
	}
	cout<<endl;
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_DOE22S20_simu);
#endif

#if LZ_DEF_compile_part(5)
#include"DOE22S20_simu.hpp"
int main_DesCritVal_simu(int argc,char *argv[]){
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int n_comp,strth,n_run_frac;
	type_int i_argv=0;
	if(argc==1){
		cout<<"n_component strength: ";
		cin>>n_comp>>strth;
	} else {
		n_comp=type_int(atoi(argv[++i_argv]));
		strth=type_int(atoi(argv[++i_argv]));
	}
	type_int n_run;
	type_mat_int mat_PofC;
	while(true){
		if(n_comp==type_int(0) && strth==type_int(0)){
			if(argc==1){
				cout<<"n_run n_component: ";
				cin>>n_run>>n_comp;
				if(n_run<=type_int(0) || n_comp<=type_int(0)){
					cout<<"The parameters are wrong!"<<endl;
					return 11;
				}
				mat_PofC.resize(n_run,n_comp);
				cout<<"design:"<<endl;
				type_int ir,ic;
				for(ir=0;ir<n_run;++ir){
					for(ic=0;ic<n_comp;++ic) cin>>mat_PofC(ir,ic);
				}
			} else {
				return 11;
			}
		} else if(n_comp<=type_int(0) || strth<=type_int(0)){
			return 11;
		} else if(strth==type_int(1)){
			mat_PofC=Get_OofADesign_cycle
				(Liuze::math::Comb_num::perm<type_int,type_int>(n_comp,strth),n_comp);
		} else if(strth==type_int(2) || strth==type_int(3)){
			if(n_comp==type_int(1)){
				mat_PofC.resize(1,1);
				mat_PofC(0,0)=type_int(0);
			} else {
				Type_ArrayTemp<type_arr_int> decomp=
					integer_prime_decomposition(strth==type_int(2) ? n_comp : n_comp-type_int(1));
				if(decomp.size()!=Type_Size(1)){
					cout<<"The parameters are wrong!"<<endl;
					return 12;
				}
				Type_ArrayTemp<type_poly> list_irrPoly=
					type_GF::type_space_act::Get_irrPoly_all(decomp[0][0],decomp[0][1]);
				type_int i_poly=list_irrPoly.size();
				if(i_poly==0 || list_irrPoly[i_poly-1].degree()!=decomp[0][1]){
					cout<<"No irr. poly."<<endl;
					return 13;
				}
				--i_poly;
				if(argc==1){
					cout<<"reduction_of_run_by_fraction: ";
					cin>>n_run_frac;
				} else {
					n_run_frac=type_int(atoi(argv[++i_argv]));
				}
				type_GF GF(typename type_GF::type_space_act
					(decomp[0][0],decomp[0][1],list_irrPoly[i_poly]));
				mat_PofC=Get_COA2frac(GF,n_run_frac);
				applyOnMat_entrywise
					(Get_COAfrac_collapser(GF,n_run_frac,Type_Stat(1)),mat_PofC);
				if(strth==type_int(3)){
					mat_PofC=Get_COA3frac(GF,n_run_frac,&mat_PofC);
				}
			}
		} else {
			return 11;
		}
		break;
	}
	type_int i_comp;
	type_int seed_rand;
	type_arr_int vec_strat_id(n_comp);
	for(i_comp=0;i_comp<n_comp;++i_comp) vec_strat_id[i_comp]=i_comp;
	type_int i_vec_strat,n_vec_strat;
	if(argc==1){
		cout<<"number_of_stratification_mappings: ";
		cin>>n_vec_strat;
	} else {
		n_vec_strat=type_int(atoi(argv[++i_argv]));
	}
	Type_ArrayTemp<type_arr_int> list_vec_strat(n_vec_strat,type_arr_int(n_comp));
	if(argc==1){
		cout<<"seed_for_random_engine: ";
		cin>>seed_rand;
		cout<<"stratification_mapping: ";
		for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
			for(i_comp=0;i_comp<n_comp;++i_comp) cin>>list_vec_strat[i_vec_strat][i_comp];
		}
	} else {
		seed_rand=type_int(atoi(argv[++i_argv]));
		for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
			for(i_comp=0;i_comp<n_comp;++i_comp){
				list_vec_strat[i_vec_strat][i_comp]=type_int(atoi(argv[++i_argv]));
			}
		}
	}
	type_int i_model,n_model_0=5,n_model=n_model_0+n_vec_strat;
	Type_ArrayTemp<string> list_model_name={"PWO","poly1","poly2","poly2var2","CP"};
	list_model_name.resize(n_model);
	Type_ArrayTemp<function<type_mat_real(type_mat_int const &)> >
		list_model_Get_mat_model(n_model);
	Type_ArrayTemp<function<type_int(type_int const &)> > list_model_Get_rank_model(n_model);
	list_model_Get_mat_model[0]=Get_mat_model_OofA_PWO<type_real,type_mat_int>;
	list_model_Get_mat_model[1]=Get_mat_model_OofA_poly1<type_real,type_mat_int>;
	list_model_Get_mat_model[2]=Get_mat_model_OofA_poly2<type_real,type_mat_int>;
	list_model_Get_mat_model[3]=Get_mat_model_OofA_poly2var2<type_real,type_mat_int>;
	list_model_Get_mat_model[4]=bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
		std::placeholders::_1,vec_strat_id);
	list_model_Get_rank_model[0]=Get_rank_model_OofA_PWO<type_int>;
	list_model_Get_rank_model[1]=Get_rank_model_OofA_poly1<type_int>;
	list_model_Get_rank_model[2]=Get_rank_model_OofA_poly2<type_int>;
	list_model_Get_rank_model[3]=Get_rank_model_OofA_poly2var2<type_int>;
	list_model_Get_rank_model[4]=[&vec_strat_id](type_int const &)->type_int{
			return Get_rank_model_OofA_SCP1<type_int>(vec_strat_id);
		};
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		list_model_name[n_model_0+i_vec_strat]="SCP";
		list_model_Get_mat_model[n_model_0+i_vec_strat]=
			bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
			std::placeholders::_1,list_vec_strat[i_vec_strat]);
		list_model_Get_rank_model[n_model_0+i_vec_strat]=
			[&list_vec_strat,i_vec_strat](type_int const &)->type_int{
				return Get_rank_model_OofA_SCP1<type_int>(list_vec_strat[i_vec_strat]);
			};
	}
	typedef std::function<type_real(type_mat_int const &)> type_crit;
	typedef Liuze::stats::crit::Discrepancy_canon_Leb_2_kernel<type_real> type_crit_disc;
	function<type_real(type_mat_int const &,type_crit_disc &)>
		loss_disc=[](type_mat_int const & mat_PofC,type_crit_disc & crit)->type_real{
			type_int ir,ic;
			type_int n_run=mat_PofC.rows(),n_comp=mat_PofC.cols();
			Type_ArrayTemp<type_mat_real> list_des_pt(n_run);
			for(ir=0;ir<n_run;++ir){
				list_des_pt[ir].resize(n_comp,1);
				for(ic=0;ic<n_comp;++ic){
					list_des_pt[ir](ic,0)=
						(type_real(mat_PofC(ir,ic))+type_real(0.5))/type_real(n_comp);
				}
			}
			Iterator_constval<type_real> iter_wt(type_real(1)/type_real(n_run));
			return crit(list_des_pt.begin(),list_des_pt.end(),iter_wt,iter_wt.end());
		};
	type_crit_disc crit_disc_centered(n_comp,type_crit_disc::val_stat_ker_centered);
	type_crit_disc crit_disc_wrap(n_comp,type_crit_disc::val_stat_ker_wrap);
	type_crit_disc crit_disc_mixture(n_comp,type_crit_disc::val_stat_ker_mixture);
	type_crit_disc::type_kernel_discrete crit_disc_kernel_discrete(n_comp);
	crit_disc_kernel_discrete.Set_level(n_comp);
	crit_disc_kernel_discrete.Set_par_ab();
	crit_disc_kernel_discrete.Set_init();
	type_crit_disc crit_disc_discrete(crit_disc_kernel_discrete);
	type_crit loss_disc_centered=
		[&loss_disc,&crit_disc_centered](type_mat_int const & mat_PofC)->type_real{
			return loss_disc(mat_PofC,crit_disc_centered);
		};
	type_crit loss_disc_wrap=
		[&loss_disc,&crit_disc_wrap](type_mat_int const & mat_PofC)->type_real{
			return loss_disc(mat_PofC,crit_disc_wrap);
		};
	type_crit loss_disc_mixture=
		[&loss_disc,&crit_disc_mixture](type_mat_int const & mat_PofC)->type_real{
			return loss_disc(mat_PofC,crit_disc_mixture);
		};
	type_crit loss_disc_discrete=
		[&loss_disc,&crit_disc_discrete](type_mat_int const & mat_PofC)->type_real{
			return loss_disc(mat_PofC,crit_disc_discrete);
		};
	type_crit negloss_mindist_Leb_2=[](type_mat_int const & mat_PofC)->type_real{
			Calculator_norm_L_stream<type_real,type_real,type_int> Cal
				(Calculator_norm_L_stream<type_real,type_real,type_int>::val_power_infinity,false);
			type_int ir,jr,ic;
			type_int nr=mat_PofC.rows(),nc=mat_PofC.cols();
			type_real dist,dist0;
			for(ir=0;ir<nr;++ir){
				for(jr=0;jr<ir;++jr){
					dist=type_real(0);
					for(ic=0;ic<nc;++ic){
						dist0=type_real(mat_PofC(ir,ic))-type_real(mat_PofC(jr,ic));
						dist+=dist0*dist0;
					}
					Cal<<(-dist);
				}
			}
			return -Cal.Get_norm();
		};
	type_crit negloss_D_PWO=[](type_mat_int const & mat_PofC)->type_real{
			return Get_Eff_PWO_D<type_real,type_int>(Get_PofC(mat_PofC));
		};
	DOE22S20_simu_LinReg<type_real> data_linreg;
	DOE22S20_data_generator<type_real>::rand.seed(seed_rand);
	DesCritVal_data_generator<type_int,type_real> data_gen
		(mat_PofC,negloss_mindist_Leb_2,DesCritVal_data_generator<type_int,type_real>::stat_leftact);
	type_mat_int mat_PofC_train=
		Get_OofARandDesign(type_int(1),n_comp,DOE22S20_data_generator<type_real>::rand);
	mat_PofC_train=act_perm_OofADesign_left(mat_PofC_train,mat_PofC);
	type_mat_int mat_PofC_test=
		Get_OofARandDesign(type_int(1),n_comp,DOE22S20_data_generator<type_real>::rand);
	mat_PofC_test=act_perm_OofADesign_left(mat_PofC_test,mat_PofC);
	type_mat_real mat_resp_train=data_gen.response(mat_PofC_train);
	type_mat_real mat_resp_test=data_gen.response(mat_PofC_test);
	type_real SSerr;
	cout<<"stratification:"<<endl;
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		cout<<list_vec_strat[i_vec_strat]<<endl;
	}
	cout<<endl;
	cout<<"model rank_model rank_mat_model est_var SS_err_train SS_err_test"<<endl;
	for(i_model=0;i_model<n_model;++i_model){
		cout<<list_model_name[i_model]<<" "<<list_model_Get_rank_model[i_model](n_comp)<<" ";
		data_linreg.compute(list_model_Get_mat_model[i_model](mat_PofC_train),mat_resp_train);
		cout<<data_linreg.rank()<<" ";
		SSerr=(mat_resp_train-data_linreg.mat_model()*data_linreg.mat_coef()).squaredNorm();
		cout<<SSerr/(type_real(mat_PofC_train.rows())-type_real(data_linreg.rank()))<<" "
			<<SSerr<<" ";
		SSerr=(mat_resp_test-
			list_model_Get_mat_model[i_model](mat_PofC_test)*data_linreg.mat_coef()).squaredNorm();
		cout<<SSerr<<endl;
	}
	cout<<endl;
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_DesCritVal_simu);
#endif

#if LZ_DEF_compile_part(6)
int main_SO_numrun_min(int argc,char *argv[]){
	type_int n_comp,strth,n_strat;
	type_arr_int vec_stratsize;
	type_int i;
	while(true){
		cout<<"n_component strength n_stratum size_strata: ";
		cin>>n_comp>>strth>>n_strat;
		vec_stratsize.resize(n_strat);
		for(i=0;i<n_strat;++i) cin>>vec_stratsize[i];
		cout<<"minimum n_run: "<<Cal_SO_numrun_min(n_comp,strth,vec_stratsize)<<endl;
		#if false
		cout<<"continue? (0 / 1) ";
		cin>>i;
		if(i==0) break;
		#endif
		cout<<endl;
	}
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_SO_numrun_min);
#endif

#if LZ_DEF_compile_part(7)
#include"DOE22S20_simu.hpp"
int main_DesCritVal_datafree(int argc,char *argv[]){
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	type_int i_des,i_comp,i_vec_strat,i_model,i_crit;
	type_int n_comp=8,n_run=28,n_des=3,n_vec_strat=1;
	type_int n_model_0=5,n_model=n_model_0+n_vec_strat;
	Type_ArrayTemp<type_mat_int> list_mat_PofC(n_des);
	list_mat_PofC[0]=Get_OofAFullDesign<type_int>(n_comp);
	list_mat_PofC[1].resize(8,28);
	list_mat_PofC[1]<<
		0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6,
		1,3,5,7,2,0,6,4,3,1,7,5,4,6,0,2,5,7,1,3,6,4,2,0,7,5,3,1,
		2,0,6,4,4,6,0,2,6,4,2,0,5,7,1,3,7,5,3,1,1,3,5,7,3,1,7,5,
		3,1,7,5,6,4,2,0,5,7,1,3,1,3,5,7,2,0,6,4,7,5,3,1,4,6,0,2,
		4,6,0,2,5,7,1,3,1,3,5,7,7,5,3,1,3,1,7,5,2,0,6,4,6,4,2,0,
		5,7,1,3,7,5,3,1,2,0,6,4,3,1,7,5,6,4,2,0,4,6,0,2,1,3,5,7,
		6,4,2,0,1,3,5,7,7,5,3,1,2,0,6,4,4,6,0,2,3,1,7,5,5,7,1,3,
		7,5,3,1,3,1,7,5,4,6,0,2,6,4,2,0,1,3,5,7,5,7,1,3,2,0,6,4;
	list_mat_PofC[1].transposeInPlace();
	list_mat_PofC[2].resize(8,28);
	list_mat_PofC[2]<<
		0,1,2,3,4,5,6,7,0,6,1,7,2,4,3,5,0,4,5,1,7,3,2,6,0,3,6,5,
		1,0,3,2,5,4,7,6,6,0,7,1,4,2,5,3,4,0,1,5,3,7,6,2,3,0,5,6,
		2,3,0,1,6,7,4,5,1,7,0,6,3,5,2,4,5,1,0,4,2,6,7,3,6,5,0,3,
		3,2,1,0,7,6,5,4,7,1,6,0,5,3,4,2,1,5,4,0,6,2,3,7,5,6,3,0,
		4,5,6,7,0,1,2,3,2,4,3,5,0,6,1,7,7,3,2,6,0,4,5,1,1,2,7,4,
		5,4,7,6,1,0,3,2,4,2,5,3,6,0,7,1,3,7,6,2,4,0,1,5,2,1,4,7,
		6,7,4,5,2,3,0,1,3,5,2,4,1,7,0,6,2,6,7,3,5,1,0,4,7,4,1,2,
		7,6,5,4,3,2,1,0,5,3,4,2,7,1,6,0,6,2,3,7,1,5,4,0,4,7,2,1;
	list_mat_PofC[2].transposeInPlace();
	type_arr_int vec_strat_id(n_comp);
	for(i_comp=0;i_comp<n_comp;++i_comp) vec_strat_id[i_comp]=i_comp;
	Type_ArrayTemp<type_arr_int> list_vec_strat(n_vec_strat,type_arr_int(n_comp));
	for(i_comp=0;i_comp<n_comp;++i_comp) list_vec_strat[0][i_comp]=i_comp/type_int(2);
	Type_ArrayTemp<string> list_model_name={"PWO","FP","QP","SP","CP"};
	list_model_name.resize(n_model);
	Type_ArrayTemp<function<type_mat_real(type_mat_int const &)> >
		list_model_Get_mat_model(n_model);
	Type_ArrayTemp<function<type_int(type_int const &)> > list_model_Get_rank_model(n_model);
	list_model_Get_mat_model[0]=Get_mat_model_OofA_PWO<type_real,type_mat_int>;
	list_model_Get_mat_model[1]=Get_mat_model_OofA_poly1<type_real,type_mat_int>;
	list_model_Get_mat_model[2]=Get_mat_model_OofA_poly2<type_real,type_mat_int>;
	list_model_Get_mat_model[3]=Get_mat_model_OofA_poly2var2<type_real,type_mat_int>;
	list_model_Get_mat_model[4]=bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
		std::placeholders::_1,vec_strat_id);
	list_model_Get_rank_model[0]=Get_rank_model_OofA_PWO<type_int>;
	list_model_Get_rank_model[1]=Get_rank_model_OofA_poly1<type_int>;
	list_model_Get_rank_model[2]=Get_rank_model_OofA_poly2<type_int>;
	list_model_Get_rank_model[3]=Get_rank_model_OofA_poly2var2<type_int>;
	list_model_Get_rank_model[4]=[&vec_strat_id](type_int const &)->type_int{
			return Get_rank_model_OofA_SCP1<type_int>(vec_strat_id);
		};
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		list_model_name[n_model_0+i_vec_strat]="SCP";
		list_model_Get_mat_model[n_model_0+i_vec_strat]=
			bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
			std::placeholders::_1,list_vec_strat[i_vec_strat]);
		list_model_Get_rank_model[n_model_0+i_vec_strat]=
			[&list_vec_strat,i_vec_strat](type_int const &)->type_int{
				return Get_rank_model_OofA_SCP1<type_int>(list_vec_strat[i_vec_strat]);
			};
	}
	typedef std::function<type_real(type_mat_int const &,type_int const &)> type_crit_model;
	type_crit_model crit_model_matmeanT=
		[&list_model_Get_mat_model,&list_model_Get_rank_model]
		(type_mat_int const & mat_PofC,type_int const & i_model)->type_real{
			type_mat_real mat_infor=list_model_Get_mat_model[i_model](mat_PofC);
			mat_infor=mat_infor.transpose()*mat_infor/type_real(mat_PofC.rows());
			return mat_infor.trace();
		};
	type_crit_model crit_model_matinforSqTraceInv=
		[&list_model_Get_mat_model,&list_model_Get_rank_model]
		(type_mat_int const & mat_PofC,type_int const & i_model)->type_real{
			type_mat_real mat_infor=list_model_Get_mat_model[i_model](mat_PofC);
			mat_infor=mat_infor.transpose()*mat_infor/type_real(mat_PofC.rows());
			mat_infor=mat_infor*mat_infor;
			return type_real(1)/mat_infor.trace();
		};
	type_crit_model crit_model_matmeanD=
		[&list_model_Get_mat_model,&list_model_Get_rank_model]
		(type_mat_int const & mat_PofC,type_int const & i_model)->type_real{
			type_mat_real mat_infor=list_model_Get_mat_model[i_model](mat_PofC);
			mat_infor=mat_infor.transpose()*mat_infor/type_real(mat_PofC.rows());
			Matrix_Inverse_MoorePenrose<type_real> MatInv(mat_infor);
			return type_int(MatInv.Get_dim())
				==list_model_Get_rank_model[i_model](type_int(mat_PofC.cols())) ?
				type_real(std::pow(MatInv.Get_det_nonzero(),
				type_real(1)/type_real(MatInv.Get_dim())))
				: type_real(0);
		};
	type_int n_crit=3;
	Type_ArrayTemp<type_crit_model> list_crit_model(n_crit);
	list_crit_model[0]=crit_model_matmeanT;
	list_crit_model[1]=crit_model_matinforSqTraceInv;
	list_crit_model[2]=crit_model_matmeanD;
	{
		//`Type_ArrayTemp<type_mat_real> list_val_crit(n_des,type_mat_real(n_crit,n_model));
		type_mat_real val_crit(n_des,n_crit*n_model);
		for(i_des=0;i_des<n_des;++i_des){
			for(i_model=0;i_model<n_model;++i_model){
				for(i_crit=0;i_crit<n_crit;++i_crit){
					val_crit(i_des,i_model*n_crit+i_crit)=
						list_crit_model[i_crit](list_mat_PofC[i_des],i_model);
				}
			}
		}
		type_mat_real eff_crit(val_crit.rows(),val_crit.cols());
		for(i_des=0;i_des<n_des;++i_des){
			eff_crit.row(i_des)=(val_crit.row(i_des).array()/val_crit.row(0).array()).matrix();
		}
		cout<<val_crit<<endl<<endl;
		cout<<eff_crit<<endl<<endl;
		cout<<val_crit.transpose()<<endl<<endl;
		cout<<eff_crit.transpose()<<endl<<endl;
	}
	{
		type_int n_rand=1000,i_rand,rand_seed=100;
		std::default_random_engine rand_eng(rand_seed);
		type_real val,val_1;
		Calculator_norm_L_stream<type_real,type_real>
			Cal_Exp(type_real(1),false),Cal_Var(type_real(2),false);
		Type_ArrayTemp<Type_ArrayTemp<type_arr_int> > list_vec_strat_des(n_des);
		type_int i_vec_strat;
		list_vec_strat_des[0].assign(1,type_arr_int(n_comp));
		list_vec_strat_des[1].assign(8,type_arr_int(n_comp));
		list_vec_strat_des[2].assign(1,type_arr_int(n_comp));
		for(i_comp=0;i_comp<n_comp;++i_comp){
			list_vec_strat_des[0][0][i_comp]=i_comp;
			//`list_vec_strat_des[1][0][i_comp]=i_comp/type_int(2);
			list_vec_strat_des[2][0][i_comp]=type_int(0);
		}
		list_vec_strat_des[1][0]={0,0,1,1,2,2,3,3};
		list_vec_strat_des[1][1]={0,0,1,1,2,2,2,2};
		list_vec_strat_des[1][2]={0,0,1,1,1,1,2,2};
		list_vec_strat_des[1][3]={0,0,0,0,1,1,2,2};
		list_vec_strat_des[1][4]={0,0,0,0,1,1,1,1};
		list_vec_strat_des[1][5]={0,0,1,1,1,1,1,1};
		list_vec_strat_des[1][6]={0,0,0,0,0,0,1,1};
		list_vec_strat_des[1][7]={0,0,0,0,0,0,0,0};
		type_int n_val_crit=0,i_val_crit;
		for(i_des=0;i_des<n_des;++i_des) n_val_crit+=type_int(list_vec_strat_des[i_des].size());
		type_mat_real
			val_crit_Exp(n_val_crit,n_crit*n_model),val_crit_SD(n_val_crit,n_crit*n_model);
		for(i_model=0;i_model<n_model;++i_model){
			for(i_crit=0;i_crit<n_crit;++i_crit){
				val_crit_Exp(0,i_model*n_crit+i_crit)=
					list_crit_model[i_crit](list_mat_PofC[0],i_model);
			}
		}
		for(i_model=0;i_model<n_model;++i_model){
			for(i_crit=0;i_crit<n_crit;++i_crit){
				val_1=val_crit_Exp(0,i_model*n_crit+i_crit);
				for(i_val_crit=1,i_des=1;i_des<n_des;++i_des){
					for(i_vec_strat=0;i_vec_strat<type_int(list_vec_strat_des[i_des].size());
						++i_vec_strat){
						Cal_Exp.Set_init();
						Cal_Var.Set_init();
						for(i_rand=0;i_rand<n_rand;++i_rand){
							val=list_crit_model[i_crit]
								(Get_mat_PofC_isomorph_symbperm_rand
								(list_mat_PofC[i_des],LevelCollapsing<type_int>
								(list_vec_strat_des[i_des][i_vec_strat]),
								rand_eng),i_model);
							val/=val_1;
							Cal_Exp<<val;
							Cal_Var<<val;
						}
						val=Cal_Exp.Get_mean();
						cout<<i_val_crit<<" "<<(i_model*n_crit+i_crit)<<" "<<val<<endl;
						val_crit_Exp(i_val_crit,i_model*n_crit+i_crit)=val;
						val_crit_SD(i_val_crit,i_model*n_crit+i_crit)=type_real(std::sqrt(
							(Cal_Var.Get_mean()-val*val)*
							(type_real(1)+type_real(1)/(type_real(n_rand)-type_real(1)))));
						++i_val_crit;
					}
				}
			}
		}
		cout<<val_crit_Exp<<endl<<endl;
		cout<<val_crit_SD<<endl<<endl;
		cout<<val_crit_Exp.transpose()<<endl<<endl;
		cout<<val_crit_SD.transpose()<<endl<<endl;
	}
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_DesCritVal_datafree);
#endif

#if LZ_DEF_compile_part(8)
#include"DOE22S20_simu.hpp"
int main_SCP_4DrugExperiment(int argc,char *argv[]){
	typedef GaloisField::Field_const_shared<GaloisField::Field_tabulated<type_int> > type_GF;
	typedef typename type_GF::type_space_act::type_poly type_poly;
	typedef GaloisField::Element<type_int,type_GF> type_eleGF;
	#if false
	type_int i_des,i_comp,i_vec_strat,i_model;
	type_int n_comp=4,n_des=3,n_vec_strat=3;
	type_int n_model_0=5,n_model=n_model_0+n_vec_strat;
	Type_ArrayTemp<type_mat_int> list_mat_PofC(n_des);
	list_mat_PofC[0].resize(4,24);
	list_mat_PofC[0]<<
		0,0,2,2,2,2,2,0,2,3,1,3,1,0,3,3,0,1,1,3,0,3,1,1,
		2,2,0,3,3,1,0,3,1,2,2,0,0,3,1,2,1,2,3,0,1,1,3,0,
		3,1,1,1,0,3,3,2,0,0,3,2,2,1,2,1,2,0,2,1,3,0,0,3,
		1,3,3,0,1,0,1,1,3,1,0,1,3,2,0,0,3,3,0,2,2,2,2,2;
	list_mat_PofC[0].transposeInPlace();
	list_mat_PofC[0]=Get_PofC(list_mat_PofC[0]);
	Type_ArrayTemp<type_mat_real> list_mat_resp(n_des);
	list_mat_resp[0].resize(list_mat_PofC[0].rows(),1);
	list_mat_resp[0]<<
		56.5,55.4,53.5,53.4,52.9,51.4,51.2,51.2,50.8,46.8,46.4,46.4,
		46.1,43.3,42.1,41.8,41.1,39.5,39.4,39.1,37.5,37.2,34.4,27.8;
	{
		type_int i_run;
		type_arr_int arr_i_run={0,2,4,5,11,13,15,16,17,18,21,23};
		list_mat_PofC[1].resize(arr_i_run.size(),n_comp);
		list_mat_resp[1].resize(arr_i_run.size(),1);
		for(i_run=0;i_run<arr_i_run.size();++i_run){
			list_mat_PofC[1].row(i_run)=list_mat_PofC[0].row(arr_i_run[i_run]);
			list_mat_resp[1](i_run,0)=list_mat_resp[0](arr_i_run[i_run],0);
		}
		arr_i_run={1,3,6,7,8,9,10,12,14,19,20,22};
		list_mat_PofC[2].resize(arr_i_run.size(),n_comp);
		list_mat_resp[2].resize(arr_i_run.size(),1);
		for(i_run=0;i_run<arr_i_run.size();++i_run){
			list_mat_PofC[2].row(i_run)=list_mat_PofC[0].row(arr_i_run[i_run]);
			list_mat_resp[2](i_run,0)=list_mat_resp[0](arr_i_run[i_run],0);
		}
	}
	type_arr_int vec_strat_id(n_comp);
	for(i_comp=0;i_comp<n_comp;++i_comp) vec_strat_id[i_comp]=i_comp;
	Type_ArrayTemp<type_arr_int> list_vec_strat(n_vec_strat,type_arr_int(n_comp));
	list_vec_strat[0]={0,1,1,2};
	list_vec_strat[1]={0,0,1,2};
	list_vec_strat[2]={0,1,2,2};
	Type_ArrayTemp<string> list_model_name={"PWO","FP","QP","SP","CP"};
	list_model_name.resize(n_model);
	Type_ArrayTemp<function<type_mat_real(type_mat_int const &)> >
		list_model_Get_mat_model(n_model);
	Type_ArrayTemp<function<type_int(type_int const &)> > list_model_Get_rank_model(n_model);
	list_model_Get_mat_model[0]=Get_mat_model_OofA_PWO<type_real,type_mat_int>;
	list_model_Get_mat_model[1]=Get_mat_model_OofA_poly1<type_real,type_mat_int>;
	list_model_Get_mat_model[2]=Get_mat_model_OofA_poly2<type_real,type_mat_int>;
	list_model_Get_mat_model[3]=Get_mat_model_OofA_poly2var2<type_real,type_mat_int>;
	list_model_Get_mat_model[4]=bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
		std::placeholders::_1,vec_strat_id);
	list_model_Get_rank_model[0]=Get_rank_model_OofA_PWO<type_int>;
	list_model_Get_rank_model[1]=Get_rank_model_OofA_poly1<type_int>;
	list_model_Get_rank_model[2]=Get_rank_model_OofA_poly2<type_int>;
	list_model_Get_rank_model[3]=Get_rank_model_OofA_poly2var2<type_int>;
	list_model_Get_rank_model[4]=[&vec_strat_id](type_int const &)->type_int{
			return Get_rank_model_OofA_SCP1<type_int>(vec_strat_id);
		};
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		list_model_name[n_model_0+i_vec_strat]="SCP";
		list_model_Get_mat_model[n_model_0+i_vec_strat]=
			bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
			std::placeholders::_1,list_vec_strat[i_vec_strat]);
		list_model_Get_rank_model[n_model_0+i_vec_strat]=
			[&list_vec_strat,i_vec_strat](type_int const &)->type_int{
				return Get_rank_model_OofA_SCP1<type_int>(list_vec_strat[i_vec_strat]);
			};
	}
	DOE22S20_simu_LinReg<type_real> data_linreg;
	type_mat_real mat_resp_pred;
	type_real SSerr;
	cout<<"model rank_model rank_mat_model est_var SS_err_train SS_err"<<endl;
	for(i_model=0;i_model<n_model;++i_model){
		for(i_des=0;i_des<n_des;++i_des){
			cout<<list_model_name[i_model]<<" "<<list_model_Get_rank_model[i_model](n_comp)<<" ";
			data_linreg.compute(list_model_Get_mat_model[i_model](list_mat_PofC[i_des]),
				list_mat_resp[i_des]);
			cout<<data_linreg.rank()<<" ";
			SSerr=(list_mat_resp[i_des]-data_linreg.mat_model()*data_linreg.mat_coef()
				).squaredNorm();
			cout<<SSerr/(type_real(list_mat_PofC[i_des].rows())-type_real(data_linreg.rank()))<<" "
				<<SSerr/type_real(list_mat_PofC[i_des].rows())<<" ";
			mat_resp_pred=
				list_model_Get_mat_model[i_model](list_mat_PofC[0])*data_linreg.mat_coef();
			SSerr=(list_mat_resp[0]-mat_resp_pred).squaredNorm();
			cout<<SSerr/type_real(list_mat_PofC[0].rows())<<endl;
			cout<<mat_resp_pred.transpose()<<endl<<endl;
		}
	}
	cout<<endl;
	#endif
	#if true
	type_int i_des,i_comp,i_vec_strat,i_model;
	type_int n_comp=5,n_des=3,n_vec_strat=4;
	type_int n_model_0=5,n_model=n_model_0+n_vec_strat;
	Type_ArrayTemp<type_mat_int> list_mat_PofC(n_des);
	list_mat_PofC[0].resize(5,40);
	list_mat_PofC[0]<<
		0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
		1,2,3,4,0,2,3,4,0,1,3,4,0,1,2,4,0,1,2,3,1,2,3,4,0,2,3,4,0,1,3,4,0,1,2,4,0,1,2,3,
		4,3,2,1,2,0,4,3,4,3,1,0,1,0,4,2,3,2,1,0,3,1,4,2,3,4,2,0,1,4,0,3,4,2,0,1,2,0,3,1,
		3,1,4,2,3,4,2,0,1,4,0,3,4,2,0,1,2,0,3,1,2,4,1,3,4,3,0,2,3,0,4,1,2,4,1,0,1,3,0,2,
		2,4,1,3,4,3,0,2,3,0,4,1,2,4,1,0,1,3,0,2,4,3,2,1,2,0,4,3,4,3,1,0,1,0,4,2,3,2,1,0;
	list_mat_PofC[0].transposeInPlace();
	list_mat_PofC[0]=Get_PofC(list_mat_PofC[0]);
	Type_ArrayTemp<type_mat_real> list_mat_resp(n_des);
	list_mat_resp[0].resize(list_mat_PofC[0].rows(),1);
	list_mat_resp[0]<<
		20.23,23.55,23.61,21.59,13.63,23.85,21.47,28.38,31.27,26.08,
		29.43,30.52,15.57,4.93,18.47,25.23,26.75,31.96,25.62,19.50,
		10.96,22.06,16.25,16.37,7.72,24.65,12.09,20.40,23.37,22.35,
		25.99,26.30,13.84,5.53,17.97,23.40,26.49,24.31,19.71,20.35;
	list_mat_PofC[1]=list_mat_PofC[0].block(0,0,20,n_comp);
	list_mat_PofC[2]=list_mat_PofC[0].block(20,0,20,n_comp);
	list_mat_resp[1]=list_mat_resp[0].block(0,0,20,1);
	list_mat_resp[2]=list_mat_resp[0].block(20,0,20,1);
	type_arr_int vec_strat_id(n_comp);
	for(i_comp=0;i_comp<n_comp;++i_comp) vec_strat_id[i_comp]=i_comp;
	Type_ArrayTemp<type_arr_int> list_vec_strat(n_vec_strat,type_arr_int(n_comp));
	list_vec_strat[0]={0,0,1,2,3};
	list_vec_strat[1]={0,1,1,2,3};
	list_vec_strat[2]={0,1,2,2,3};
	list_vec_strat[3]={0,1,2,3,3};
	Type_ArrayTemp<string> list_model_name={"PWO","FP","QP","SP","CP"};
	list_model_name.resize(n_model);
	Type_ArrayTemp<function<type_mat_real(type_mat_int const &)> >
		list_model_Get_mat_model(n_model);
	Type_ArrayTemp<function<type_int(type_int const &)> > list_model_Get_rank_model(n_model);
	list_model_Get_mat_model[0]=Get_mat_model_OofA_PWO<type_real,type_mat_int>;
	list_model_Get_mat_model[1]=Get_mat_model_OofA_poly1<type_real,type_mat_int>;
	list_model_Get_mat_model[2]=Get_mat_model_OofA_poly2<type_real,type_mat_int>;
	list_model_Get_mat_model[3]=Get_mat_model_OofA_poly2var2<type_real,type_mat_int>;
	list_model_Get_mat_model[4]=bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
		std::placeholders::_1,vec_strat_id);
	list_model_Get_rank_model[0]=Get_rank_model_OofA_PWO<type_int>;
	list_model_Get_rank_model[1]=Get_rank_model_OofA_poly1<type_int>;
	list_model_Get_rank_model[2]=Get_rank_model_OofA_poly2<type_int>;
	list_model_Get_rank_model[3]=Get_rank_model_OofA_poly2var2<type_int>;
	list_model_Get_rank_model[4]=[&vec_strat_id](type_int const &)->type_int{
			return Get_rank_model_OofA_SCP1<type_int>(vec_strat_id);
		};
	for(i_vec_strat=0;i_vec_strat<n_vec_strat;++i_vec_strat){
		list_model_name[n_model_0+i_vec_strat]="SCP";
		list_model_Get_mat_model[n_model_0+i_vec_strat]=
			bind(Get_mat_model_OofA_SCP1<type_real,type_mat_int,type_arr_int>,
			std::placeholders::_1,list_vec_strat[i_vec_strat]);
		list_model_Get_rank_model[n_model_0+i_vec_strat]=
			[&list_vec_strat,i_vec_strat](type_int const &)->type_int{
				return Get_rank_model_OofA_SCP1<type_int>(list_vec_strat[i_vec_strat]);
			};
	}
	DOE22S20_simu_LinReg<type_real> data_linreg;
	type_mat_real mat_resp_pred;
	type_real SSerr;
	cout<<"i_model i_des model rank_model rank_mat_model est_var SS_err_train SS_err"<<endl;
	for(i_model=0;i_model<n_model;++i_model){
		for(i_des=0;i_des<n_des;++i_des){
			cout<<i_model<<" "<<i_des<<" ";
			cout<<list_model_name[i_model]<<" "<<list_model_Get_rank_model[i_model](n_comp)<<" ";
			data_linreg.compute(list_model_Get_mat_model[i_model](list_mat_PofC[i_des]),
				list_mat_resp[i_des]);
			cout<<data_linreg.rank()<<" ";
			SSerr=(list_mat_resp[i_des]-data_linreg.mat_model()*data_linreg.mat_coef()
				).squaredNorm();
			cout<<SSerr/(type_real(list_mat_PofC[i_des].rows())-type_real(data_linreg.rank()))<<" "
				<<SSerr/type_real(list_mat_PofC[i_des].rows())<<" ";
			mat_resp_pred=
				list_model_Get_mat_model[i_model](list_mat_PofC[0])*data_linreg.mat_coef();
			SSerr=(list_mat_resp[0]-mat_resp_pred).squaredNorm();
			cout<<SSerr/type_real(list_mat_PofC[0].rows())<<endl;
		}
	}
	cout<<endl;
	#endif
	return 0;
}
LZ_DEF_compile_part_setfun_auto(main_SCP_4DrugExperiment);
#endif

int main(int argc,char *argv[]){
	LZ_DEF_run_mainfun(LZ_DEF_compile_part_val);
}