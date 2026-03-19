#ifndef LZ_DEF_hpp_OofA
#define LZ_DEF_hpp_OofA

#include<Eigen/LU>

//`#include<Eigen/Eigenvalues>
#include<Eigen/SparseCore>
#include"SymmetricGroup.hpp"
#include"temp220726_utility.hpp"
//`using namespace Liuze;
//`using namespace Liuze::math;

/****************************************************************************************************/

template<typename VecType>
Liuze::Type_Bool order_lex_less(VecType const & x0,VecType const & x1){
	Liuze::Type_Size L=std::min(x0.size(),x1.size()),i;
	for(i=0;i<L;++i){
		if(x0[i]<x1[i]) return true;
		if(x0[i]>x1[i]) return false;
	}
	return false;
}

template<typename IntType=Liuze::Type_UInt,typename RealType=Liuze::Type_Real,
	typename vec_permutation,
	typename=typename std::enable_if<std::is_integral<IntType>::value &&
	std::is_arithmetic<RealType>::value>::type>
RealType FiniteSymmetricGroup_distance_position
(vec_permutation const & perm0,vec_permutation const & perm1,RealType power=2){
	typedef IntType type_int;
	typedef RealType type_real;
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(perm0.size());
	vec_permutation pos0=SymGrp.inverse(perm0),pos1=SymGrp.inverse(perm1);
	Calculator_norm_L_stream<type_real,type_real,type_int> Cal_norm(power,true);
	type_int i_comp;
	for(i_comp=0;i_comp<perm0.size();++i_comp){
		Cal_norm<<(type_real(pos0[i_comp])-type_real(pos1[i_comp]));
	}
	return Cal_norm.Get_norm();
} //fun: FiniteSymmetricGroup_distance_position;

template<typename IntType=Liuze::Type_UInt,typename RealType=Liuze::Type_Real,
	typename vec_permutation,
	typename=typename std::enable_if<std::is_integral<IntType>::value &&
	std::is_arithmetic<RealType>::value>::type>
RealType FiniteSymmetricGroup_distance_permMr_oprL2norm
(vec_permutation const & perm0,vec_permutation const & perm1){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<type_real> type_matrix;
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(perm0.size());
	/*
	type_matrix diff=SymGrp.matrix_permutation_row(perm0)-SymGrp.matrix_permutation_row(perm1);
	diff=diff.transpose()*diff;
	diff=Eigen::SelfAdjointEigenSolver<type_matrix>(diff).eigenvalues();
	return type_real(sqrt(
		std::max(Liuze::math::abs(diff(0,0)),Liuze::math::abs(diff(diff.rows()-1,0)))
		));
	*/
	//*
	Eigen::VectorX<type_real> singular=
		Eigen::BDCSVD<type_matrix>(
		SymGrp.matrix_permutation_row(perm0)-SymGrp.matrix_permutation_row(perm1)
		).singularValues();
	return std::max(Liuze::math::abs(singular[0]),Liuze::math::abs(singular[singular.size()-1]));
	//*/
} //fun: FiniteSymmetricGroup_distance_permMr_oprL2norm;

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_PofC(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp<=type_int(1)) return mat_des;
	type_des res(n_run,n_comp);
	Liuze::math::FiniteSymmetricGroup<type_int,Liuze::Type_Real> SymGrp(n_comp);
	type_int i_run;
	for(i_run=0;i_run<n_run;++i_run){
		res.row(i_run)=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
	}
	return res;
} //fun: Get_PofC;

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> act_perm_OofADesign_left
(Liuze::Type_MatTemp<IntType> const & mat_perm,Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	Liuze::math::FiniteSymmetricGroup<type_int,Liuze::Type_Real> SymGrp(mat_des.cols());
	type_res res(mat_perm.rows()*mat_des.rows(),mat_des.cols());
	type_int i_perm,i_run,i_run0=0;
	for(i_perm=0;i_perm<mat_perm.rows();++i_perm){
		for(i_run=0;i_run<mat_des.rows();++i_run){
			res.row(i_run0)=SymGrp.compose(Eigen::RowVectorX<type_int>(mat_perm.row(i_perm)),
				Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
			++i_run0;
		}
	}
	return res;
} //fun: act_perm_OofADesign_left;
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> act_invperm_OofADesign_right
(Liuze::Type_MatTemp<IntType> const & mat_perm,Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	Liuze::math::FiniteSymmetricGroup<type_int,Liuze::Type_Real> SymGrp(mat_des.cols());
	type_res res(mat_perm.rows()*mat_des.rows(),mat_des.cols());
	type_int i_perm,i_run,i_run0=0;
	for(i_perm=0;i_perm<mat_perm.rows();++i_perm){
		for(i_run=0;i_run<mat_des.rows();++i_run){
			res.row(i_run0)=SymGrp.compose(Eigen::RowVectorX<type_int>(mat_des.row(i_run)),
				Eigen::RowVectorX<type_int>(mat_perm.row(i_perm)));
			++i_run0;
		}
	}
	return res;
} //fun: act_invperm_OofADesign_right;

/****************************************************************************************************/

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_OofAFullDesign(IntType const & n_component){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component==type_int(0)) return type_res(0,0);
	type_int i_run,n_run=1,i_comp;
	for(i_run=2;i_run<=n_component;++i_run) n_run*=i_run;
	type_res res(n_run,n_component);
	i_run=0;
	Liuze::math::Comb_enum_perm<type_int> iter_perm(n_component);
	while(i_run<n_run){
		for(i_comp=0;i_comp<n_component;++i_comp){
			res(i_run,i_comp)=(*iter_perm)[i_comp];
		}
		++iter_perm;
		++i_run;
	}
	return res;
} //fun: Get_OofAFullDesign;

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_OofAFullDesign_recursive
(IntType const & n_component){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component==type_int(0)) return type_res(0,0);
	type_int i_run,n_run=1,i_comp,k_comp,k_run,i_cols;
	for(i_run=2;i_run<=n_component;++i_run) n_run*=i_run;
	type_res res(n_run,n_component);
	//the design for 1 component:
	res(0,0)=type_int(0);
	k_run=1;
	for(k_comp=1;k_comp<n_component;){ //the design for (k_comp+1) components;
		i_run=0;
		res.block(i_run,k_comp,k_run,1).setConstant(k_comp);
		for(i_comp=1;i_comp<=k_comp;++i_comp){
			i_run+=k_run;
			i_cols=k_comp-i_comp;
			if(i_cols>type_int(0)){
				res.block(i_run,0,k_run,i_cols)=res.block(0,0,k_run,i_cols);
			}
			res.block(i_run,i_cols,k_run,1).setConstant(k_comp);
			res.block(i_run,i_cols+type_int(1),k_run,i_comp)=res.block(0,i_cols,k_run,i_comp);
		}
		k_run*=(++k_comp);
	}
	return res;
} //fun: Get_OofAFullDesign_recursive;

template<typename IntType,typename RandEngineType>
Liuze::Type_MatTemp<IntType> Get_OofARandDesign
(IntType const & n_run,IntType const & n_component,RandEngineType & rand_eng){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_run==type_int(0) || n_component==type_int(0)) return type_res(0,0);
	type_int i_run,i_comp;
	type_res res(n_run,n_component);
	Liuze::Type_ArrayTemp<type_int> perm(n_component);
	for(i_comp=0;i_comp<n_component;++i_comp) perm[i_comp]=i_comp;
	for(i_run=0;i_run<n_run;++i_run){
		std::shuffle(perm.begin(),perm.end(),rand_eng);
		for(i_comp=0;i_comp<n_component;++i_comp) res(i_run,i_comp)=perm[i_comp];
	}
	return res;
} //fun: Get_OofARandDesign;

template<typename FieldType,typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<typename FieldType::value_type>::value &&
	std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA_augrow
(IntType const & n_run,FieldType const & field,bool perm_col=true){
	typedef typename FieldType::value_type type_val;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_comp=field.size();
	if(n_comp<=type_int(0) || n_run<=type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return type_res::Zero(n_run,1);
	type_res res(n_run,n_comp);
	type_int ir,ic,ir0=0,i_part;
	type_int n_run_rest=n_run,n_run_part,nr0=n_comp;
	for(i_part=1;i_part<n_comp && n_run_rest>type_int(0);++i_part){
		n_run_part=std::min(nr0,n_run_rest);
		//a (n_comp-1)-cycle of {1,...,n_comp-1}:
		for(ic=0;ic<n_comp;++ic){
			res(ir0,ic)=type_int(field(i_part)*field(ic));
		}
		//composation of the above cycle and permutations generated by (1,...,n_comp-1,0):
		for(ir=1;ir<n_run_part;++ir){
			for(ic=0;ic<n_comp;++ic){
				res(ir0+ir,ic)=type_int(field(ir)+field(res(ir0,ic)));
			}
		}
		n_run_rest-=n_run_part;
		ir0+=n_run_part;
	}
	//permuting the last (n_comp-2) columns or components:
	if(n_run_rest>type_int(0) && n_comp>=type_int(4)){
		Liuze::math::Comb_enum_perm<type_int> iter_perm(n_comp-type_int(2));
		Liuze::Type_ArrayTemp<type_int> perm(n_comp);
		for(ic=0;ic<type_int(2);++ic) perm[ic]=ic;
		nr0=ir0;
		for(++iter_perm;n_run_rest>type_int(0) && iter_perm.Is_stat_deref();++iter_perm){
			n_run_part=std::min(nr0,n_run_rest);
			for(ic=type_int(2);ic<n_comp;++ic) perm[ic]=(*iter_perm)[ic-type_int(2)]+type_int(2);
			if(perm_col){
				for(ir=0;ir<n_run_part;++ir){
					for(ic=0;ic<n_comp;++ic) res(ir0+ir,ic)=res(ir,perm[ic]);
				}
			} else {
				for(ir=0;ir<n_run_part;++ir){
					for(ic=0;ic<n_comp;++ic) res(ir0+ir,ic)=perm[res(ir,ic)];
				}
			}
			n_run_rest-=n_run_part;
			ir0+=n_run_part;
		}
	}
	nr0=ir0;
	//repeating the full design:
	while(n_run_rest>type_int(0)){
		n_run_part=std::min(nr0,n_run_rest);
		res.block(ir0,0,n_run_part,n_comp)=res.block(0,0,n_run_part,n_comp);
		n_run_rest-=n_run_part;
		ir0+=n_run_part;
	}
	return res;
} //fun: Get_COA_augrow;

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA_addComponent_Latin
(Liuze::Type_MatTemp<IntType> const & design_init){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_run_0=design_init.rows();
	type_int n_comp=type_int(design_init.cols())+type_int(1);
	if(n_comp<=type_int(0) || n_run_0<=type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return type_res::Zero(1,1);
	type_int n_run=n_run_0*n_comp;
	type_res res(n_run,n_comp);
	type_int i_block,i_run=0,i_comp;
	for(i_block=0;i_block<n_run_0;++i_block){
		res.block(i_run,0,1,design_init.cols())=design_init.row(i_block);
		res(i_run,design_init.cols())=type_int(design_init.cols());
		++i_run;
		for(i_comp=1;i_comp<n_comp;++i_comp){
			res(i_run,design_init.cols())=res(i_run-type_int(1),0);
			res.block(i_run,0,1,design_init.cols())=
				res.block(i_run-type_int(1),1,1,design_init.cols());
			++i_run;
		}
	}
	return res;
} //fun: Get_COA_addComponent_Latin;

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA_addComponent_custom
(Liuze::Type_MatTemp<IntType> const & design_init,
Liuze::Type_MatTemp<IntType> const & augment_scheme){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_run_0=design_init.rows();
	type_int n_comp_0=design_init.cols();
	type_int n_comp=n_comp_0+type_int(1);
	type_int n_block=augment_scheme.rows();
	type_int n_col=augment_scheme.cols();
	if(n_block<=type_int(0) || n_block%n_comp!=type_int(0) || n_col>n_comp){
		return type_res(0,0);
	}
	if(n_comp<=type_int(0) || n_run_0<=type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return augment_scheme;
	type_int i_col,i_block,i_run,i_comp;
	type_res res(n_run_0*n_block,n_col);
	for(i_run=i_block=0;i_block<n_block;++i_block){
		for(i_col=0;i_col<n_col;++i_col){
			i_comp=augment_scheme(i_block,i_col);
			if(i_comp==n_comp_0){
				res.block(i_run,i_col,n_run_0,1)=type_res::Constant(n_run_0,1,n_comp_0);
			} else {
				res.block(i_run,i_col,n_run_0,1)=design_init.col(i_comp);
			}
		}
		i_run+=n_run_0;
	}
	return res;
} //fun: Get_COA_addComponent_custom;
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA_addComponent_custom
(Liuze::Type_MatTemp<IntType> const & design_init){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int ir,ic;
	type_int n_comp=design_init.cols()+type_int(1);
	type_res augment_scheme(n_comp,n_comp);
	for(ic=0;ic<n_comp;++ic) augment_scheme(0,ic)=ic;
	for(ir=1;ir<n_comp;++ir){
		augment_scheme(ir,n_comp-type_int(1))=augment_scheme(ir-type_int(1),0);
		for(ic=1;ic<n_comp;++ic){
			augment_scheme(ir,ic-type_int(1))=augment_scheme(ir-type_int(1),ic);
		}
	}
	return Get_COA_addComponent_custom<type_int>(design_init,augment_scheme);
} //fun: Get_COA_addComponent_custom;

template<typename FieldType,typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<typename FieldType::value_type>::value &&
	std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA3
(IntType const & n_run_init,FieldType const & field_init){
	typedef typename FieldType::value_type type_val;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_comp_init=field_init.size();
	if(n_comp_init<=type_int(1) || n_run_init<=type_int(0)) return type_res(0,0);
	type_int i_row,i_col;
	type_int n_comp=n_comp_init+type_int(1);
	type_res des_init=Get_COA_augrow(n_run_init,field_init,false);
	type_res augment_scheme(n_comp,n_comp);
	for(i_row=0;i_row<n_comp;++i_row) augment_scheme(i_row,i_row)=n_comp_init;
	for(i_col=1;i_col<n_comp;++i_col){
		augment_scheme(0,i_col)=type_int(field_init(i_col-type_int(1)));
	}
	augment_scheme(1,0)=type_int(field_init(0));
	for(i_col=2;i_col<n_comp;++i_col){
		augment_scheme(1,i_col)=type_int(field_init(i_col-type_int(1)).inv());
	}
	for(i_row=2;i_row<n_comp;++i_row){
		for(i_col=0;i_col<type_int(2);++i_col){
			augment_scheme(i_row,i_col)=type_int(field_init(i_col));
		}
		for(;i_col<n_comp;++i_col){
			if(i_col==i_row) continue;
			augment_scheme(i_row,i_col)=type_int(
				field_init(i_row-type_int(1))*
				Liuze::math::inverse(field_init(i_row-type_int(1))-field_init(i_col-type_int(1))));
		}
	}
	return Get_COA_addComponent_custom(des_init,augment_scheme);
} //fun: Get_COA3;

template<typename IntType>
Liuze::Type_MatTemp<IntType> Get_COA_halffull
(IntType const & n_comp){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_comp==type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return type_res::Zero(1,1);
	type_int i_run,n_run=1,i_comp;
	for(i_run=3;i_run<=n_comp;++i_run) n_run*=i_run;
	type_res res(n_run,n_comp);
	i_run=0;
	Liuze::math::Comb_enum_perm<type_int> iter_perm(n_comp);
	while(i_run<n_run){
		if(Liuze::math::FiniteSymmetricGroup<type_int,type_int>::weight_Cayley(*iter_perm,n_comp)
		%type_int(2)==type_int(0)){
			for(i_comp=0;i_comp<n_comp;++i_comp){
				res(i_run,i_comp)=(*iter_perm)[i_comp];
			}
			++i_run;
		}
		++iter_perm;
	}
	return res;
} //fun: Get_COA_halffull;

template<typename IntType>
Liuze::Type_MatTemp<IntType> Get_OofADesign_cycle
(IntType const & n_run,IntType const & n_component,bool perm_col=true){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component<=type_int(0) || n_run<=type_int(0)) return type_res(0,0);
	if(n_component==type_int(1)) return type_res::Zero(n_run,1);
	type_res res(n_run,n_component);
	type_int ir,ic,ir0=0,i_comp=n_component,i_part,ir00;
	type_int n_run_rest=n_run,n_run_part,nr0;
	Eigen::RowVectorX<type_int> perm(n_component);
	Liuze::math::FiniteSymmetricGroup<type_int,Liuze::Type_Real> SymGrp(n_component);
	//the first run:
	for(ic=0;ic<n_component;++ic) res(ir0,ic)=ic;
	--n_run_rest;
	++ir0;
	nr0=ir0;
	while(i_comp>type_int(1) && n_run_rest>type_int(0)){
		//2 for's: the cycle (1,2,...,i_comp-1,0, i_comp,...,n_component-1):
		for(ic=0;ic<i_comp;++ic) perm[ic]=Liuze::math::mod(ic+type_int(1),i_comp);
		for(;ic<n_component;++ic) perm[ic]=ic;
		ir00=0;
		//compose the previous block and perm^i_part:
		for(i_part=1;i_part<i_comp && n_run_rest>type_int(0);++i_part){
			n_run_part=std::min(n_run_rest,nr0); //the next block size;
			//compose the last block and perm:
			for(ir=0;ir<n_run_part;++ir){
				res.row(ir0+ir)=perm_col ?
					SymGrp.compose(Eigen::RowVectorX<type_int>(res.row(ir00+ir)),perm) :
					SymGrp.compose(perm,Eigen::RowVectorX<type_int>(res.row(ir00+ir)));
			}
			ir00+=nr0;
			n_run_rest-=n_run_part;
			ir0+=n_run_part;
		}
		nr0=ir0;
		--i_comp;
	}
	//repeating the full design:
	while(n_run_rest>type_int(0)){
		n_run_part=std::min(n_run_rest,nr0);
		res.block(ir0,0,n_run_part,n_component)=res.block(0,0,n_run_part,n_component);
		n_run_rest-=n_run_part;
		ir0+=n_run_part;
	}
	return res;
} //fun: Get_OofADesign_cycle;

template<typename RealType=Liuze::Type_Real,typename UIntType=Liuze::Type_UInt,
	int StorageMajor=Eigen::ColMajor>
Eigen::SparseMatrix<RealType> Get_MatTuple
(Liuze::Type_MatTemp<UIntType> const & mat_des,UIntType const & n_comp,UIntType strength){
	typedef RealType type_real;
	typedef UIntType type_int;
	typedef Eigen::SparseMatrix<RealType,StorageMajor> type_mat_sparse;
	type_int n_run=mat_des.rows(),n_col=mat_des.cols();
	if(n_run==type_int(0) || n_col==type_int(0) || n_comp<=type_int(0) || n_comp<n_col ||
	strength<=type_int(0) || strength>n_col){
		return type_mat_sparse(0,0);
	}
	type_int i_run,i_comb,i_comp;
	type_int n_comb=Liuze::math::Comb_num::comb<type_int>(n_col,strength);
	type_int n_perm=Liuze::math::Comb_num::perm<type_int>(n_comp,strength);
	type_int n_perm_1=Liuze::math::Comb_num::perm<type_int>(strength,strength);
	type_mat_sparse mat_tuple(n_run,n_comb*n_perm);
	Liuze::Type_ArrayTemp<Eigen::Triplet<type_real> > list_mat_tuple;
	list_mat_tuple.reserve(n_comb*n_run);
	Liuze::math::Comb_enum_comb<type_int> iter_comb(n_col,strength);
	Liuze::Type_ArrayTemp<type_int> run_part(strength),run_part_rank(strength);
	Liuze::Type_ArrayTemp<type_int> run_part_comb(strength),run_part_perm(strength);
	auto fun_less_run_part=[&run_part](type_int const & x0,type_int const & x1)->bool{
			return run_part[x0]<run_part[x1];
		};
	for(i_comb=0;i_comb<n_comb;++i_comb){
		for(i_run=0;i_run<n_run;++i_run){
			for(i_comp=0;i_comp<strength;++i_comp){
				run_part[i_comp]=mat_des(i_run,(*iter_comb)[i_comp]);
				run_part_rank[i_comp]=i_comp;
			}
			std::sort(run_part_rank.begin(),run_part_rank.end(),fun_less_run_part);
			for(i_comp=0;i_comp<strength;++i_comp){
				run_part_comb[i_comp]=run_part[run_part_rank[i_comp]];
				run_part_perm[run_part_rank[i_comp]]=i_comp;
			}
			i_comp=i_comb*n_perm
				+Comb_enum_comb_index<type_int,type_int>(n_comp,run_part_comb)*n_perm_1
				+Comb_enum_perm_index<type_int,type_int>(run_part_perm);
			list_mat_tuple.push_back(Eigen::Triplet<type_real>(i_run,i_comp,type_real(1)));
		}
		++iter_comb;
	}
	mat_tuple.setFromTriplets(list_mat_tuple.begin(),list_mat_tuple.end());
	return mat_tuple;
} //fun: Get_MatTuple;
template<typename RealType=Liuze::Type_Real,typename UIntType=Liuze::Type_UInt,
	int StorageMajor=Eigen::ColMajor>
Eigen::SparseMatrix<RealType> Get_MatTuple
(Liuze::Type_MatTemp<UIntType> const & mat_des,UIntType strength=2){
	return Get_MatTuple<RealType,UIntType,StorageMajor>(mat_des,mat_des.cols(),strength);
} //fun: Get_MatTuple;

template<typename UIntType=Liuze::Type_UInt>
Liuze::Type_Bool Is_COA
(Liuze::Type_MatTemp<UIntType> const & mat_des,UIntType strength=2,UIntType * ptr_index=NULL){
	typedef UIntType type_int;
	typedef Eigen::SparseMatrix<type_int,Eigen::ColMajor> type_mat_sparse;
	type_int n_row=mat_des.rows(),n_col=mat_des.cols();
	if(ptr_index) *ptr_index=0;
	if(strength<type_int(0) || strength>n_col) return false;
	if(n_col==type_int(0)) return true;
	if(strength==type_int(0)) return true;
	if(n_row==type_int(0)) return false;
	type_mat_sparse mat_tuple=Get_MatTuple<type_int,type_int,Eigen::ColMajor>(mat_des,strength);
	n_col=mat_tuple.cols();
	if(n_col==type_int(0)) return false;
	type_int i_col;
	n_row=mat_tuple.col(0).nonZeros();
	if(n_row<=type_int(0)) return false;
	for(i_col=1;i_col<n_col;++i_col){
		if(type_int(mat_tuple.col(i_col).nonZeros())!=n_row) return false;
	}
	if(ptr_index) *ptr_index=n_row;
	return true;
} //fun: Is_COA;

template<typename UIntType=Liuze::Type_UInt>
Liuze::Type_Bool Is_OAtypeI
(Liuze::Type_MatTemp<UIntType> const & mat_des,UIntType const & n_comp,UIntType strength){
	typedef UIntType type_int;
	type_int n_run=mat_des.rows(),n_col=mat_des.cols();
	if(n_run==type_int(0) || n_col==type_int(0) || n_comp<=type_int(0) || n_comp<n_col ||
	strength<type_int(0) || strength>n_col){
		return false;
	}
	if(strength==type_int(0)) return true;
	type_int i_run,i_comp;
	type_int n_perm=Liuze::math::Comb_num::perm<type_int>(n_comp,strength);
	if(n_run%n_perm!=type_int(0)) return false;
	type_int index=n_run/n_perm;
	Liuze::Type_ArrayTemp<type_int> list_tuple(n_perm);
	type_int n_perm_1=Liuze::math::Comb_num::perm<type_int>(strength,strength);
	Liuze::math::Comb_enum_comb<type_int> iter_comb(n_col,strength);
	Liuze::Type_ArrayTemp<type_int> run_part(strength),run_part_rank(strength);
	Liuze::Type_ArrayTemp<type_int> run_part_comb(strength),run_part_perm(strength);
	auto fun_less_run_part=[&run_part](type_int const & x0,type_int const & x1)->bool{
			return run_part[x0]<run_part[x1];
		};
	while(iter_comb){
		list_tuple.assign(n_perm,type_int(0));
		for(i_run=0;i_run<n_run;++i_run){
			for(i_comp=0;i_comp<strength;++i_comp){
				run_part[i_comp]=mat_des(i_run,(*iter_comb)[i_comp]);
				run_part_rank[i_comp]=i_comp;
			}
			std::sort(run_part_rank.begin(),run_part_rank.end(),fun_less_run_part);
			for(i_comp=0;i_comp<strength;++i_comp){
				run_part_comb[i_comp]=run_part[run_part_rank[i_comp]];
				run_part_perm[run_part_rank[i_comp]]=i_comp;
			}
			i_comp=Comb_enum_comb_index<type_int,type_int>(n_comp,run_part_comb)*n_perm_1
				+Comb_enum_perm_index<type_int,type_int>(run_part_perm);
			++(list_tuple[i_comp]);
			if(list_tuple[i_comp]>index) return false;
		}
		++iter_comb;
	}
	return true;
} //fun: Is_OAtypeI;

template<typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<IntType> Get_PWOOD_even(IntType const & n_component){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component<type_int(4) || n_component%type_int(2)!=type_int(0)){
		return type_res(0,0);
	}
	type_int i0_run,i_run,i_comp,ti0,ti1;
	type_int n_comp_half=n_component/type_int(2);
	type_int n_run=n_component;
	for(i_run=n_component-1;i_run>n_comp_half;--i_run) n_run*=i_run;
	type_int n_run_halfgroup=1;
	//n_run_halfgroup=frac_prod(n_comp_half):
	for(i_run=2;i_run<=n_comp_half;++i_run) n_run_halfgroup*=i_run;
	type_res res(n_run,n_component);
	Liuze::math::Comb_enum_comb<type_int> iter_comb(n_component,n_comp_half);
	typename std::remove_cv<typename Liuze::math::Comb_enum_comb<type_int>::value_type>::type
		iter_comb_rest(n_comp_half);
	Liuze::math::Comb_enum_perm<type_int> iter_perm(n_comp_half);
	i0_run=0;
	while((*iter_comb)[0]==type_int(0)){
		//filling iter_comb_rest{
		i_comp=0;
		ti1=1;
		ti0=1;
		while(i_comp<n_comp_half){
			while(i_comp<n_comp_half && (ti1>=n_comp_half || ti0<(*iter_comb)[ti1])){
				iter_comb_rest[i_comp]=ti0;
				++i_comp;
				++ti0;
			}
			++ti0;
			++ti1;
		}
		//}(filling iter_comb_rest)
		iter_perm.Set_begin();
		i_run=0;
		while(i_run<n_run_halfgroup){
			ti0=i0_run+n_run_halfgroup;
			//`ti1=i0_run-i_run+n_run_halfgroup+n_run_halfgroup-type_int(1)-i_run;
			for(i_comp=0;i_comp<n_comp_half;++i_comp){
				res(i0_run,i_comp)=(*iter_comb)[(*iter_perm)[i_comp]];
				res(i0_run,n_comp_half+i_comp)=iter_comb_rest[(*iter_perm)[i_comp]];
				res(ti0,n_comp_half+i_comp)=res(i0_run,i_comp);
				//`res(ti1,i_comp)=res(i0_run,n_comp_half+i_comp);
				res(ti0,n_comp_half-type_int(1)-i_comp)=res(i0_run,n_comp_half+i_comp);
			}
			++i_run;
			++i0_run;
			++iter_perm;
		}
		i0_run+=n_run_halfgroup;
		++iter_comb;
	}
	return res;
} //fun: Get_PWOOD_even;

template<typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<IntType> Get_PWOOD_odd(IntType const & n_component){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component<type_int(5) || n_component%type_int(2)==type_int(0)){
		return type_res(0,0);
	}
	type_int i_comp,i_run,i_col;
	type_int c_component=n_component-type_int(1);
	type_res des_sub=Get_PWOOD_even(c_component);
	type_int n_run=n_component*des_sub.rows();
	type_res res(n_run,n_component);
	for(i_col=0;i_col<n_component;++i_col){
		i_run=i_col*des_sub.rows();
		for(i_comp=0;i_comp<i_col;++i_comp){
			res.block(i_run,i_comp,des_sub.rows(),1)=des_sub.col(i_comp);
		}
		res.block(i_run,i_comp,des_sub.rows(),1).setConstant(c_component);
		for(;i_comp<c_component;++i_comp){
			res.block(i_run,i_comp+1,des_sub.rows(),1)=des_sub.col(i_comp);
		}
	}
	return res;
} //fun: Get_PWOOD_odd;

template<typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<IntType> Get_PWOOD(IntType const & n_component){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<IntType> type_res;
	if(n_component<type_int(4)) return Get_OofAFullDesign(n_component);
	if(n_component%type_int(2)==type_int(0)) return Get_PWOOD_even(n_component);
	else return Get_PWOOD_odd(n_component);
} //fun: Get_PWOOD;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<RealType> Get_MatModel_PWO(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp==type_int(0)) return type_matrix(0,0);
	type_int n_fac=1+n_comp*(n_comp-type_int(1))/type_int(2);
	type_matrix mat_model(n_run,n_fac);
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(n_comp);
	Eigen::RowVectorX<type_int> run_inv;
	type_int i_run,i_comp,j_comp,i_fac;
	for(i_run=0;i_run<n_run;++i_run){
		run_inv=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
		mat_model(i_run,0)=type_real(1);
		i_fac=1;
		for(i_comp=0;i_comp<n_comp;++i_comp){
			for(j_comp=i_comp+1;j_comp<n_comp;++j_comp){
				mat_model(i_run,i_fac)= run_inv[i_comp]<run_inv[j_comp] ? type_real(1) : -type_real(1);
				++i_fac;
			}
		}
	}
	return mat_model;
} //fun: Get_MatModel_PWO;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
RealType Get_Eff_PWO_D(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_real r_n_run=mat_des.rows(),r_n_comp=mat_des.cols();
	if(r_n_comp==type_real(0)) return type_real(1);
	if(r_n_run==type_real(0)) return type_real(0);
	type_matrix mat_model=Get_MatModel_PWO(mat_des);
	type_real r_n_fac=mat_model.cols();
	type_real eff=type_real(pow(r_n_comp+type_real(1),r_n_comp-type_real(1)))/
		type_real(pow(type_real(3),r_n_fac-type_real(1)));
	return type_real(pow((mat_model.transpose()*mat_model/r_n_run).determinant()/eff,
		type_real(1)/r_n_fac));
} // fun: Get_Eff_PWO_D;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<RealType> Get_MatModel_CP_baseline(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp==type_int(0)) return type_matrix(0,0);
	type_int n_fac=n_comp-type_int(1);
	n_fac=n_fac*n_fac+type_int(1);
	type_matrix mat_model=type_matrix::Zero(n_run,n_fac);
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(n_comp);
	Eigen::RowVectorX<type_int> run_inv;
	type_int i_run,i_comp,k_fac=n_comp-type_int(1);
	for(i_run=0;i_run<n_run;++i_run){
		run_inv=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
		mat_model(i_run,0)=type_real(1);
		for(i_comp=1;i_comp<n_comp;++i_comp){
			if(run_inv[i_comp]==type_int(0)) continue;
			mat_model(i_run,(i_comp-1)*k_fac+run_inv[i_comp])=type_real(1);
		}
	}
	return mat_model;
} //fun: Get_MatModel_CP_baseline;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
RealType Get_Eff_CP_D(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_real r_n_run=mat_des.rows(),r_n_comp=mat_des.cols();
	if(r_n_comp==type_real(0)) return type_real(1);
	if(r_n_run==type_real(0)) return type_real(0);
	type_matrix mat_model=Get_MatModel_CP_baseline(mat_des);
	type_real eff=r_n_comp-type_real(2);
	eff*=eff;
	eff=(r_n_comp-type_real(1))*type_real(pow(r_n_comp,eff));
	return r_n_comp*(r_n_comp-type_real(1))*
		type_real(pow((mat_model.transpose()*mat_model/r_n_run).determinant()/eff,
		type_real(1)/type_real(mat_model.cols())));
} // fun: Get_Eff_CP_D;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<RealType> Get_MatModel_CPP_deg1_baseline
(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp==type_int(0)) return type_matrix(0,0);
	type_matrix mat_model=type_matrix::Zero(n_run,n_comp);
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(n_comp);
	Eigen::RowVectorX<type_int> run_pos;
	type_int i_run,i_comp;
	for(i_run=0;i_run<n_run;++i_run){
		run_pos=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
		mat_model(i_run,0)=type_real(1);
		for(i_comp=1;i_comp<n_comp;++i_comp){
			mat_model(i_run,i_comp)=type_real(run_pos[i_comp]);
		}
	}
	return mat_model;
} //fun: Get_MatModel_CPP_deg1_baseline;
template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<RealType> Get_MatModel_CPP_deg2indep_baseline
(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp==type_int(0)) return type_matrix(0,0);
	type_matrix mat_model=type_matrix::Zero(n_run,2*n_comp-1);
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(n_comp);
	Eigen::RowVectorX<type_int> run_pos;
	type_int i_run,i_comp;
	for(i_run=0;i_run<n_run;++i_run){
		run_pos=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
		mat_model(i_run,0)=type_real(1);
		for(i_comp=1;i_comp<n_comp;++i_comp){
			mat_model(i_run,i_comp)=type_real(run_pos[i_comp]);
			mat_model(i_run,n_comp-1+i_comp)=mat_model(i_run,i_comp)*mat_model(i_run,i_comp);
		}
	}
	return mat_model;
} //fun: Get_MatModel_CPP_deg2indep_baseline;
template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt>
Liuze::Type_MatTemp<RealType> Get_MatModel_CPP_deg2_baseline
(Liuze::Type_MatTemp<IntType> const & mat_des){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	typedef Liuze::Type_MatTemp<RealType> type_matrix;
	type_int n_run=mat_des.rows(),n_comp=mat_des.cols();
	if(n_run==type_int(0) || n_comp==type_int(0)) return type_matrix(0,0);
	type_matrix mat_model=type_matrix::Zero(n_run,(n_comp-type_int(1))*(n_comp+type_int(2))/type_int(2));
	Liuze::math::FiniteSymmetricGroup<type_int,type_real> SymGrp(n_comp);
	Eigen::RowVectorX<type_int> run_pos;
	type_int i_run,i_comp,i_fac,j_comp;
	for(i_run=0;i_run<n_run;++i_run){
		run_pos=SymGrp.inverse(Eigen::RowVectorX<type_int>(mat_des.row(i_run)));
		mat_model(i_run,0)=type_real(1);
		i_fac=n_comp-type_int(2);
		for(i_comp=1;i_comp<n_comp;++i_comp){
			mat_model(i_run,i_comp)=type_real(run_pos[i_comp]);
			if(i_comp==type_int(1)) continue;
			mat_model(i_run,i_fac+i_comp)=mat_model(i_run,i_comp)*mat_model(i_run,i_comp);
		}
		i_fac=type_int(2)*n_comp-type_int(2);
		for(i_comp=1;i_comp<n_comp;++i_comp){
			for(j_comp=i_comp+type_int(1);j_comp<n_comp;++j_comp){
				mat_model(i_run,i_fac)=mat_model(i_run,i_comp)*mat_model(i_run,j_comp);
				++i_fac;
			}
		}
	}
	return mat_model;
} //fun: Get_MatModel_CPP_deg2_baseline;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt,
	typename CritType=std::function<RealType(Liuze::Type_MatTemp<IntType> const &)> >
Liuze::Type_MatTemp<IntType> Get_OofAOptDesign_bubbleEx
(Liuze::Type_MatTemp<IntType> des_init,CritType & loss,RealType * ptr_val_crit=NULL,IntType n_round=100){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	type_int n_run=des_init.rows(),n_comp=des_init.cols();
	if(n_run==type_int(0) || n_comp<=type_int(1)){
		if(ptr_val_crit) *ptr_val_crit=loss(des_init);
		return des_init;
	}
	if(n_round<=type_int(0)) n_round=type_int(100);
	type_int i_round,i_run,i_comp,i0_comp;
	type_int ti_swap;
	type_real crit,crit0=loss(des_init);
	for(i_round=0;i_round<n_round;++i_round){
		for(i_run=0;i_run<n_run;++i_run){
			for(i0_comp=n_comp-type_int(1);i0_comp>type_int(0);--i0_comp){
				for(i_comp=0;i_comp<i0_comp;++i_comp){
					ti_swap=des_init(i_run,i_comp);
					des_init(i_run,i_comp)=des_init(i_run,i_comp+1);
					des_init(i_run,i_comp+1)=ti_swap;
					crit=loss(des_init);
					if(crit<crit0){
						crit0=crit;
					} else if(crit>crit0){
						des_init(i_run,i_comp+1)=des_init(i_run,i_comp);
						des_init(i_run,i_comp)=ti_swap;
					}
				}
			}
		}
	}
	if(ptr_val_crit) *ptr_val_crit=crit0;
	return des_init;
} //fun: Get_OofAOptDesign_bubbleEx;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt,
	typename CritType=std::function<RealType(Liuze::Type_MatTemp<IntType> const &)> >
Liuze::Type_MatTemp<IntType> Get_COA_improve_bubbleEx
(Liuze::Type_MatTemp<IntType> COA_init,CritType & loss,RealType * ptr_val_crit=NULL){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	type_int n_run=COA_init.rows(),n_comp=COA_init.cols();
	if(n_run==type_int(0) || n_comp<=type_int(3)){
		if(ptr_val_crit) *ptr_val_crit=loss(COA_init);
		return COA_init;
	}
	type_int i_comp,i0_comp;
	Eigen::VectorX<type_int> tv_swap;
	type_real crit,crit0=loss(COA_init);
	for(i0_comp=n_comp-type_int(3);i0_comp>type_int(0);--i0_comp){
		for(i_comp=0;i_comp<i0_comp;++i_comp){
			tv_swap=COA_init.col(i_comp);
			COA_init.col(i_comp)=COA_init.col(i_comp+1);
			COA_init.col(i_comp+1)=tv_swap;
			crit=loss(COA_init);
			if(crit<crit0){
				crit0=crit;
			} else if(crit>crit0){
				COA_init.col(i_comp+1)=COA_init.col(i_comp);
				COA_init.col(i_comp)=tv_swap;
			}
		}
	}
	if(ptr_val_crit) *ptr_val_crit=crit0;
	return COA_init;
} //fun: Get_COA_improve_bubbleEx;

template<typename RealType=Liuze::Type_Real,typename IntType=Liuze::Type_UInt,
	typename CritType=std::function<RealType(Liuze::Type_MatTemp<IntType> const &)>,
	typename IterShiftType,
	typename=typename std::enable_if<
	std::is_integral<IntType>::value &&
	std::is_convertible<typename IterShiftType::value_type,Liuze::Type_MatTemp<IntType> >::value
	>::type>
Liuze::Type_MatTemp<IntType> Get_OofAOpt_shiftBubbleEx
(CritType & loss,Liuze::Type_MatTemp<IntType> const & des,IterShiftType iter_shift_begin,
IntType const & n_shift,Liuze::Type_Bool perm_col=true,
IntType n_bubbleEx=1,RealType * ptr_val_loss=NULL){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<IntType> type_des;
	type_int n_comp=des.cols();
	if(n_comp==type_int(0) || n_shift<=type_int(0) || n_bubbleEx<type_int(0)){
		return type_des(0,0);
	}
	if(n_comp==type_int(1)) return type_des::Ones(1,1);
	type_int i_shift,i_run,i_comp,i_comp_0,i_ex;
	type_int ti_swap;
	type_des sol(n_comp,1),sol_local(n_comp,1);
	type_real val_loss,val_loss_sol,val_loss_sol_local;
	type_int n_run_des=des.rows();
	Liuze::Type_Bool tag_n_run_des=true,tag_bubbleEx;
	if(n_run_des==type_int(0)){
		n_run_des=type_int(1);
		tag_n_run_des=false;
	}
	type_des des_shift(n_run_des,n_comp);
	for(i_shift=0;i_shift<n_shift;++i_shift){
		if(!tag_n_run_des){
			des_shift=*iter_shift_begin;
		} else if(perm_col){
			for(i_run=0;i_run<n_run_des;++i_run){
				des_shift.row(i_run)=
					Liuze::math::FiniteSymmetricGroup<type_int,type_real>::compose
					(Eigen::RowVectorX<type_int>(des.row(i_run)),
					Eigen::RowVectorX<type_int>(*iter_shift_begin));
			}
		} else {
			for(i_run=0;i_run<n_run_des;++i_run){
				des_shift.row(i_run)=
					Liuze::math::FiniteSymmetricGroup<type_int,type_real>::compose
					(Eigen::RowVectorX<type_int>(*iter_shift_begin),
					Eigen::RowVectorX<type_int>(des.row(i_run)));
			}
		}
		for(i_run=0;i_run<n_run_des;++i_run){
			val_loss=loss(des_shift.row(i_run).transpose());
			if(i_run==type_int(0) || val_loss<val_loss_sol_local){
				i_comp=i_run;
				val_loss_sol_local=val_loss;
			}
		}
		sol_local=des_shift.row(i_comp).transpose();
		for(i_ex=0;i_ex<n_bubbleEx;++i_ex){
			tag_bubbleEx=false;
			for(i_comp=n_comp-type_int(1);i_comp>type_int(0);--i_comp){
				for(i_comp_0=0;i_comp_0<i_comp;++i_comp_0){
					ti_swap=sol_local(i_comp_0,0);
					sol_local(i_comp_0,0)=sol_local(i_comp_0+type_int(1),0);
					sol_local(i_comp_0+type_int(1),0)=ti_swap;
					val_loss=loss(sol_local);
					if(val_loss<val_loss_sol_local){
						tag_bubbleEx=true;
						val_loss_sol_local=val_loss;
					} else {
						sol_local(i_comp_0+type_int(1),0)=sol_local(i_comp_0,0);
						sol_local(i_comp_0,0)=ti_swap;
					}
				}
			}
			if(!tag_bubbleEx) break;
		}
		if(i_shift==type_int(0) || val_loss_sol_local<val_loss_sol){
			val_loss_sol=val_loss_sol_local;
			sol=sol_local;
		}
		++iter_shift_begin;
	}
	if(ptr_val_loss) *ptr_val_loss=val_loss_sol;
	return sol;
} //fun: Get_OofAOpt_shiftBubbleEx;

/****************************************************************************************************/

#endif //#ifndef LZ_DEF_hpp_OofA