#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_Real_Vec_CPP
#define LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_Real_Vec_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//Discrete distribution for real vectors:
//class Distr_Discrete_Real_Vec{
//public:
	//Constructor:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::Distr_Discrete_Real_Vec
	(void):
	type_base(1,type_this::stat_category_discrete,(unsigned char)15),
	supp(1,type_res::Zero(1,1)),distr_index(1),pre(LZ_DEF_const_default_precision){
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::Distr_Discrete_Real_Vec
	(type_res_seq const & support):
	type_base(1,type_this::stat_category_discrete,(unsigned char)15),
	supp(support.size()==0 ? 1 : support.size()),distr_index(support.size()),
	pre(LZ_DEF_const_default_precision){
		if(support.size()==0){
			supp[0]=type_res::Zero(this->dim,1);
		} else {
			Type_UInt j=0,i;
			typename type_res_seq::iterator ii=support.begin();
			for(;ii<support.end();++ii){
				if(ii->cols()>0 && ii->rows()>j) j=ii->rows();
			}
			if(j==0){
				for(i=0;i<supp.size();++i){
					supp[i]=(type_real)i*type_res::Ones(1,1);
				}
			} else {
				this->dim=j;
				for(i=0;i<supp.size();++i){
					if(support[i].cols()==0){
						supp[i]=type_res::Zero(this->dim,1);
					} else {
						supp[i].resize(this->dim,1);
						for(j=0;j<support[i].rows();++j) supp[i](j,0)=support[i](j,0);
						for(;j<this->dim;++j) supp[i](j,0)=(type_real)0;
					}
				}
			}
		}
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::Distr_Discrete_Real_Vec
	(type_res_seq const & support,type_real_prob_seq const & prob):
	type_base(1,type_this::stat_category_discrete,(unsigned char)15),
	supp(),distr_index(),pre(LZ_DEF_const_default_precision){
		if(support.size()==0){
			supp.resize(1);
			supp[0]=type_res::Zero(1,1);
		} else {
			Type_UInt j=0,i;
			supp.resize(prob.size()==0 ? support.size() : std::min(support.size(),prob.size()));
			for(i=0;i<supp.size();++i){
				if(support[i].cols()>0 && support[i].rows()>j) j=support[i].rows();
			}
			if(j==0){
				for(i=0;i<supp.size();++i){
					supp[i]=(type_real)i*type_res::Ones(1,1);
				}
			} else {
				this->dim=j;
				for(i=0;i<supp.size();++i){
					if(support[i].cols()==0){
						supp[i]=type_res::Zero(this->dim,1);
					} else {
						supp[i].resize(this->dim,1);
						for(j=0;j<support[i].rows();++j) supp[i](j,0)=support[i](j,0);
						for(;j<this->dim;++j) supp[i](j,0)=(type_real)0;
					}
				}
			}
			if(prob.size()==0 ||
				distr_index.Set_prob(std::vector<type_real_prob>(prob.begin(),prob.begin()+supp.size()))!=
				(typename type_distr_index::type_stat)0){
				distr_index.Set_prob(supp.size());
			}
		}
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::Distr_Discrete_Real_Vec
	(type_this const & distr):
	type_base((type_base const &)distr),supp(distr.supp),distr_index(distr.distr_index),pre(distr.pre){
	}
	
	//Destructor:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::~Distr_Discrete_Real_Vec
	(void){
	}
	
	//operator:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_this&
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::operator =
	(type_this const & distr){
		if(&distr==this) return *this;
		((type_base*)this)->operator =((type_base const &)distr);
		supp=distr.supp;
		distr_index=distr.distr_index;
		return *this;
	}
	
	//virtual from type_base:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	inline
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_res
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::rand
	(type_rand_engine & engine) const {
		return supp[distr_index.rand(engine)];
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_res_seq
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::sample
	(Type_UInt const num,type_rand_engine & engine) const {
		if(num==0) return type_res_seq();
		auto res_index=distr_index.sample(num,engine);
		Type_UInt i;
		type_res_seq res(num);
		for(i=0;i<num;++i) res[i]=supp[res_index[i]];
		return res;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::PDF
	(type_res const & x) const {
		if(!this->Is_type_res(x)) return (type_real_prob)0;
		Type_UInt i;
		for(i=0;i<this->Get_size();++i){
			if(x.isApprox(supp[i],pre)) return distr_index.Get_prob(i);
		}
		return (type_real_prob)0;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::CDF
	(Type_MatTemp<type_real> const & x) const {
		if(x.rows()!=this->dim || x.cols()!=1) return (type_real_prob)0;
		Type_UInt i,j;
		type_real_prob res=0;
		//CDF_def: the definition of CDF: 0: P(X<x); 1: P(X<=x):
		if(distr_index.Get_CDF_def()==(typename type_distr_index::type_stat)0){
			for(i=0;i<this->Get_size();++i){
				for(j=0;j<this->dim && supp[i](j,0)<x(j,0);++j);
				if(j==this->dim) res+=this->Get_prob(i);
			}
		} else {
			for(i=0;i<this->Get_size();++i){
				for(j=0;j<this->dim && supp[i](j,0)<=x(j,0);++j);
				if(j==this->dim) res+=this->Get_prob(i);
			}
		}
		return res;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::type_comp
	Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType>::CF
	(Type_MatTemp<type_real> const & x) const {
		type_comp res((typename type_comp::value_type)0,(typename type_comp::value_type)0);
		if(x.rows()!=this->dim || x.cols()!=1) return res;
		Type_UInt i;
		for(i=0;i<this->Get_size();++i){
			res+=exp(type_comp((type_real)0,(x.transpose()*supp[i])(0,0)))*(this->Get_prob(i));
		}
		return res;
	}
//}(class Distr_Discrete_Real_Vec)
//(Discrete distribution for real vectors)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_Real_Vec_CPP
#endif //#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L