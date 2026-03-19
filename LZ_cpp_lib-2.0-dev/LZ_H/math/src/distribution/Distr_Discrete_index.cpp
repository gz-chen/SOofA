#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_index_CPP
#define LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_index_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//General discrete distribution through index:
//class Distr_Discrete_index{
//public:
	//Constructor:
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Distr_Discrete_index
	(void):
	type_base((type_stat)0),uniform((type_real_prob)0,(type_real_prob)1),
	size(1),pdf(1,(type_real_prob)1),cdf(1,(type_real_prob)1),CDF_def((type_stat)0){
		this->stat_able=(type_stat)31;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Distr_Discrete_index
	(type_res const & number):
	type_base((type_stat)0),uniform((type_real_prob)0,(type_real_prob)1),
	size(number>(type_res)0 ? number : (type_res)1),pdf(size),cdf(size),CDF_def((type_stat)0){
		type_res i;
		type_real_prob prob=(type_real_prob)1/(type_real_prob)size;
		cdf[0]=pdf[0]=prob;
		for(i=1;i<size;++i){
			pdf[i]=prob;
			cdf[i]=cdf[i-1]+prob;
		}
		this->stat_able=(type_stat)31;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Distr_Discrete_index
	(type_real_prob_seq const & prob):
	type_base((type_stat)0),uniform((type_real_prob)0,(type_real_prob)1),
	size(prob.size()==0 ? (type_res)1 : prob.size()),pdf(size),cdf(size),CDF_def((type_stat)0){
		if(prob.size()==0){
			cdf[0]=pdf[0]=(type_real_prob)1;
		} else {
			type_res i;
			type_real_prob pdf_total=0;
			for(i=0;i<size;++i){
				pdf[i]= prob[i]>(type_real_prob)0 ? prob[i] : (type_real_prob)0;
				pdf_total+=pdf[i];
			}
			if(pdf_total!=(type_real_prob)1){
				if(pdf_total==(type_real_prob)0){
					std::fill_n(std::begin(pdf),size,(type_real_prob)1/(type_real_prob)size);
				} else {
					for(i=0;i<size;++i) pdf[i]/=pdf_total; //normalising pdf;
				}
			}
			std::partial_sum(std::begin(pdf),std::end(pdf),std::begin(cdf)); //calculate cdf;
		}
		this->stat_able=(type_stat)31;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Distr_Discrete_index
	(type_this const & distr):
	type_base((type_base&)distr),uniform((type_real_prob)0,(type_real_prob)1),
	size(distr.size),pdf(distr.pdf),cdf(distr.cdf),CDF_def(distr.CDF_def){
	}
	
	//operator:
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_this&
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::operator =
	(type_this const & distr){
		if(&distr==this) return *this;
		((type_base*)this)->operator =((type_base&)distr);
		size=distr.size;
		pdf=distr.pdf;
		cdf=distr.cdf;
		CDF_def=distr.CDF_def;
		return *this;
	}
	
	//virtual from type_base:
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_res
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::rand
	(type_rand_engine & engine) const {
		type_real_prob prob=uniform(engine);
		//the following is the same as this->quantile(prob).
		if(CDF_def==(type_stat)0){
			if(size==1 || prob<cdf[0]) return (type_res)0;
			if(prob>=cdf[size-2]) return size-1;
			Type_UInt iL=0,iR=size-2,iM;
			while(iL+1<iR){
				iM=(iL+iR)/2;
				if(prob<cdf[iM]) iR=iM;
				else iL=iM;
			}
			return iR;
		} else {
			if(size==1 || prob<=cdf[0]) return (type_res)0;
			if(prob>cdf[size-2]) return size-1;
			Type_UInt iL=0,iR=size-2,iM;
			while(iL+1<iR){
				iM=(iL+iR)/2;
				if(prob<=cdf[iM]) iR=iM;
				else iL=iM;
			}
			return iR;
		}
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_res_seq
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::sample
	(Type_UInt const num,type_rand_engine & engine) const {
		if(num==0) return type_res_seq();
		if(size==1) return type_res_seq(num,(type_res)0);
		type_res_seq res(num);
		Type_UInt i;
		type_real_prob prob;
		//the following is similar to this->quantile(prob).
		if(CDF_def==(type_stat)0){
			for(i=0;i<num;++i){
				prob=uniform(engine);
				if(prob<cdf[0]) res[i]=(type_res)0;
				if(prob>=cdf[size-2]) res[i]=size-1;
				Type_UInt iL=0,iR=size-2,iM;
				while(iL+1<iR){
					iM=(iL+iR)/2;
					if(prob<cdf[iM]) iR=iM;
					else iL=iM;
				}
				res[i]=(type_res)iR;
			}
		} else {
			for(i=0;i<num;++i){
				prob=uniform(engine);
				if(prob<=cdf[0]) res[i]=(type_res)0;
				if(prob>cdf[size-2]) res[i]=size-1;
				Type_UInt iL=0,iR=size-2,iM;
				while(iL+1<iR){
					iM=(iL+iR)/2;
					if(prob<=cdf[iM]) iR=iM;
					else iL=iM;
				}
				res[i]=(type_res)iR;
			}
		}
		return res;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	inline
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::PDF
	(type_res const & x) const {
		return (x>=(type_res)0 && x<size) ? pdf[x] : (type_real_prob)0;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::CDF
	(type_real const & x) const {
		if(x<(type_real)0) return (type_real_prob)0;
		if(x>(type_real)(size-1)) return (type_real_prob)1;
		Type_UInt i=(Type_UInt)floor(x);
		if(!Is_integer(x)) return cdf[i];
		if(CDF_def==(type_stat)1) return cdf[i];
		return i==0 ? (type_real_prob)0 : cdf[i-1];
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_comp
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::CF
	(type_real const & x) const {
		type_comp res((type_real)0,(type_real)0);
		Type_UInt i;
		for(i=0;i<size;++i){
			res+=pdf[i]*exp(type_comp((type_real)0,x*(type_real)i));
		}
		return res;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_res
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::quantile
	(type_real_prob const & prob) const {
		//When CDF is P(X<x), the quantile is "[0)[1)...[size-1)":
		if(CDF_def==(type_stat)0){
			if(size==1 || prob<cdf[0]) return (type_res)0;
			if(prob>=cdf[size-2]) return size-1;
			Type_UInt iL=0,iR=size-2,iM;
			//the loop is the binary search.
			while(iL+1<iR){
				iM=(iL+iR)/2;
				if(prob<cdf[iM]) iR=iM;
				else iL=iM;
			}
			return iR;
		} else {
			if(size==1 || prob<=cdf[0]) return (type_res)0;
			if(prob>cdf[size-2]) return size-1;
			Type_UInt iL=0,iR=size-2,iM;
			//the loop is the binary search.
			while(iL+1<iR){
				iM=(iL+iR)/2;
				if(prob<=cdf[iM]) iR=iM;
				else iL=iM;
			}
			return iR;
		}
	}
	
	//other functions:
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_stat
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Set_prob
	(type_res const & number){
		if(number==(type_res)0) return (type_stat)1;
		size=number;
		pdf.resize(size);
		cdf.resize(size);
		type_res i;
		type_real_prob prob=(type_real_prob)1/(type_real_prob)size;
		cdf[0]=pdf[0]=prob;
		for(i=1;i<size;++i){
			pdf[i]=prob;
			cdf[i]=cdf[i-1]+prob;
		}
		return (type_stat)0;
	}
	
	template<typename UIntType,typename ProbRealType,typename RandUIntType>
	typename Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::type_stat
	Distr_Discrete_index<UIntType,ProbRealType,RandUIntType>::Set_prob
	(type_real_prob_seq const & prob){
		if(prob.size()==0) return (type_stat)1;
		type_res i;
		type_real_prob pdf_total=0;
		for(i=0;i<prob.size();++i){
			pdf_total+=(prob[i]>(type_real_prob)0 ? prob[i] : (type_real_prob)0);
		}
		if(pdf_total==(type_real_prob)0) return (type_stat)1;
		size=prob.size();
		pdf.resize(size);
		cdf.resize(size);
		cdf[0]=pdf[0]=(prob[0]>(type_real_prob)0 ? prob[0]/pdf_total : (type_real_prob)0);
		for(i=1;i<size;++i){
			pdf[i]=(prob[i]>(type_real_prob)0 ? prob[i]/pdf_total : (type_real_prob)0);
			cdf[i]=cdf[i-1]+pdf[i];
		}
		return (type_stat)0;
	}
//}(class Distr_Discrete_index)
//(General discrete distribution through index)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Discrete_index_CPP
#endif //#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L