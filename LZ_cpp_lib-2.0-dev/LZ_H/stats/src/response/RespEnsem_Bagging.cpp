#if LZ_DEF_LZ_H_stats_src_response_response_method_ensemble!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_response_RespEnsem_Bagging_CPP
#define LZ_DEF_LZ_H_stats_src_response_RespEnsem_Bagging_CPP 202105L

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//ensemble method: bagging{
//class RespEnsem_Bagging{
//public:
	//Constructor:
	template<typename RespType,typename RandUIntType>
	RespEnsem_Bagging<RespType,RandUIntType>::RespEnsem_Bagging
	(void):
	type_base(),rand_engine(NULL){
	}
	
	template<typename RespType,typename RandUIntType>
	RespEnsem_Bagging<RespType,RandUIntType>::RespEnsem_Bagging
	(type_rand_engine & engine):
	type_base(),rand_engine(std::addressof(engine)){
	}
	
	template<typename RespType,typename RandUIntType>
	RespEnsem_Bagging<RespType,RandUIntType>::RespEnsem_Bagging
	(type_model_ptseq const & model_ptseq):
	type_base(model_ptseq),rand_engine(NULL){
	}
	
	template<typename RespType,typename RandUIntType>
	RespEnsem_Bagging<RespType,RandUIntType>::RespEnsem_Bagging
	(type_model_ptseq const & model_ptseq,type_rand_engine & engine):
	type_base(model_ptseq),rand_engine(std::addressof(engine)){
	}
	
	//Destructor:
	template<typename RespType,typename RandUIntType>
	RespEnsem_Bagging<RespType,RandUIntType>::~RespEnsem_Bagging
	(void){
	}
	
	//virtual: from type_base:
	
	template<typename RespType,typename RandUIntType>
	typename RespEnsem_Bagging<RespType,RandUIntType>::type_out_weight
	RespEnsem_Bagging<RespType,RandUIntType>::Get_response_wt
	(type_in const & x) const {
		if(this->num_model_base==(Type_UInt)0) return type_out_weight();
		Type_UInt i;
		type_out_weight res=this->model_base[0]->Get_response_wt(x);
		for(i=1;i<this->num_model_base;++i){
			res+=(this->model_base[i]->Get_response_wt(x));
		}
		return res/(type_weight)(this->num_model_base);
	}
	
	template<typename RespType,typename RandUIntType>
	typename RespEnsem_Bagging<RespType,RandUIntType>::type_stat
	RespEnsem_Bagging<RespType,RandUIntType>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X){
		if(this->num_model_base==0 || rand_engine==NULL){
			return (type_stat)11; //cannot train;
		}
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=Y.size(),i_model,i;
		if(n_sam!=X.size()){
			res=(type_stat)1;
			n_sam=std::min(n_sam,X.size());
		}
		if(n_sam==0) return (type_stat)2; //no need to train;
		//the random chooser for bootstrap:
		std::uniform_int_distribution<Type_UInt> distr_uni(0,n_sam-1);
		type_weight_seq wt_b(n_sam); //the weight sequence for bagging;
		for(i_model=0;i_model<this->num_model_base;++i_model){
			wt_b.assign(n_sam,(type_weight)0); //initialise wt_b;
			for(i=0;i<n_sam;++i) wt_b[distr_uni(*rand_engine)]+=(type_weight)1; //bootstrap;
			//train the i-th base model:
			res1=this->model_base[i_model]->Cal_train(
				type_out_seq(Y.begin(),Y.begin()+n_sam),type_in_seq(X.begin(),X.begin()+n_sam),wt_b);
			if(res1!=(type_stat)0) res=res1;
		}
		return res;
	}
	
	template<typename RespType,typename RandUIntType>
	typename RespEnsem_Bagging<RespType,RandUIntType>::type_stat
	RespEnsem_Bagging<RespType,RandUIntType>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt){
		if(this->num_model_base==0 || rand_engine==NULL){
			return (type_stat)11; //cannot train;
		}
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=wt.size(),i_model,i,j;
		if(X.size()<n_sam || Y.size()<n_sam){
			res=(type_stat)1;
			n_sam=std::min(Y.size(),X.size());
		}
		if(n_sam==0) return (type_stat)2; //no need to train;
		//the random chooser for bootstrap:
		std::uniform_int_distribution<Type_UInt> distr_uni(0,n_sam-1);
		type_weight_seq wt_b(n_sam); //the weight sequence for bagging;
		for(i_model=0;i_model<this->num_model_base;++i_model){
			wt_b.assign(n_sam,(type_weight)0); //initialise wt_b;
			for(i=0;i<n_sam;++i){
				j=distr_uni(*rand_engine); //bootstrap;
				wt_b[j]+=wt[j];
			}
			//train the i-th base model:
			res1=this->model_base[i_model]->Cal_train(
				type_out_seq(Y.begin(),Y.begin()+n_sam),type_in_seq(X.begin(),X.begin()+n_sam),wt_b);
			if(res1!=(type_stat)0) res=res1;
		}
		return res;
	}
	
	template<typename RespType,typename RandUIntType>
	void
	RespEnsem_Bagging<RespType,RandUIntType>::Set_reset
	(void){
		for(Type_UInt i=0;i<this->num_model_base;++i) this->model_base[i]->Set_reset();
	}
	
	//other functions:
	
	template<typename RespType,typename RandUIntType>
	inline
	void
	RespEnsem_Bagging<RespType,RandUIntType>::Set_rand_engine
	(type_rand_engine & engine){
		rand_engine=std::addressof(engine);
	}
//}(class RespEnsem_Bagging)
//(ensemble method: bagging)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_response_RespEnsem_Bagging_CPP
#endif //#if LZ_DEF_LZ_H_stats_src_response_response_method_ensemble!=202105L