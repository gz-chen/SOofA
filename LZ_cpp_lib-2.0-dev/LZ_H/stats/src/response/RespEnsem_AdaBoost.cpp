#if LZ_DEF_LZ_H_stats_src_response_response_method_ensemble!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_response_RespEnsem_AdaBoost_CPP
#define LZ_DEF_LZ_H_stats_src_response_RespEnsem_AdaBoost_CPP 202105L

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//ensemble method: AdaBoost{
//class RespEnsem_AdaBoost{
//public:
	//Constructor:
	template<typename RespType,Type_Bool LossByInput>
	RespEnsem_AdaBoost<RespType,LossByInput>::RespEnsem_AdaBoost
	(void):
	type_base(),weight_model_base(),loss(NULL){
	}
	
	template<typename RespType,Type_Bool LossByInput>
	RespEnsem_AdaBoost<RespType,LossByInput>::RespEnsem_AdaBoost
	(type_loss & loss_func):
	type_base(),weight_model_base(),loss(std::addressof(loss_func)){
	}
	
	template<typename RespType,Type_Bool LossByInput>
	RespEnsem_AdaBoost<RespType,LossByInput>::RespEnsem_AdaBoost
	(type_model_ptseq const & model_ptseq):
	type_base(model_ptseq),weight_model_base(),loss(NULL){
		if(this->num_model_base>0){
			weight_model_base.assign(
				this->num_model_base,(type_weight)1/(type_weight)(this->num_model_base));
		}
	}
	
	template<typename RespType,Type_Bool LossByInput>
	RespEnsem_AdaBoost<RespType,LossByInput>::RespEnsem_AdaBoost
	(type_model_ptseq const & model_ptseq,type_loss & loss_func):
	type_base(model_ptseq),weight_model_base(),loss(std::addressof(loss_func)){
		if(this->num_model_base>0){
			weight_model_base.assign(
				this->num_model_base,(type_weight)1/(type_weight)(this->num_model_base));
		}
	}
	
	//Destructor:
	
	template<typename RespType,Type_Bool LossByInput>
	RespEnsem_AdaBoost<RespType,LossByInput>::~RespEnsem_AdaBoost
	(void){
	}
	
	//virtual: from type_base:
	
	template<typename RespType,Type_Bool LossByInput>
	typename RespEnsem_AdaBoost<RespType,LossByInput>::type_out_weight
	RespEnsem_AdaBoost<RespType,LossByInput>::Get_response_wt
	(type_in const & x) const {
		if(this->num_model_base==(Type_UInt)0) return type_out_weight(); //error: no model;
		Type_UInt i=0;
		while(i<this->num_model_base && weight_model_base[i]==(type_weight)0) ++i;
		if(i==this->num_model_base) return type_out_weight(); //error: all model weights are 0;
		type_out_weight res=weight_model_base[i]*(this->model_base[i]->Get_response_wt(x));
		for(++i;i<this->num_model_base;++i){
			if(weight_model_base[i]==(type_weight)0) continue;
			res+=(weight_model_base[i]*(this->model_base[i]->Get_response_wt(x)));
		}
		return res;
	}
	
	template<typename RespType,Type_Bool LossByInput>
	typename RespEnsem_AdaBoost<RespType,LossByInput>::type_stat
	RespEnsem_AdaBoost<RespType,LossByInput>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X){
		if(this->num_model_base==0 || loss==NULL){
			return (type_stat)11; //cannot train;
		}
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=Y.size(),i_model,i;
		if(n_sam!=X.size()){
			res=(type_stat)1;
			n_sam=std::min(n_sam,X.size());
		}
		if(n_sam==0) return (type_stat)2; //no need to train;
		type_weight wt_sum=(type_weight)n_sam; //wt[0]+...+wt[n_sam-1];
		//the weight sequence of samples,
		//wt_train is wt*wt_coef, and then normalised to sum to wt_sum:
		//here wt=(1,...,1) and wt_train is proportional to wt_coef:
		//wt_train should be normalised to sum to 1 in order to avoid overflow:
		type_weight_seq wt_train(n_sam,(type_weight)1);
		type_weight_seq val_loss(n_sam); //the values of the loss function at the samples;
		type_weight val_loss_sum;
		for(i_model=0;i_model<this->num_model_base;++i_model){
			//train the i-th base model:
			res1=this->model_base[i_model]->Cal_train(
				type_out_seq(Y.begin(),Y.begin()+n_sam),type_in_seq(X.begin(),X.begin()+n_sam),wt_train);
			if(res1!=(type_stat)0) res=res1;
			//calculate the loss:
			this->infun_Cal_train_cal_loss
				(Y,X,wt_train,wt_sum,n_sam,this->model_base[i_model],*loss,val_loss,val_loss_sum);
			//renew wt_coef and weight_model_base:
			//here wt_train will be normalised so that its sum is 1:
			this->infun_Cal_train_cal_weight
				(n_sam,i_model,val_loss,val_loss_sum,wt_train,weight_model_base);
			//re-normalise wt_train:
			for(i=0;i<n_sam;++i) wt_train[i]*=wt_sum;
		}
		//the following normalises weight_model_base.
		wt_sum=std::accumulate(weight_model_base.begin(),weight_model_base.end(),(type_weight)0);
		if(wt_sum!=(type_weight)0){
			for(i_model=0;i_model<this->num_model_base;++i_model) weight_model_base[i]/=wt_sum;
		}
		return res;
	}
	
	template<typename RespType,Type_Bool LossByInput>
	typename RespEnsem_AdaBoost<RespType,LossByInput>::type_stat
	RespEnsem_AdaBoost<RespType,LossByInput>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt){
		if(this->num_model_base==0 || loss==NULL){
			return (type_stat)11; //cannot train;
		}
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=wt.size(),i_model,i;
		if(X.size()<n_sam || Y.size()<n_sam){
			res=(type_stat)1;
			n_sam=std::min(Y.size(),X.size());
		}
		if(n_sam==0) return (type_stat)2; //no need to train;
		type_weight wt_sum=std::accumulate(wt.begin(),wt.begin()+n_sam,(type_weight)0); //wt[0]+...+wt[n_sam-1];
		if(wt_sum<=(type_weight)0) return (type_stat)3; //no weight to train;
		//the weight sequence of samples,
		//wt_train is wt*wt_coef, and then normalised to sum to wt_sum:
		//wt_train should be normalised to sum to 1 in order to avoid overflow:
		type_weight_seq wt_coef(n_sam,(type_weight)1/(type_weight)n_sam),wt_train(n_sam);
		type_weight wt_train_sum;
		type_weight_seq val_loss(n_sam); //the values of the loss function at the samples;
		type_weight val_loss_sum;
		for(i_model=0;i_model<this->num_model_base;++i_model){
			//begin to calculate the training weights of the samples for the i-th base model.
			wt_train_sum=(type_weight)0;
			for(i=0;i<n_sam;++i) wt_train_sum+=(wt_train[i]=wt[i]*wt_coef[i]);
			wt_train_sum=wt_sum/wt_train_sum;
			for(i=0;i<n_sam;++i) wt_train[i]*=wt_train_sum;
			//train the i-th base model:
			res1=this->model_base[i_model]->Cal_train(
				type_out_seq(Y.begin(),Y.begin()+n_sam),type_in_seq(X.begin(),X.begin()+n_sam),wt_train);
			if(res1!=(type_stat)0) res=res1;
			//calculate the loss:
			this->infun_Cal_train_cal_loss
				(Y,X,wt_train,wt_sum,n_sam,this->model_base[i_model],*loss,val_loss,val_loss_sum);
			//renew wt_coef and weight_model_base:
			this->infun_Cal_train_cal_weight
				(n_sam,i_model,val_loss,val_loss_sum,wt_coef,weight_model_base);
		}
		//the following normalises weight_model_base.
		wt_sum=std::accumulate(weight_model_base.begin(),weight_model_base.end(),(type_weight)0);
		if(wt_sum!=(type_weight)0){
			for(i_model=0;i_model<this->num_model_base;++i_model) weight_model_base[i]/=wt_sum;
		}
		return res;
	}
	
	template<typename RespType,Type_Bool LossByInput>
	void
	RespEnsem_AdaBoost<RespType,LossByInput>::Set_reset
	(void){
		if(this->num_model_base!=weight_model_base.size()){
			weight_model_base.resize(this->num_model_base);
		}
		type_weight wt=(type_weight)1/(type_weight)(this->num_model_base);
		for(Type_UInt i=0;i<this->num_model_base;++i){
			this->model_base[i]->Set_reset();
			weight_model_base[i]=wt;
		}
	}
	
	//other functions:
	
	template<typename RespType,Type_Bool LossByInput>
	inline
	void
	RespEnsem_AdaBoost<RespType,LossByInput>::Set_loss
	(type_loss & loss_func){
		loss=std::addressof(loss_func);
	}
	
//protected:
	
	//other functions:
	
	template<typename RespType,Type_Bool LossByInput>
	void
	RespEnsem_AdaBoost<RespType,LossByInput>::infun_Cal_train_cal_loss
	(type_out_seq const & Y,type_in_seq const & X,
	type_weight_seq const & wt_train,type_weight const & wt_sum,
	Type_UInt const & n_sam,type_model * const & model,type_loss & loss_func,
	type_weight_seq & val_loss,type_weight & val_loss_sum){
		type_weight val_loss_max=(type_weight)0;
		val_loss_sum=(type_weight)0;
		Type_UInt i;
		for(i=0;i<n_sam;++i){
			if(wt_train[i]<=(type_weight)0) continue;
			//calculate the values of the loss function at the samples:
			val_loss[i]=
				incls_RespEnsem_AdaBoost_loss_adaptor<value_loss_by_input>::Cal_loss
				(loss_func,model->Get_response(X[i]),Y[i],X[i]);
			if(val_loss[i]>val_loss_max) val_loss_max=val_loss[i];
		}
		if(val_loss_max==(type_weight)0) return;
		for(i=0;i<n_sam;++i){
			if(wt_train[i]<=(type_weight)0) continue;
			val_loss[i]/=val_loss_max; //normalise val_loss into [0,1];
			val_loss_sum+=(wt_train[i]*val_loss[i]); //val_loss_sum is in (0,wt_sum];
		}
		val_loss_sum/=wt_sum; //val_loss_sum is in (0,1];
		return;
	}
	
	template<typename RespType,Type_Bool LossByInput>
	void
	RespEnsem_AdaBoost<RespType,LossByInput>::infun_Cal_train_cal_weight
	(Type_UInt const & n_sam,Type_UInt const & i_model,
	type_weight_seq const & val_loss,type_weight const & val_loss_sum,
	type_weight_seq & wt_coef,type_weight_seq & seq_weight_model_base){
		type_weight wt_new=this->infun_safeunit(val_loss_sum),wt_coef_sum=(type_weight)0;
		wt_new=sqrt(((type_weight)1-wt_new)/wt_new);
		Type_UInt i;
		for(i=0;i<n_sam;++i){
			wt_coef[i]*=pow(wt_new,(type_weight)2*infun_safeunit(val_loss[i])-(type_weight)1);
			wt_coef_sum+=wt_coef[i];
		}
		for(i=0;i<n_sam;++i) wt_coef[i]/=wt_coef_sum; //normalise wt_coef to avoid overflow;
		seq_weight_model_base[i_model]*=wt_new;
	}
//}(class RespEnsem_AdaBoost)
//(ensemble method: AdaBoost)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_response_RespEnsem_AdaBoost_CPP
#endif //#if LZ_DEF_LZ_H_stats_src_response_response_method_ensemble!=202105L