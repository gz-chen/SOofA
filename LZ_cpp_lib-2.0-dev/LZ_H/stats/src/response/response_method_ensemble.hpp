#if LZ_DEF_LZ_H_stats_response!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_response_response_method_ensemble
#define LZ_DEF_LZ_H_stats_src_response_response_method_ensemble 202105L

#include"LZ_H/math/distribution.hpp"

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//Bagging (bootstrap aggregation):
//ensemble method: bagging{
template<typename RespType=
	Response<Type_Real,Eigen::MatrixX<Type_Real>,Type_Real,Type_Real>,
	typename RandUIntType=typename math::Type_RandEngine::result_type>
class RespEnsem_Bagging: public Response_Ensemble<RespType>{
	LZ_DEF_func_check_traits(std::is_unsigned_integral<RandUIntType>::value);
	public:
		//Type: this class:
		typedef RespEnsem_Bagging<RespType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef Response_Ensemble<RespType> type_base;
		//Type: input:
		typedef typename type_base::type_in type_in;
		//Type: output:
		typedef typename type_base::type_out type_out;
		//Type: weight-able output:
		typedef typename type_base::type_out_weight type_out_weight;
		//Type: the real number, or the weights:
		typedef typename type_base::type_weight type_weight;
		//Type: the type of state:
		typedef typename type_base::type_stat type_stat;
		//Type: the sequence of inputs:
		typedef typename type_base::type_in_seq type_in_seq;
		//Type: the sequence of outputs:
		typedef typename type_base::type_out_seq type_out_seq;
		//Type: the sequence of weights:
		typedef typename type_base::type_weight_seq type_weight_seq;
		//Type: the base model:
		typedef typename type_base::type_model type_model;
		//Type: the sequence of pointers to base models:
		typedef typename type_base::type_model_ptseq type_model_ptseq;
		//Type: the random engine:
		typedef math::RandomEngine<RandUIntType> type_rand_engine;
		
		//Constructor:
		RespEnsem_Bagging(void);
		RespEnsem_Bagging(type_rand_engine & engine);
		RespEnsem_Bagging(type_model_ptseq const & model_ptseq);
		RespEnsem_Bagging(type_model_ptseq const & model_ptseq,type_rand_engine & engine);
		
		//Destructor:
		virtual ~RespEnsem_Bagging(void);
		
		//virtual: from type_base:
		virtual type_out_weight Get_response_wt(type_in const & x) const;
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X);
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt);
		virtual void Set_reset(void);
		
		//Set: the random engine:
		void Set_rand_engine(type_rand_engine & engine);
	protected:
		type_rand_engine * rand_engine; //the random engine for bootstrap;
}; //class RespEnsem_Bagging;
//}(ensemble method: bagging)
/****************************************************************************************************/
//AdaBoost (adaptive boosting):
//ensemble method: AdaBoost{
//template parameters{
//RespType: the type of base models.
//LossByInput: true: the loss function depends on the input, i.e. in the form L(decision,true_val,input);
//false(default): the loss function is of the form L(decision,true_value).
//}
template<typename RespType=
	Response<Type_Real,Eigen::MatrixX<Type_Real>,Type_Real,Type_Real>,
	Type_Bool LossByInput=(Type_Bool)false>
class RespEnsem_AdaBoost: public Response_Ensemble<RespType>{
	public:
		//template parameter value:
		static Type_Bool constexpr value_loss_by_input=LossByInput;
		
		//Type: this class:
		typedef RespEnsem_AdaBoost<RespType,LossByInput> type_this;
		//Type: type of the base class:
		typedef Response_Ensemble<RespType> type_base;
		//Type: input:
		typedef typename type_base::type_in type_in;
		//Type: output:
		typedef typename type_base::type_out type_out;
		//Type: weight-able output:
		typedef typename type_base::type_out_weight type_out_weight;
		//Type: the real number, or the weights:
		typedef typename type_base::type_weight type_weight;
		//Type: the type of state:
		typedef typename type_base::type_stat type_stat;
		//Type: the sequence of inputs:
		typedef typename type_base::type_in_seq type_in_seq;
		//Type: the sequence of outputs:
		typedef typename type_base::type_out_seq type_out_seq;
		//Type: the sequence of weights:
		typedef typename type_base::type_weight_seq type_weight_seq;
		//Type: the base model:
		typedef typename type_base::type_model type_model;
		//Type: the sequence of pointers to base models:
		typedef typename type_base::type_model_ptseq type_model_ptseq;
		//Type: the class of the loss function:
		typedef
			typename std::conditional<value_loss_by_input,
				std::function<type_weight(type_out const &,type_out const &,type_in const &)>,
				std::function<type_weight(type_out const &,type_out const &)>
			>::type
			type_loss;
		
		//Constructor:
		RespEnsem_AdaBoost(void);
		RespEnsem_AdaBoost(type_loss & loss_func);
		RespEnsem_AdaBoost(type_model_ptseq const & model_ptseq);
		RespEnsem_AdaBoost(type_model_ptseq const & model_ptseq,type_loss & loss_func);
		
		//Destructor:
		virtual ~RespEnsem_AdaBoost(void);
		
		//virtual: from type_base:
		virtual type_out_weight Get_response_wt(type_in const & x) const;
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X);
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt);
		virtual void Set_reset(void);
		
		//Set: the loss function:
		void Set_loss(type_loss & loss_func);
	protected:
		//members:
		type_weight_seq weight_model_base; //weights of the base models;
		type_loss * loss; //the loss function;
		
		//inner function: calculate the weighted loss of the base models during the training:
		void infun_Cal_train_cal_loss
		(type_out_seq const & Y,type_in_seq const & X,
			type_weight_seq const & wt_train,type_weight const & wt_sum,
			Type_UInt const & n_sam,type_model * const & model,type_loss & loss_func,
		 	type_weight_seq & val_loss,type_weight & val_loss_sum);
		//inner function: renew the weights of the samples and the base models during the training:
		void infun_Cal_train_cal_weight(Type_UInt const & n_sam,Type_UInt const & i_model,
			type_weight_seq const & val_loss,type_weight const & val_loss_sum,
			type_weight_seq & wt_coef,type_weight_seq & seq_weight_model_base);
		//inner function: affine transformation from [0,1] to [pre,1-pre]:
		type_weight constexpr infun_safeunit
		(type_weight const & x,type_weight const pre=(type_weight)(1e-5)) const {
			return x*((type_weight)1-pre)+pre/(type_weight)2;
		}
}; //class RespEnsem_AdaBoost;

//class incls_RespEnsem_AdaBoost_loss_adaptor:
//an assistant class to calculate the two different types of loss functions:
//template parameter:
//LossByInput: false: loss(decision,val_true); true: loss(decision,val_true,input);
//the general template is for LossByInput==false:
template<Type_Bool LossByInput>
class incls_RespEnsem_AdaBoost_loss_adaptor{
	public:
		template<typename LossType,typename OutType,typename InType>
		static typename LossType::result_type Cal_loss
		(LossType & loss_func,
		OutType const & decision,OutType const & val_true,InType const & input){
			return loss_func(decision,val_true);
		}
}; //class incls_RespEnsem_AdaBoost_loss_adaptor;
//class incls_RespEnsem_AdaBoost_loss_adaptor:
//specialisation for LossByInput==true:
template<>
class incls_RespEnsem_AdaBoost_loss_adaptor<(Type_Bool)true>{
	public:
		template<typename LossType,typename OutType,typename InType>
		static typename LossType::result_type Cal_loss
		(LossType & loss_func,
		OutType const & decision,OutType const & val_true,InType const & input){
			return loss_func(decision,val_true,input);
		}
}; //class incls_RespEnsem_AdaBoost_loss_adaptor;

//}(ensemble method: AdaBoost)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/response/RespEnsem_Bagging.cpp"
#include"LZ_H/stats/src/response/RespEnsem_AdaBoost.cpp"

#endif //#ifndef LZ_DEF_LZ_H_stats_src_response_response_method_ensemble
#endif //#if LZ_DEF_LZ_H_stats_response!=202105L