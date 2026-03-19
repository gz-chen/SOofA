#ifndef LZ_DEF_extLIB_Eigen
	#error "It requires the library: Eigen."
#else

#ifndef LZ_DEF_LZ_H_stats_response
#define LZ_DEF_LZ_H_stats_response 202105L

#include"LZ_H/stats/basic.hpp"

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//base class for input-response models{
template<typename OutType=Type_Real,typename InType=Eigen::MatrixX<Type_Real>,
	typename WeightOutType=OutType,typename WeightType=Type_Real>
class Response{
	LZ_DEF_func_check_traits(std::is_floating_point<WeightType>::value);
	public:
		//Type: this class:
		typedef Response<OutType,InType,WeightOutType,WeightType> type_this;
		//Type: type of the base class:
		typedef void type_base;
		//Type: input:
		typedef InType type_in;
		//Type: output:
		typedef OutType type_out;
		//Type: weight-able output:
		typedef WeightOutType type_out_weight;
		//Type: the real number, or the weights:
		typedef WeightType type_weight;
		//Type: the type of state:
		typedef Type_Stat type_stat;
		//Type: the sequence of inputs:
		typedef std::vector<type_in> type_in_seq;
		//Type: the sequence of outputs:
		typedef std::vector<type_out> type_out_seq;
		//Type: the sequence of weights:
		typedef std::vector<type_weight> type_weight_seq;
		
		//Constructor:
		Response()=default;
		//Destructor:
		virtual ~Response()=default;
		
		//Get: the predicted value:
		virtual type_out Get_response(type_in const & x) const =0;
		virtual type_out operator ()(type_in const & x) const final {
			return this->Get_response(x);
		}
		//Get: the weight-able prediction:
		virtual type_out_weight Get_response_wt(type_in const & x) const =0;
		//Get: get the prediction from the weight-able prediction:
		virtual type_out Get_response_by_wt(type_out_weight const & resp) const =0;
		//Cal: train the model with data (and equal weights):
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X)=0;
		//Cal: train the model with data and weights:
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt)=0;
		//Set: reset:
		virtual void Set_reset(void)=0;
}; //class Response;
//}(base class for input-response models)
/****************************************************************************************************/
//base class for real-vector input and category output models{
template<typename OutType=Type_UInt,typename InRealType=Type_Real>
class Response_Category_RealVec:
public Response<OutType,Eigen::MatrixX<InRealType>,Eigen::MatrixX<InRealType>,InRealType>{
	LZ_DEF_func_check_traits(std::is_integral<OutType>::value);
	LZ_DEF_func_check_traits(std::is_floating_point<InRealType>::value);
	public:
		//Type: this class:
		typedef Response_Category_RealVec<OutType,InRealType> type_this;
		//Type: type of the base class:
		typedef Response<OutType,Eigen::MatrixX<InRealType>,Eigen::MatrixX<InRealType>,InRealType>
			type_base;
		//Type: input:
		typedef typename type_base::type_in type_in;
		//Type: output:
		typedef typename type_base::type_out type_out;
		//Type: weight-able output:
		typedef typename type_base::type_out_weight type_out_weight;
		//Type: the weights:
		typedef typename type_base::type_weight type_weight;
		//Type: the real number:
		typedef InRealType type_real; //the same as type_weight;
		//Type: the type of state:
		typedef typename type_base::type_stat type_stat;
		//Type: the sequence of inputs:
		typedef typename type_base::type_in_seq type_in_seq;
		//Type: the sequence of outputs:
		typedef typename type_base::type_out_seq type_out_seq;
		//Type: the sequence of weights:
		typedef typename type_base::type_weight_seq type_weight_seq;
		
		//Constructor:
		Response_Category_RealVec(void):
		type_base(),dim_out(2),dim_in(1){
		}
		Response_Category_RealVec(Type_UInt const & dim_input):
		type_base(),dim_out(2),dim_in(dim_input>0 ? dim_input : 1){
		}
		Response_Category_RealVec(type_out const & num_category,Type_UInt const & dim_input):
		type_base(),dim_out(num_category>(type_out)0 ? num_category : (type_out)2),
		dim_in(dim_input>0 ? dim_input : 1){
		}
		Response_Category_RealVec(type_this const & model):
		type_base((type_base&)model),dim_out(model.dim_out),dim_in(model.dim_in){
		}
		
		//Destructor:
		~Response_Category_RealVec(void){}
		
		//virtual: from type_base:
		virtual type_out Get_response(type_in const & x) const =0;
		virtual type_out_weight Get_response_wt(type_in const & x) const =0;
		virtual type_out Get_response_by_wt(type_out_weight const & resp) const =0;
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X)=0;
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt)=0;
		virtual void Set_reset(void)=0;
		
		//Is: a matrix of size dim_in*1:
		Type_Bool Is_type_in(type_in const & x) const {
			return x.rows()==this->dim_in && x.cols()==(Type_UInt)1;
		}
		
		//Get: the number of categories:
		type_out Get_dim_out() const {
			return dim_out;
		}
		//Get: the dimension of the input:
		Type_UInt Get_dim_in() const {
			return dim_in;
		}
		
		//Set: dim_out and dim_in:
		type_stat Set_dim(type_out const & num_category,Type_UInt const & dim_input){
			if(num_category<=(type_out)0 || dim_input<=0) return (type_stat)1;
			dim_out=num_category;
			dim_in=dim_input;
			this->Set_reset();
			return (type_stat)0;
		}
	protected:
		type_out dim_out; //the number of candidate values of the response;
		Type_UInt dim_in; //the dimension of the input;
}; //class Response_Category_RealVec;
//}(base class for real-vector input and category output models)
/****************************************************************************************************/
//base class for real-vector input and 1-dim real output models{
template<typename RealType=Type_Real>
class Response_Real_RealVec:
public Response<RealType,Eigen::MatrixX<RealType>,RealType,RealType>{
	LZ_DEF_func_check_traits(std::is_floating_point<RealType>::value);
	public:
		//Type: this class:
		typedef Response_Real_RealVec<RealType> type_this;
		//Type: type of the base class:
		typedef Response<RealType,Eigen::MatrixX<RealType>,RealType,RealType> type_base;
		//Type: input:
		typedef typename type_base::type_in type_in;
		//Type: output:
		typedef typename type_base::type_out type_out;
		//Type: weight-able output:
		typedef typename type_base::type_out_weight type_out_weight; //the same as type_out;
		//Type: the weights:
		typedef typename type_base::type_weight type_weight;
		//Type: the real number:
		typedef RealType type_real; //the same as type_weight, type_out;
		//Type: the type of state:
		typedef typename type_base::type_stat type_stat;
		//Type: the sequence of inputs:
		typedef typename type_base::type_in_seq type_in_seq;
		//Type: the sequence of outputs:
		typedef typename type_base::type_out_seq type_out_seq;
		//Type: the sequence of weights:
		typedef typename type_base::type_weight_seq type_weight_seq;
		
		//Constructor:
		Response_Real_RealVec(void):
		type_base(),dim_out(1),dim_in(1){
		}
		Response_Real_RealVec(Type_UInt const & dim_input):
		type_base(),dim_out(1),dim_in(dim_input>0 ? dim_input : 1){
		}
		Response_Real_RealVec(type_this const & model):
		type_base((type_base&)model),dim_out(1),dim_in(model.dim_in){
		}
		
		//Destructor:
		~Response_Real_RealVec(void){}
		
		//virtual: from type_base:
		virtual type_out Get_response(type_in const & x) const =0;
		virtual type_out_weight Get_response_wt(type_in const & x) const final {
			return this->operator ()(x);
		}
		virtual type_out Get_response_by_wt(type_out_weight const & resp) const final {
			return resp;
		}
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X)=0;
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt)=0;
		virtual void Set_reset(void)=0;
		
		//Is: a matrix of size dim_in*1:
		Type_Bool Is_type_in(type_in const & x) const {
			return x.rows()==this->dim_in && x.cols()==(Type_UInt)1;
		}
		
		//Get: the dimension of the output:
		Type_UInt Get_dim_out() const {
			return dim_out;
		}
		//Get: the dimension of the input:
		Type_UInt Get_dim_in() const {
			return dim_in;
		}
		
		//Set: dim_in:
		type_stat Set_dim(Type_UInt const & dim_input){
			if(dim_input<=0) return (type_stat)1;
			dim_in=dim_input;
			this->Set_reset();
			return (type_stat)0;
		}
	protected:
		Type_UInt const dim_out; //the dimensions of output;
		Type_UInt dim_in; //the dimensions of input;
}; //class Response_Real_RealVec;
//}(base class for real-vector input and category output models)
/****************************************************************************************************/
//base class for ensemble models{
template<typename RespType=
	Response<Type_Real,Eigen::MatrixX<Type_Real>,Type_Real,Type_Real> >
class Response_Ensemble:
public Response<typename RespType::type_out,typename RespType::type_in,
	typename RespType::type_out_weight,typename RespType::type_weight>{
	LZ_DEF_func_check_traits((std::is_base_of<
		Response<typename RespType::type_out,typename RespType::type_in,
		typename RespType::type_out_weight,typename RespType::type_weight>,RespType>::value));
	public:
		//Type: this class:
		typedef Response_Ensemble<RespType> type_this;
		//Type: type of the base class:
		typedef
			Response<typename RespType::type_out,typename RespType::type_in,
				typename RespType::type_out_weight,typename RespType::type_weight>
			type_base;
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
		typedef RespType type_model;
		//Type: the sequence of pointers to base models:
		typedef std::vector<type_model*> type_model_ptseq;
		
		//Constructor:
		Response_Ensemble(void):
		num_model_base(0),model_base(),model_out(NULL){
		}
		Response_Ensemble(type_model_ptseq const & model_ptseq):
		num_model_base(0),model_base(),model_out(NULL){
			this->Set_model(model_ptseq);
		}
		Response_Ensemble(type_model_ptseq const & model_ptseq,type_model const & model_output):
		num_model_base(0),model_base(),
		model_out(std::addressof(model_output)==this ? NULL : std::addressof(model_output)){
			this->Set_model(model_ptseq);
		}
		//Destructor:
		virtual ~Response_Ensemble(void){}
		
		//virtual: from type_base:
		virtual type_out Get_response(type_in const & x) const final {
			return this->Get_response_by_wt(this->Get_response_wt(x));
		}
		virtual type_out_weight Get_response_wt(type_in const & x) const =0;
		virtual type_out Get_response_by_wt(type_out_weight const & resp) const final {
			return (const_cast<type_this*>(this)->Get_model_output())==NULL ? type_out() :
				const_cast<type_this*>(this)->Get_model_output()->Get_response_by_wt(resp);
		}
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X)=0;
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt)=0;
		virtual void Set_reset(void)=0;
		
		//Get: the number of base models:
		Type_UInt Get_num_model() const {
			return num_model_base;
		}
		//Get: the sequence of the pointers to base models:
		type_model_ptseq Get_model(){
			return model_base;
		}
		//Get: the reference to a base model:
		type_model* Get_model(Type_UInt const & id_model){
			return id_model<num_model_base ? (model_base[id_model]) :
				(num_model_base>0 ? (model_base[0]) : NULL);
		}
		//Get: the model for generating output from the weight-able output:
		type_model* Get_model_output(){
			return model_out==NULL ? this->Get_model((Type_UInt)0) : model_out;
		}
		//Set: the sequence of the pointers to base models:
		type_stat Set_model(type_model_ptseq const & model_ptseq){
			if(model_ptseq.size()==0) return (type_stat)1;
			Type_UInt i;
			for(i=0;i<model_ptseq.size();++i){
				if(model_ptseq[i]==NULL) return (type_stat)2;
			}
			model_base=model_ptseq;
			num_model_base=model_base.size();
			return (type_stat)0;
		}
		//Set: set the model for generating output from the weight-able output:
		void Set_model_output(type_model * model_output=NULL){
			if((type_base)model_output!=(type_base)this) model_out=model_output;
		}
		void Set_model_output(type_model const & model_output){
			if((type_base)(std::addressof(model_output))!=(type_base)this){
				model_out=std::addressof(model_output);
			}
		}
	protected:
		type_model_ptseq model_base; //a sequence of base models;
		//use this model to generate the output from the weight-able output:
		//NULL(default): use model_base[0]; otherwise: pointer to a model:
		type_model* model_out;
		Type_UInt num_model_base; //number of base models;
}; //class Response_Ensemble;
//}(base class for ensemble models)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

#include"LZ_H/stats/src/response/response_method_category.hpp"
#include"LZ_H/stats/src/response/response_method_ensemble.hpp"

#endif //#ifndef LZ_DEF_LZ_H_stats_response
#endif //#ifndef LZ_DEF_extLIB_Eigen