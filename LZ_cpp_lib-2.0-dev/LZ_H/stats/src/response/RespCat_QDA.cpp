#if LZ_DEF_LZ_H_stats_src_response_response_method_category!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_response_RespCat_QDA_CPP
#define LZ_DEF_LZ_H_stats_src_response_RespCat_QDA_CPP 202105L

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//quadratic discrimination analysis{
//class RespCat_QDA{
//public:
	//Constructor:
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::RespCat_QDA
	(void):
	type_base(2,1),size(0),cat_info(this->dim_out,this->dim_in),mode(0){
	}
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::RespCat_QDA
	(Type_UInt const & dim_input):
	type_base(2,(dim_input>0 ? dim_input : 1)),
	size(0),cat_info(this->dim_out,this->dim_in),mode(0){
	}
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::RespCat_QDA
	(type_out const & num_category,Type_UInt const & dim_input):
	type_base((num_category>(type_out)0 ? num_category : (type_out)2),(dim_input>0 ? dim_input : 1)),
	size(0),cat_info(this->dim_out,this->dim_in),mode(0){
	}
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::RespCat_QDA
	(type_this const & model):
	type_base((type_base&)model),size(model.size),cat_info(model.cat_info),mode(model.mode){
	}
	
	//Destructor:
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::~RespCat_QDA
	(void){
	}
	
	//virtual from type_base:
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_out
	RespCat_QDA<OutType,InRealType>::Get_response
	(type_in const & x) const {
		if(!(this->Is_type_in(x))) return (type_out)0;
		type_out res(0),i;
		type_real crit,critm; //the criterion value and its maximum;
		Type_Bool tag=true; //true: no criterion value has been calculated;
		for(i=0;i<this->dim_out;++i){
			if(cat_info.size[i]<=(type_real)0) continue;
			crit=this->Get_crit(i,x);
			if(tag || critm<crit){
				tag=false;
				res=i; //the category res has the max criterion value;
				critm=crit;
			}
		}
		return res;
	}
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_out_weight
	RespCat_QDA<OutType,InRealType>::Get_response_wt
	(type_in const & x) const {
		type_out_weight res=type_out_weight::Zero(this->dim_out,1);
		if(!(this->Is_type_in(x))) return res;
		type_real sum_res=0; //the sum of res;
		type_out i;
		for(i=0;i<this->dim_out;++i){
			if(cat_info.size[i]>(type_real)0){
				res(i,0)=exp(this->Get_crit(i,x)/(type_real)2);
				sum_res+=res(i,0);
			}
		}
		if(sum_res>(type_real)0) res/=sum_res; //res is the posterior probability;
		return res;
	}
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_out
	RespCat_QDA<OutType,InRealType>::Get_response_by_wt
	(type_out_weight const & resp) const {
		if(resp.rows()==0 || resp.cols()==0) return (type_out)0;
		type_out res=0,i,L=std::min((type_out)(resp.rows()),this->dim_out);
		//resp.block(0,0,L,1) takes its maximum at (res,0):
		for(i=1;i<L;++i){
			if(resp(i,0)>resp(res,0)) res=i;
		}
		return res;
	}
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_stat
	RespCat_QDA<OutType,InRealType>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X){
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=Y.size(),i;
		if(X.size()!=n_sam){
			res=(type_stat)11;
			n_sam=std::min(Y.size(),X.size());
		}
		if(n_sam==0) return (type_stat)2; //no sample to train;
		for(i=0;i<n_sam;++i){
			if(Y[i]<(type_out)0 || Y[i]>=this->dim_out ||
				X[i].rows()!=this->dim_in || X[i].cols()==0){
				res=(type_stat)12;
				continue;
			}
			cat_info.recal[Y[i]]=(type_stat)1; //category Y[i] has been changed;
			cat_info.size[Y[i]]+=(type_real)1;
			cat_info.sum1[Y[i]]+=X[i].col(0);
			cat_info.sum2[Y[i]]+=X[i].col(0)*X[i].col(0).transpose();
		}
		res1=cat_info.Cal_para(*this); //re-calculation of cat_info;
		return res1==(type_stat)0 ? res : res1;
	}
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_stat
	RespCat_QDA<OutType,InRealType>::Cal_train
	(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt){
		type_stat res=0,res1; //the return state;
		Type_UInt n_sam=wt.size(),i;
		if(X.size()<n_sam || Y.size()<n_sam){
			res=(type_stat)11;
			n_sam=std::min(Y.size(),X.size());
		}
		if(n_sam==0) return (type_stat)2; //no sample to train;
		for(i=0;i<n_sam;++i){
			if(wt[i]==(type_real)0) continue;
			if(Y[i]<(type_out)0 || Y[i]>=this->dim_out ||
				X[i].rows()!=this->dim_in || X[i].cols()==0){
				res=(type_stat)12;
				continue;
			}
			cat_info.recal[Y[i]]=(type_stat)1; //category Y[i] has been changed;
			cat_info.size[Y[i]]+=wt[i];
			cat_info.sum1[Y[i]]+=wt[i]*(X[i].col(0));
			cat_info.sum2[Y[i]]+=wt[i]*(X[i].col(0)*X[i].col(0).transpose());
		}
		res1=cat_info.Cal_para(*this); //re-calculation of cat_info;
		return res1==(type_stat)0 ? res : res1;
	}
	
	template<typename OutType,typename InRealType>
	inline
	void
	RespCat_QDA<OutType,InRealType>::Set_reset
	(void){
		size=(type_real)0;
		cat_info.Set_reset(*this);
	}
	
	//other functions:
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_stat
	RespCat_QDA<OutType,InRealType>::Set_mode
	(type_stat const mode_new){
		switch(mode_new){
			case (type_stat)0: //QDA;
			case (type_stat)1: //QDA with uniform priori probabilities;
			case (type_stat)2: //by the Mahalanobis distance;
				mode=mode_new;
				return (type_stat)0;
			default:
				return (type_stat)1;
		}
	}
//protected:
	//member functions of type_para:
	
	template<typename OutType,typename InRealType>
	RespCat_QDA<OutType,InRealType>::type_para::type_para
	(type_out const & num_category,Type_UInt const & dim_input):
	pre((type_real)LZ_DEF_const_default_precision),
	Exp(num_category,type_in::Zero(dim_input,1)),
	sum1(num_category,type_in::Zero(dim_input,1)),
	//a very "small" initial variance:
	Var(num_category,
	pre*type_matrix::Identity(dim_input,dim_input)),
	sum2(num_category,type_matrix::Zero(dim_input,dim_input)),
	Var_inv(num_category,
	type_matrix::Identity(dim_input,dim_input)/pre),
	size(num_category,(type_real)0),
	logsize2(num_category,(type_real)0),logdetVar(num_category),
	recal(num_category,(type_stat)0){
	}
	
	template<typename OutType,typename InRealType>
	void
	RespCat_QDA<OutType,InRealType>::type_para::Set_reset
	(type_this const & model){
		type_out L=model.Get_dim_out(),i;
		Type_UInt p=model.Get_dim_in();
		Exp.resize(L);
		sum1.resize(L);
		Var.resize(L);
		sum2.resize(L);
		Var_inv.resize(L);
		size.resize(L);
		logsize2.resize(L);
		logdetVar.resize(L);
		recal.resize(L);
		//a very "small" initial variance:
		type_matrix Var_init=
			pre*type_matrix::Identity(p,p);
		type_matrix Var_inv_init=
			type_matrix::Identity(p,p)/pre;
		for(i=0;i<L;++i){
			Exp[i]=sum1[i]=type_in::Zero(p,1);
			Var[i]=Var_init;
			sum2[i]=type_matrix::Zero(p,p);
			Var_inv[i]=Var_inv_init;
			logsize2[i]=size[i]=(type_real)0;
			recal[i]=(type_stat)0;
		}
	}
	
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_stat
	RespCat_QDA<OutType,InRealType>::type_para::Cal_para
	(type_this const & model){
		type_stat res=0;
		type_out i;
		Type_UInt p=model.Get_dim_in(),j;
		Eigen::SelfAdjointEigenSolver<type_matrix> Var_eigen(p);
		type_real var_max,var_eq0;
		type_matrix Var_eigen_val;
		for(i=0;i<model.Get_dim_out();++i){
			//need re-calculation only when recal[i]==(type_stat)1,
			//and need calculation only when size[i]>(type_real)0:
			if(size[i]<=(type_real)0 || recal[i]==(type_stat)0) continue;
			recal[i]=(type_stat)0; //reset the tag;
			logsize2[i]=log(size[i])*(type_real)2;
			Exp[i]=sum1[i]/size[i];
			//the following are calculating Var[i], Var_inv[i] and logdetVar[i].
			//initial calculation of Var[i]:
			Var[i]=(sum2[i]-size[i]*Exp[i]*Exp[i].transpose())/size[i];
			Var_eigen.compute(Var[i]); //eigenvalue decomposition of Var[i];
			Var_eigen_val=Var_eigen.eigenvalues();
			var_max=Var_eigen_val(p-1,0); //the "scale" of Var[i];
			if(var_max<pre){
				var_max=pre;
				res=(type_stat)21;
			}
			var_eq0=var_max*pre; //the "zero" of Var[i];
			logdetVar[i]=(type_real)1;
			for(j=0;j<p && Var_eigen_val(j,0)<var_eq0;++j){
				res=(type_stat)22;
				Var_eigen_val(j,0)=var_eq0; //make Var[i] invertible;
				logdetVar[i]*=Var_eigen_val(j,0);
			}
			for(;j<p;++j) logdetVar[i]*=Var_eigen_val(j,0);
			Var[i]=Var_eigen.eigenvectors()*
				Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(Var_eigen_val)*
				Var_eigen.eigenvectors().transpose();
			Var_inv[i]=Var_eigen.eigenvectors()*
				Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(Var_eigen_val.cwiseInverse())*
				Var_eigen.eigenvectors().transpose();
		}
		return res;
	}
	
	//other function:
	template<typename OutType,typename InRealType>
	typename RespCat_QDA<OutType,InRealType>::type_real
	RespCat_QDA<OutType,InRealType>::Get_crit
	(type_out const & y,type_in const & x) const {
		type_in x0=x-cat_info.Exp[y]; //the centred x;
		type_real dis=(x0.transpose()*cat_info.Var_inv[y]*x0)(0,0);
		switch(mode){
			case 0:
				return cat_info.logsize2[y]-cat_info.logdetVar[y]-dis;
			case 1:
				return -cat_info.logdetVar[y]-dis;
			case 2:
				return -dis;
			default:
				return (type_real)0;
		}
	}
//}(class RespCat_QDA)
//(quadratic discrimination analysis)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_response_RespCat_QDA_CPP
#endif //#if LZ_DEF_LZ_H_stats_src_response_response_method_category!=202105L