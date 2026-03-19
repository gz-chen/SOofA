#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Normal_multi_CPP
#define LZ_DEF_LZ_H_math_src_distribution_Distr_Normal_multi_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//Multi-variate normal distribution:
//class Distr_Normal_multi{
//public:
	//Constructor:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Distr_Normal_multi
	():
	type_base((Type_UInt)1,type_this::stat_category_continuous),
	std_normal((type_real)0,(type_real)1){
		para.Exp=type_res::Zero(1,1);
		para.Var=type_matrix::Identity(1,1);
		para_in.Var_inv=para_in.Var_tran=type_matrix::Identity(1,1);
		para_in.coef_PDF=(type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi;
		para_in.pre=(type_real)LZ_DEF_const_default_precision;
		para_in.rank=(Type_UInt)1;
		this->stat_able=(unsigned char)15;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Distr_Normal_multi
	(Type_UInt const & dimension):
	type_base(dimension,type_this::stat_category_continuous),
	std_normal((type_real)0,(type_real)1){
		para.Exp=type_res::Zero(this->dim,1);
		para.Var=type_matrix::Identity(this->dim,this->dim);
		para_in.Var_inv=para_in.Var_tran=type_matrix::Identity(this->dim,this->dim);
		para_in.coef_PDF=
			pow((type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi,(type_real_prob)(this->dim));
		para_in.pre=(type_real)LZ_DEF_const_default_precision;
		para_in.rank=this->dim;
		this->stat_able=(unsigned char)15;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Distr_Normal_multi
	(type_res const & mean,type_real var):
	type_base(mean.rows(),type_this::stat_category_continuous),
	std_normal((type_real)0,(type_real)1){
		para_in.pre=(type_real)LZ_DEF_const_default_precision;
		para.Exp= (mean.cols()==0 || mean.rows()==0) ?
			(type_res const &)(type_res::Zero(this->dim,1)) : (type_res const &)(mean.col(0));
		//if(var<=0){}else{}:
		if(var<=para_in.pre){
			this->stat_category=type_this::stat_category_discrete;
			para.Var=type_matrix::Zero(this->dim,this->dim);
			para_in.Var_tran=type_matrix();
			para_in.rank=(Type_UInt)0;
			this->stat_able=(unsigned char)9;
		} else {
			para.Var=type_matrix::Identity(this->dim,this->dim)*var;
			type_real var_sqrt=sqrt(var);
			para_in.Var_tran=type_matrix::Identity(this->dim,this->dim)*var_sqrt;
			para_in.Var_inv=type_matrix::Identity(this->dim,this->dim)/var;
			para_in.coef_PDF=
				pow((type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi/var_sqrt,(type_real_prob)(this->dim));
			para_in.rank=this->dim;
			this->stat_able=(unsigned char)15;
		}
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Distr_Normal_multi
	(type_res const & mean,type_matrix const & var):
	type_base(mean.rows(),type_this::stat_category_continuous),
	std_normal((type_real)0,(type_real)1){
		para_in.pre=(type_real)LZ_DEF_const_default_precision;
		para.Exp= (mean.cols()==0 || mean.rows()==0) ?
			(type_res const &)(type_res::Zero(this->dim,1)) : (type_res const &)(mean.col(0));
		type_stat tag_var; //var: 0: zero; 1: positive; 2: semi-positive and not positive; 10: illegal.
		if(var.rows()!=this->dim || var.cols()!=this->dim ||
			!var.isApprox(var.transpose().eval(),para_in.pre)){
			tag_var=(type_stat)10; //illegal var;
		} else {
			//eigen decomposition of var:
			Eigen::SelfAdjointEigenSolver<type_matrix> var_eigen(var);
			type_matrix var_eigen_val=var_eigen.eigenvalues();
			//calculate the max component-wise var:
			type_real var_max=std::max(abs(var_eigen_val(0,0)),abs(var_eigen_val(this->dim-1,0)));
			//var is zero ?:
			if(var_max<=para_in.pre){
				tag_var=(type_stat)0; //var is zero;
				this->stat_category=type_this::stat_category_discrete;
				para.Var=type_matrix::Zero(this->dim,this->dim);
				para_in.Var_tran=type_matrix();
				para_in.rank=(Type_UInt)0;
				this->stat_able=(unsigned char)9;
			} else {
				type_real var_eq0=para_in.pre*var_max; //approx 0 of var;
				//var is not semi-positive ?:
				if(var_eigen_val(0,0)<-var_eq0){
					tag_var=(type_stat)10; //illegal var;
				} else {
					///type_matrix var_eigen_D=Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val);
					Type_UInt i=0;
					para_in.rank=this->dim;
					while(i<this->dim && var_eigen_val(i,0)<=var_eq0){
						var_eigen_val(i,0)=(type_real)0;
						para_in.rank--;
						++i;
					}
					//var is not of full rank ?:
					if((para_in.rank)<this->dim){
						tag_var=(type_stat)2;
						this->stat_category=type_this::stat_category_singular;
						para.Var=Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val);
						para.Var=var_eigen.eigenvectors()*para.Var*var_eigen.eigenvectors().transpose();
						para_in.Var_tran.resize(this->dim,para_in.rank);
						Type_UInt i0=this->dim-para_in.rank;
						for(i=0;i<para_in.rank;++i){
							para_in.Var_tran.col(i)=
								var_eigen.eigenvectors().col(i0+i)*sqrt(var_eigen_val(i0+i,0));
						}
						this->stat_able=(unsigned char)9;
					} else {
						tag_var=(type_stat)1;
						///this->stat_category=type_this::stat_category_continuous;
						para.Var=var;
						para_in.Var_tran=var_eigen.eigenvectors()*
							Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val.cwiseSqrt());
						para_in.Var_inv=var_eigen.eigenvectors()*
							Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val.cwiseInverse())*
							var_eigen.eigenvectors().transpose();
						para_in.coef_PDF=
							pow((type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi,
							(type_real_prob)(this->dim))/
							sqrt(var_eigen_val.prod());
						this->stat_able=(unsigned char)15;
					}
				}
			}
		}
		if(tag_var==(type_stat)10){
			//this->stat_category=type_this::stat_category_continuous;
			para_in.Var_tran=para_in.Var_inv=para.Var=
				type_matrix::Identity(this->dim,this->dim);
			para_in.coef_PDF=
				pow((type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi,(type_real_prob)(this->dim));
			para_in.rank=this->dim;
			this->stat_able=(unsigned char)15;
		}
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Distr_Normal_multi
	(type_this const & normal):
	type_base((type_base const &)normal),para(normal.para),para_in(normal.para_in),
	std_normal((type_real)0,(type_real)1){
	}
	
	//Destructor:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::~Distr_Normal_multi(){
	}
	
	//operator:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_this&
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::operator =
	(type_this const & distr){
		if(&distr==this) return *this;
		((type_base*)this)->operator =((type_base const &)distr);
		para=distr.para;
		para_in=distr.para_in;
		return *this;
	}
	
	//virtual from type_base:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_res
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::rand
	(type_rand_engine & engine) const {
		if(!this->Get_able_rand()) return type_res();
		if(para_in.rank==0) return para.Exp; //constant;
		type_res sam(para_in.rank,1); //the standard normal sample;
		Type_UInt i;
		for(i=0;i<para_in.rank;++i) sam(i,0)=std_normal(engine); //generating the standard normal variable;
		return para_in.Var_tran*sam+para.Exp;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_res_seq
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::sample
	(Type_UInt const num,type_rand_engine & engine) const {
		type_res_seq res; //the result;
		if(!this->Get_able_rand() || num==0) return res;
		Type_UInt i,i_num;
		res.resize(num);
		//constant:
		if(para_in.rank==0){
			for(i_num=0;i_num<num;++i_num) res[i_num]=para.Exp;
			return res;
		}
		//following is the case rank>0.
		type_res sam(para_in.rank,1); //the standard normal sample;
		for(i_num=0;i_num<num;++i_num){
			for(i=0;i<para_in.rank;++i) sam(i,0)=std_normal(engine); //generating sam;
			res[i_num]=para_in.Var_tran*sam+para.Exp;
		}
		return res;
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	inline
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::PDF
	(type_res const & x) const {
		if(!this->Get_able_PDF()) return (type_real_prob)0;
		type_res x0=x-para.Exp;
		return para_in.coef_PDF*type_real_prob(exp((x0.transpose()*para_in.Var_inv*x0)(0,0)/type_real(-2)));
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_real_prob
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::CDF
	(Type_MatTemp<type_real> const & x) const {
		if(!this->Get_able_CDF()) return (type_real_prob)0;
		Type_UInt i;
		type_real coef_0=sqrt(-log(para_in.pre)*(type_real)2); ////the univariate N(0,1)'s PDF is quite small at coef_0;
		type_real_prob res=(type_real_prob)0; //the result;
		type_res x_min(this->dim,1),x0=x-para.Exp; //x_min: the min x being considered; x0: centered x;
		for(i=0;i<this->dim;++i){
			x_min(i,0)=-sqrt(para.Var(i,i))*coef_0;
			if(x_min(i,0)>=x0(i,0)) return res;
		} //cwise-checking x_min<x0;
		//The following is the numerical integration, using the composed 2-knot Gauss-Legendre method.
		Type_ArrayTemp<type_real> std_knot(2); //the roots of the Legendre polynomial in [0,1];
		std_knot[0]=(type_real)1/((type_real)2*sqrt((type_real)3));
		std_knot[1]=(type_real)0.5+std_knot[0];
		std_knot[0]=(type_real)1-std_knot[1];
		//information about the knots:
		//cwise-min(x0,-x_min), regular segment width, 2 segment points, 2 knots, real segment width:
		type_matrix knots(this->dim,7);
		//temp values in calculating the integration:
		//the measure of one segment of the integration region,
		//the partial sum in one segment;
		Type_ArrayTemp<type_real> inte_temp(2);
		inte_temp[0]=(type_real)1;
		//initialising:
		for(i=0;i<this->dim;++i){
			knots(i,0)= x0(i,0)<=-(x_min(i,0)) ? x0(i,0) : -(x_min(i,0));
			knots(i,1)=pow(para_in.pre,(type_real)1/(type_real)3)*(-x_min(i,0));
			knots(i,2)=x_min(i,0);
			knots(i,3)=knots(i,2)+knots(i,1);
			knots(i,4)=knots(i,2)*std_knot[1]+knots(i,3)*std_knot[0];
			knots(i,5)=knots(i,2)*std_knot[0]+knots(i,3)*std_knot[1];
			knots(i,6)=knots(i,1);
			inte_temp[0]*=knots(i,6);
		}
		//iteration:
		Type_UInt i_seg;
		Type_ArrayTemp<Type_UInt> i_knot(this->dim);
		type_res x_inte(this->dim,1);
		while(true){
			//calculating the integration in one segment:
			inte_temp[1]=(type_real)0;
			std::fill_n(std::begin(i_knot),this->dim,(Type_UInt)0);
			while(true){
				for(i=0;i<this->dim;++i) x_inte(i,0)=knots(i,4+i_knot[i]);
				inte_temp[1]+=exp((x_inte.transpose()*para_in.Var_inv*x_inte)(0,0)/(type_real)(-2));
				i=this->dim-1;
				while(i_knot[i]==1){
					if(i==0) goto GTS_iter_end_knot;
					--i;
				}
				i_knot[i]=1;
				while(++i<this->dim) i_knot[i]=0;
			}
			//exit of the iteration for knots:
			GTS_iter_end_knot:
			res+=(type_real_prob)(inte_temp[1]*inte_temp[0]);
			//move to the next segment:
			i_seg=this->dim-1;
			while(knots(i_seg,3)==x0(i_seg,0)){
				if(i_seg==0) goto GTS_iter_end_seg;
				--i_seg;
			}
			knots(i_seg,2)+=knots(i_seg,1);
			knots(i_seg,3)+=knots(i_seg,1);
			if(knots(i_seg,3)>=knots(i_seg,0)){
				knots(i_seg,3)=x0(i_seg,0);
				knots(i_seg,6)=knots(i_seg,3)-knots(i_seg,2);
			}
			knots(i_seg,4)=knots(i_seg,2)*std_knot[1]+knots(i_seg,3)*std_knot[0];
			knots(i_seg,5)=knots(i_seg,2)*std_knot[0]+knots(i_seg,3)*std_knot[1];
			inte_temp[0]=(type_real)1;
			for(i=0;i<=i_seg;++i) inte_temp[0]*=knots(i,6);
			while(++i_seg<this->dim){
				knots(i_seg,2)=x_min(i_seg,0);
				knots(i_seg,3)=knots(i_seg,2)+knots(i_seg,1);
				knots(i_seg,4)=knots(i_seg,2)*std_knot[1]+knots(i_seg,3)*std_knot[0];
				knots(i_seg,5)=knots(i_seg,2)*std_knot[0]+knots(i_seg,3)*std_knot[1];
				knots(i_seg,6)=knots(i_seg,1);
				inte_temp[0]*=knots(i_seg,6);
			}
		}
		//exit of the iteration for segments:
		GTS_iter_end_seg:
		i_seg=1;
		for(i=0;i<this->dim;++i) i_seg*=2;
		res*=(para_in.coef_PDF/(type_real_prob)i_seg);
		return res<(type_real_prob)0 ? (type_real_prob)0 :
			(res>(type_real_prob)1 ? (type_real_prob)1 : res);
	}
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	inline
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_comp
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::CF
	(Type_MatTemp<type_real> const & x) const {
		if(!this->Is_type_res(x)) return type_comp();
		return exp(type_comp(-(x.transpose()*para.Var*x)(0,0)/(type_real)2,(x.transpose()*para.Exp)(0,0)));
	}
	
	//other functions:
	
	template<typename RealType,typename ProbRealType,typename RandUIntType>
	typename Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::type_stat
	Distr_Normal_multi<RealType,ProbRealType,RandUIntType>::Set_para
	(type_res const & mean,type_matrix const & var){
		if(mean.cols()==0 || mean.rows()==0) return (type_stat)1; //illegal mean;
		Type_UInt dim_new=mean.rows(); //new dim;
		if(var.rows()!=dim_new || var.cols()!=dim_new ||
			!var.isApprox(var.transpose().eval(),para_in.pre)){
			return (type_stat)2; //illegal var;
		}
		//eigen decomposition of var:
		Eigen::SelfAdjointEigenSolver<type_matrix> var_eigen(var);
		type_matrix var_eigen_val=var_eigen.eigenvalues();
		//calculate the max component-wise var:
		type_real var_max=std::max(abs(var_eigen_val(0,0)),abs(var_eigen_val(dim_new-1,0)));
		//var is zero ?:
		if(var_max<=para_in.pre){
			this->stat_category=type_this::stat_category_discrete;
			this->dim=dim_new;
			para.Exp=mean;
			para.Var=type_matrix::Zero(this->dim,this->dim);
			para_in.Var_tran=type_matrix();
			para_in.rank=(Type_UInt)0;
			this->stat_able=(unsigned char)9;
		} else {
			type_real var_eq0=para_in.pre*var_max; //approx 0 of var;
			//var is not semi-positive ?:
			if(var_eigen_val(0,0)<-var_eq0){
				return (type_stat)2; //illegal var;
			}
			this->dim=dim_new; //set dim;
			para.Exp=mean; //set para.Exp;
			Type_UInt i=0;
			para_in.rank=this->dim;
			while(i<this->dim && var_eigen_val(i,0)<=var_eq0){
				var_eigen_val(i,0)=(type_real)0;
				--(para_in.rank);
				++i;
			}
			//var is not of full rank ?:
			if((para_in.rank)<this->dim){
				this->stat_category=type_this::stat_category_singular;
				para.Var=Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val);
				para.Var=var_eigen.eigenvectors()*para.Var*var_eigen.eigenvectors().transpose();
				para_in.Var_tran.resize(this->dim,para_in.rank);
				Type_UInt i0=this->dim-para_in.rank;
				for(i=0;i<para_in.rank;++i){
					para_in.Var_tran.col(i)=
						var_eigen.eigenvectors().col(i0+i)*sqrt(var_eigen_val(i0+i,0));
				}
				this->stat_able=(unsigned char)9;
			} else {
				this->stat_category=type_this::stat_category_continuous;
				para.Var=var;
				para_in.Var_tran=var_eigen.eigenvectors()*
					Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val.cwiseSqrt());
				para_in.Var_inv=var_eigen.eigenvectors()*
					Eigen::DiagonalMatrix<type_real,Eigen::Dynamic>(var_eigen_val.cwiseInverse())*
					var_eigen.eigenvectors().transpose();
				para_in.coef_PDF=
					pow((type_real_prob)1/(type_real_prob)LZ_DEF_const_math_sqrt_2pi,(type_real_prob)(this->dim))/
					sqrt(var_eigen_val.prod());
				this->stat_able=(unsigned char)15;
			}
		}
		return (type_stat)0; //succeed;
	}
//}(class Distr_Normal_multi)
//(Multi-variate normal distribution)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_distribution_Distr_Normal_multi_CPP
#endif //#if LZ_DEF_LZ_H_math_src_distribution_distr!=202105L