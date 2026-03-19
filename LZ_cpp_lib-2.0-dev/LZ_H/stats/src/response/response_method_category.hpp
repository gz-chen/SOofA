#if LZ_DEF_LZ_H_stats_response!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_response_response_method_category
#define LZ_DEF_LZ_H_stats_src_response_response_method_category 202105L

#include<Eigen/Eigenvalues>

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//quadratic discrimination analysis{
template<typename OutType=Type_UInt,typename InRealType=Type_Real>
class RespCat_QDA:
public Response_Category_RealVec<OutType,InRealType>{
	public:
		//Type: this class:
		typedef RespCat_QDA<OutType,InRealType> type_this;
		//Type: type of the base class:
		typedef Response_Category_RealVec<OutType,InRealType> type_base;
		//Type: input:
		typedef typename type_base::type_in type_in;
		//Type: output:
		typedef typename type_base::type_out type_out;
		//Type: weight-able output:
		typedef typename type_base::type_out_weight type_out_weight;
		//Type: the weights:
		typedef typename type_base::type_weight type_weight;
		//Type: the real number:
		typedef typename type_base::type_real type_real; //the same as type_weight;
		//Type: the type of state:
		typedef typename type_base::type_stat type_stat;
		//Type: the sequence of inputs:
		typedef typename type_base::type_in_seq type_in_seq;
		//Type: the sequence of outputs:
		typedef typename type_base::type_out_seq type_out_seq;
		//Type: the sequence of weights:
		typedef typename type_base::type_weight_seq type_weight_seq;
		//Type: the matrix:
		typedef Eigen::MatrixX<type_real> type_matrix;
		
		//state values:
		static type_stat constexpr stat_mode_QDA=0;
		static type_stat constexpr stat_mode_condQDA=1;
		static type_stat constexpr stat_mode_distMah=2;
		static type_stat constexpr stat_mode_default=0;
		
		//Constructor:
		RespCat_QDA(void);
		RespCat_QDA(Type_UInt const & dim_input);
		RespCat_QDA(type_out const & num_category,Type_UInt const & dim_input);
		RespCat_QDA(type_this const & model);
		
		//Destructor:
		~RespCat_QDA(void);
		
		//virtual: from type_base:
		virtual type_out Get_response(type_in const & x) const;
		virtual type_out_weight Get_response_wt(type_in const & x) const;
		virtual type_out Get_response_by_wt(type_out_weight const & resp) const;
		virtual type_stat Cal_train(type_out_seq const & Y,type_in_seq const & X);
		virtual type_stat Cal_train
		(type_out_seq const & Y,type_in_seq const & X,type_weight_seq const & wt);
		virtual void Set_reset(void);
		
		//Set: the mode of this classifier:
		type_stat Set_mode(type_stat const mode_new=stat_mode_QDA);
	protected:
		//class: the parameter:
		class type_para{
			public:
				std::vector<type_in> Exp,sum1; //the expectations, the sum of inputs;
				//the variances, the sum of "squared" inputs, the inverse of variances:
				std::vector<type_matrix> Var,sum2,Var_inv;
				//the sizes of the categories, 2 times of their log values,
				//the log-det values of this->Var:
				std::vector<type_real> size,logsize2,logdetVar;
				//states: 0: the parameters keep unchanged; 1: the parameters need re-calculation:
				std::vector<type_stat> recal;
				type_real const pre; //the precision;
				
				//Constructor:
				type_para(type_out const & num_category,Type_UInt const & dim_input);
				type_para(type_para const & para)=default;
				
				//Set: reset:
				void Set_reset(type_this const & model);
				//Cal: calculate other members from sum1, sum2 and size:
				type_stat Cal_para(type_this const & model);
		};
		
		type_para cat_info; //the information of the categories;
		type_real size; //the total sample size;
		//the mode of this classifier:
		//0(default): QDA (Bayesian classification with normal populations with different variances);
		//1: conditional QDA: QDA, but using the same marginal probabilities for each category;
		//2: the distance discrimination using the Mahalanobis distance:
		type_stat mode;
		
		//Get: the criterion value for none-weight-able prediction:
		type_real Get_crit(type_out const & y,type_in const & x) const; //without check for y and x;
}; //class RespCat_QDA;
//}(quadratic discrimination analysis)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/response/RespCat_QDA.cpp"

#endif //#ifndef LZ_DEF_LZ_H_stats_src_response_response_method_category
#endif //#if LZ_DEF_LZ_H_stats_response!=202105L