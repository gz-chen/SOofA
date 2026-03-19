#if LZ_DEF_LZ_H_math_distribution!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_distribution_distr
#define LZ_DEF_LZ_H_math_src_distribution_distr 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//declarations of distributions{
//General discrete distribution through index:
template<typename UIntType,typename ProbRealType,typename RandUIntType>
class Distr_Discrete_index;
#ifdef LZ_DEF_extLIB_Eigen
//Discrete distribution for real vectors:
template<typename RealType,typename ProbRealType,typename RandUIntType>
class Distr_Discrete_Real_Vec;
//Multi-variate normal distribution:
template<typename RealType,typename ProbRealType,typename RandUIntType>
class Distr_Normal_multi;
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(declarations of distributions)
/****************************************************************************************************/
/****************************************************************************************************/
//General discrete distribution through index{
template<typename UIntType=Type_UInt,typename ProbRealType=Type_Real,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distr_Discrete_index: public Distribution_Real<UIntType,ProbRealType,RandUIntType>{
	LZ_DEF_func_check_traits(std::is_integral<UIntType>::value);
	public:
		//Type: this class:
		typedef Distr_Discrete_index<UIntType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef Distribution_Real<UIntType,ProbRealType,RandUIntType> type_base;
		//Type: type of the sample:
		typedef typename type_base::type_res type_res;
		//Type: type of the real number for probability:
		typedef typename type_base::type_real_prob type_real_prob;
		//Type: type of the real number:
		typedef typename type_base::type_real type_real;
		//Type: type of the complex number:
		typedef typename type_base::type_comp type_comp;
		//Type: StatType:
		typedef typename type_base::type_stat type_stat;
		//Type: type of the sample sequence:
		typedef typename type_base::type_res_seq type_res_seq;
		//Type: type of the std random engine:
		typedef typename type_base::type_rand_engine type_rand_engine;
		//Type: the abilities of this class:
		typedef typename type_base::type_able type_able;
		//Type: type of the sequence of probabilities:
		typedef Type_ArrayTemp<type_real_prob> type_real_prob_seq;
		
		//Constructor:
		Distr_Discrete_index(void);
		Distr_Discrete_index(type_res const & number);
		Distr_Discrete_index(type_real_prob_seq const & prob);
		Distr_Discrete_index(type_this const & distr);
		
		//Destructor:
		~Distr_Discrete_index(void){}
		
		//operator:
		type_this& operator =(type_this const & distr);
		
		//virtual: from type_base:
		virtual type_res rand(type_rand_engine & engine) const;
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const;
		virtual type_real_prob PDF(type_res const & x) const;
		virtual type_real_prob CDF(type_real const & x) const;
		virtual type_comp CF(type_real const & x) const;
		virtual type_res quantile(type_real_prob const & prob) const;
		
		//Get: the size of the support:
		type_res Get_size() const {
			return size;
		}
		//Get: the probabilities:
		type_real_prob_seq Get_prob() const {
			return pdf;
		}
		type_real_prob Get_prob(Type_UInt const & num) const {
			return num<size ? pdf[num] : (type_real_prob)0;
		}
		//Set: the probabilities:
		type_stat Set_prob(type_res const & number);
		type_stat Set_prob(type_real_prob_seq const & prob);
		//Get: the definition of CDF:
		type_stat Get_CDF_def() const {
			return CDF_def;
		}
		//Set: the definition of CDF:
		type_stat Set_CDF_def(type_stat const def=(type_stat)0){
			switch(def){
				case (type_stat)0: case (type_stat)1:
					CDF_def=def;
					return (type_stat)0;
				default:
					return (type_stat)1;
			}
		}
	protected:
		type_res size; //the number of possible indices, 0<size<infinity;
		type_real_prob_seq pdf,cdf; //the probabilities;
		type_stat CDF_def; //the definition of CDF: 0: P(X<x); 1: P(X<=x);
		std::uniform_real_distribution<type_real_prob> mutable uniform; //the U(0,1) variable;
}; //class Distr_Discrete_index;
//}(General discrete distribution through index)
/****************************************************************************************************/
//Uniform on a 1-dim interval{

//}(Uniform on a 1-dim interval)
/****************************************************************************************************/
//Discrete distribution for real vectors{
#ifdef LZ_DEF_extLIB_Eigen
template<typename RealType=Type_Real,
	typename ProbRealType=
	typename std::conditional<std::is_floating_point<RealType>::value,RealType,Type_Real>::type,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distr_Discrete_Real_Vec: public Distribution_Real_Vec<RealType,ProbRealType,RandUIntType>{
	public:
		//Type: this class:
		typedef Distr_Discrete_Real_Vec<RealType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef Distribution_Real_Vec<RealType,ProbRealType,RandUIntType> type_base;
		//Type: type of the sample:
		typedef typename type_base::type_res type_res;
		//Type: type of the real number for probability:
		typedef typename type_base::type_real_prob type_real_prob;
		//Type: type of the real number:
		typedef typename type_base::type_real type_real;
		//Type: type of the complex number:
		typedef typename type_base::type_comp type_comp;
		//Type: StatType:
		typedef typename type_base::type_stat type_stat;
		//Type: type of the sample sequence:
		typedef typename type_base::type_res_seq type_res_seq;
		//Type: type of the std random engine:
		typedef typename type_base::type_rand_engine type_rand_engine;
		//Type: the abilities of this class:
		typedef typename type_base::type_able type_able;
		//Type: class Distr_Discrete_index:
		typedef Distr_Discrete_index<Type_UInt,ProbRealType,RandUIntType> type_distr_index;
		//Type: type of the sequence of probabilities:
		typedef typename type_distr_index::type_real_prob_seq type_real_prob_seq;
		
		//Constructor:
		Distr_Discrete_Real_Vec(void);
		Distr_Discrete_Real_Vec(type_res_seq const & support);
		Distr_Discrete_Real_Vec(type_res_seq const & support,type_real_prob_seq const & prob);
		Distr_Discrete_Real_Vec(type_this const & distr);
		
		//Destructor:
		virtual ~Distr_Discrete_Real_Vec(void);
		
		//operator:
		type_this& operator =(type_this const & distr);
		
		//virtual: from type_base:
		virtual type_res rand(type_rand_engine & engine) const;
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const;
		virtual type_real_prob PDF(type_res const & x) const;
		virtual type_real_prob CDF(Type_MatTemp<type_real> const & x) const;
		virtual type_comp CF(Type_MatTemp<type_real> const & x) const;
		
		//Get: the size of the support:
		Type_UInt Get_size() const {
			return distr_index.Get_size();
		}
		//Get: the support:
		type_res_seq Get_support() const {
			return supp;
		}
		type_res Get_support(Type_UInt const & num) const {
			return num<distr_index.Get_size() ? supp[num] : type_res();
		}
		//Get: the probabilities:
		type_real_prob_seq Get_prob() const {
			return distr_index.Get_prob();
		}
		type_real_prob Get_prob(Type_UInt const & num) const {
			return distr_index.Get_prob(num);
		}
		//Get: the definition of CDF:
		type_stat Get_CDF_def() const {
			return distr_index.Get_CDF_def();
		}
		//Set: the definition of CDF:
		type_stat Set_CDF_def(type_stat const def=(type_stat)0){
			return distr_index.Set_CDF_def(def);
		}
		//Set: precision:
		type_stat Set_precision(type_real const precis=(type_real)(-1)){
			if(precis<(type_real)0){
				pre=(type_real)LZ_DEF_const_default_precision;
				return (type_stat)1;
			}
			pre=precis;
			return (type_stat)0;
		}
	protected:
		type_res_seq supp; //the support;
		type_distr_index distr_index; //the distribution of the indices;
		type_real pre; //the precision;
}; //class Distr_Discrete_Real_Vec;
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(Discrete distribution for real vectors)
/****************************************************************************************************/
//Multi-variate normal distribution{
#ifdef LZ_DEF_extLIB_Eigen
template<typename RealType=Type_Real,typename ProbRealType=RealType,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distr_Normal_multi: public Distribution_Real_Vec<RealType,ProbRealType,RandUIntType>{
	LZ_DEF_func_check_traits(std::is_floating_point<RealType>::value);
	public:
		//Type: this class:
		typedef Distr_Normal_multi<RealType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef Distribution_Real_Vec<RealType,ProbRealType,RandUIntType> type_base;
		//Type: type of the sample:
		typedef typename type_base::type_res type_res;
		//Type: type of the real number for probability:
		typedef typename type_base::type_real_prob type_real_prob;
		//Type: type of the real number:
		typedef typename type_base::type_real type_real;
		//Type: type of the complex number:
		typedef typename type_base::type_comp type_comp;
		//Type: StatType:
		typedef typename type_base::type_stat type_stat;
		//Type: type of the sample sequence:
		typedef typename type_base::type_res_seq type_res_seq;
		//Type: type of the std random engine:
		typedef typename type_base::type_rand_engine type_rand_engine;
		//Type: the abilities of this class:
		typedef typename type_base::type_able type_able;
		//Type: matrix:
		typedef Type_MatTemp<type_real> type_matrix;
		
		//class: parameters:
		class type_para{
			public:
				type_para()=default;
				//expectation and variance:
				type_res Exp;
				type_matrix Var;
		};
		
		//Constructor:
		Distr_Normal_multi();
		Distr_Normal_multi(Type_UInt const & dimension);
		Distr_Normal_multi(type_res const & mean,type_real var=(type_real)1);
		Distr_Normal_multi(type_res const & mean,type_matrix const & var);
		Distr_Normal_multi(type_this const & normal);
		
		//Destructor:
		virtual ~Distr_Normal_multi(void);
		
		//operator:
		type_this& operator =(type_this const & distr);
		
		//virtual: from type_base:
		virtual type_res rand(type_rand_engine & engine) const;
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const;
		virtual type_real_prob PDF(type_res const & x) const;
		virtual type_real_prob CDF(Type_MatTemp<type_real> const & x) const;
		virtual type_comp CF(Type_MatTemp<type_real> const & x) const;
		
		//Get: para:
		type_para Get_para() const {
			return para;
		}
		//Get: the expectation:
		type_res Get_para_Exp() const {
			return para.Exp;
		}
		//Get: the variance:
		type_matrix Get_para_Var() const {
			return para.Var;
		}
		//Get: the rank:
		Type_UInt Get_para_rank() const {
			return para_in.rank;
		}
		//Set: precision:
		type_stat Set_precision(type_real const precis=(type_real)(-1)){
			if(precis<=(type_real)0 || precis>=(type_real)1){
				para_in.pre=(type_real)LZ_DEF_const_default_precision;
				return (type_stat)1;
			}
			para_in.pre=precis;
			return (type_stat)0;
		}
		//Set: parameter (Expectation and Variance):
		type_stat Set_para(type_res const & mean,type_matrix const & var);
		type_stat Set_para(type_para const & para_new){
			return this->Set_para(para_new.Exp,para_new.Var);
		}
	protected:
		//class: inside parameters:
		class type_para_inside{
			public:
				type_para_inside()=default;
				//Var_tran= non-zero cols of U*D^(1/2) (Var=U*D*U.transpose()),Var^(-1):
				type_matrix Var_tran,Var_inv;
				//(2*pi)^(-p/2)*(det(Var))^(-1/2):
				type_real_prob coef_PDF;
				//the precision:
				type_real pre;
				//the rank:
				Type_UInt rank;
		};
		
		//members:
		type_para para;
		type_para_inside para_in;
		std::normal_distribution<type_real> mutable std_normal;
}; //class Distr_Normal_multi;
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(Multi-variate normal distribution)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//implementation:
#include"LZ_H/math/src/distribution/Distr_Discrete_index.cpp"
#ifdef LZ_DEF_extLIB_Eigen
	#include"LZ_H/math/src/distribution/Distr_Discrete_Real_Vec.cpp"
	#include"LZ_H/math/src/distribution/Distr_Normal_multi.cpp"
#endif

#endif //#ifndef LZ_DEF_LZ_H_math_src_distribution_distr
#endif //#if LZ_DEF_LZ_H_math_distribution!=202105L