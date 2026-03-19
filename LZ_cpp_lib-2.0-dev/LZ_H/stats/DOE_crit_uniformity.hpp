#ifndef LZ_DEF_LZ_H_stats_DOE_crit_uniformity
#define LZ_DEF_LZ_H_stats_DOE_crit_uniformity 202105L

#include"LZ_H/stats/basic.hpp"

#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace stats{
namespace crit{ //a namespace for criteria;

/****************************************************************************************************/
template<typename RealType=Type_Real>
class Discrepancy_Leb_infty{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		typedef Discrepancy_Leb_infty<RealType> type_this; //this class;
		typedef Type_Size type_size;
		typedef RealType type_real; //the real number;
		typedef Type_MatTemp<type_real> type_matrix; //the real number matrix;
		typedef std::function<type_real(type_matrix)> type_cdf; //the CDF;
		typedef Type_Stat type_stat;
		
		//class for CDF(U([0,1]^dim)):
		class CDF_Uniform_std{
			public:
				//Constructor:
				CDF_Uniform_std(type_size dimension=1);
				//CDF:
				type_real operator ()(type_matrix const & x) const;
			protected:
				type_size dim;
		}; //class CDF_Uniform_std;
		
		//static value:
		static type_stat constexpr val_CDF_def_less=0;
		static type_stat constexpr val_CDF_def_lesseq=1;
		
		//Constructor:
		//the target distribution is Uniform([0,1]^dimension):
		Discrepancy_Leb_infty
		(type_size dimension=1,type_real precision=-1,type_stat cdf_def=type_this::val_CDF_def_less);
		//the target distribution function is func_cdf:
		template<typename FuncDistrType>
		Discrepancy_Leb_infty
		(FuncDistrType & func_distr,type_real precision=-1,type_stat cdf_def=type_this::val_CDF_def_less);
		
		//Calculate:
		//this operator returns the L^infty discrepancy
		//between the target CDF and the ECDF of the samples,
		//where for example, for common IID samples, the weights should be 1/(sample size):
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const;
		
		//Get:
		type_stat Get_CDF_def() const;
		type_cdf Get_CDF() const;
		type_real Get_precision() const;
		
		//Set:
		type_stat Set_CDF_def(type_stat cdf_def=type_this::val_CDF_def_less);
		type_cdf Set_CDF(type_size dimension=1);
		template<typename FuncDistrType>
		type_cdf Set_CDF(FuncDistrType & func_distr);
		type_real Set_precision(type_real precision=-1);
	protected:
		type_stat CDF_def; //the definition of CDF: 0: F(x)=P(X<x); 1: F(x)=P(X<=x);
		type_cdf fun_cdf; //CDF of the target distribution;
		type_real pre; //the precision;
}; //class Discrepancy_Leb_infty;
/****************************************************************************************************/
//this is the approximated L^infinity discrepancy towards U([0,1]^dim),
//based on p97 of Hua L.-K. and Wang Y. (1981):
//this criterion has some problems and it should not be used!
template<typename RealType=Type_Real>
class Discrepancy_canon_Leb_infty_approxLambda{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		typedef Discrepancy_canon_Leb_infty_approxLambda<RealType> type_this; //this class;
		typedef Type_Size type_size;
		typedef RealType type_real; //the real number;
		typedef Type_MatTemp<type_real> type_matrix; //the real number matrix;
		
		//Constructor:
		Discrepancy_canon_Leb_infty_approxLambda()=default;
		
		//Calculate:
		//this operator returns the approximated L^infty discrepancy
		//between CDF(U([0,1]^dim)) and the ECDF of the samples,
		//where for example, for common IID samples, the weights should be 1/(sample size):
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const;
}; //class Discrepancy_canon_Leb_infty_approxLambda
/****************************************************************************************************/
template<typename RealType=Type_Real>
class Discrepancy_canon_Leb_2_kernel{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		typedef Discrepancy_canon_Leb_2_kernel<RealType> type_this; //this class;
		typedef Type_Size type_size;
		typedef RealType type_real; //the real number;
		typedef Type_MatTemp<type_real> type_matrix; //the real number matrix;
		typedef Type_Stat type_stat;
		
		//the base class for kernel functions:
		class type_kernel{
			public:
				typedef type_kernel type_this;
				//Constructor:
				type_kernel(type_size const & dimension);
				//Destructor:
				~type_kernel()=default;
				//Get:
				type_size Get_dim() const;
				//the kernel:
				virtual type_real operator ()
				(type_matrix const & x0,type_matrix const & x1) const=0;
				//the integration of the kernel with respect to one of its parameter:
				virtual type_real operator ()(type_matrix const & x) const=0;
				//the integration of the kernel:
				virtual type_real operator ()() const=0;
			protected:
				type_size dim; //the dimension;
		};
		class type_kernel_centered; //centered; constructor: (dimension,weight_projective=1);
		class type_kernel_wrap; //wrap-around; constructor: (dimension,weight_projective=1);
		class type_kernel_mixture; //mixture; constructor: (dimension,weight_projective=1);
		class type_kernel_discrete; //mixture; constructor: (dimension); need other initialisations;
		
		//static value:
		static type_stat constexpr val_stat_ker_centered=11; //centered;
		static type_stat constexpr val_stat_ker_wrap=12; //wrap-around;
		static type_stat constexpr val_stat_ker_mixture=13; //mixture;
		//`static type_stat constexpr val_stat_ker_discrete=14; //discrete;
		
		//Constructor:
		Discrepancy_canon_Leb_2_kernel(type_kernel & kernel);
		Discrepancy_canon_Leb_2_kernel
		(type_size dimension,type_stat kernel,type_real weight_projective=1);
		
		//Destructor:
		~Discrepancy_canon_Leb_2_kernel();
		
		//Calculate:
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const;
		
		//Set:
		type_stat Set_kernel(type_kernel & kernel);
		type_stat Set_kernel(type_size dimension,type_stat kernel,type_real weight_projective=1);
	protected:
		type_kernel * pt_ker; //the pointer to kernel;
		type_stat stat_ker; //the state of kernel; see the static values: val_stat_kernel_...;
		
		//static value:
		static type_stat constexpr val_stat_ker_nonew=0; //the kernel is obtained from outside;
		static type_stat constexpr val_stat_ker_NULL=1; //pt_ker is NULL;
}; //class Discrepancy_canon_Leb_2_kernel;
/****************************************************************************************************/

} //namespace crit;
} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/DOE/DOE_crit_uniformity.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_stats_DOE_crit_uniformity