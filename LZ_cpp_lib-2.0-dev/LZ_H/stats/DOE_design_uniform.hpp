#ifndef LZ_DEF_LZ_H_stats_DOE_design_uniform
#define LZ_DEF_LZ_H_stats_DOE_design_uniform 202105L

#include"LZ_H/stats/basic.hpp"
#include"LZ_H/math/distribution.hpp"

#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//Iterator for sequential canonical design points{
//the base class: Iterator_design_canon_sequential:
template<typename RealType=Type_Real>
class Iterator_design_canon_sequential{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//Type: this class:
		typedef Iterator_design_canon_sequential<RealType> type_this;
		//Type: the base class:
		typedef void type_base;
		//Type: the integer for size:
		typedef Type_Size type_size;
		//Type: the real number:
		typedef RealType type_real;
		//Type: the matrix:
		typedef Type_MatTemp<type_real> type_matrix;
		//Type: the state:
		typedef Type_Stat type_stat;
		//Type: its iterator:
		typedef type_this iterator;
		//Types: for iterator:
		typedef std::iterator<std::forward_iterator_tag,type_matrix const> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_sequential();
		Iterator_design_canon_sequential(type_size const & dimension);
		Iterator_design_canon_sequential(type_this const & iter);
		Iterator_design_canon_sequential(type_this && iter);
		
		//Destructor:
		virtual ~Iterator_design_canon_sequential()=default;
		
		//operator:
		virtual value_type & operator *() const final;
		virtual value_type * operator ->() const final;
		virtual type_this& operator ++()=0;
		
		//begin-end:
		virtual void Set_begin()=0;
		
		//Is:
		Type_Bool Is_deref() const;
		Type_Bool Is_end() const;
		Type_Bool Is_null() const;
		
		//Get:
		type_size Get_dim() const;
	protected:
		type_size dim;
		type_matrix point; //the current design point;
		type_stat stat; //the state;
		
		//static values:
		static type_stat constexpr val_stat_deref=0;
		static type_stat constexpr val_stat_end=1;
		static type_stat constexpr val_stat_null=10;
		
		//operator:
		type_this& operator =(type_this const & iter);
		type_this& operator =(type_this && iter);
		Type_Bool operator ==(type_this const & iter) const;
}; //class Iterator_design_canon_sequential;
//random U([0,1]^dim) samples: Iterator_design_canon_RandomUniform:
template<typename RealType=Type_Real,
	typename RandUIntType=typename Liuze::math::Type_RandEngine::result_type>
class Iterator_design_canon_RandomUniform:
public Iterator_design_canon_sequential<RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_RandomUniform<RealType,RandUIntType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<RealType> type_base;
		//Type: the integer for size:
		typedef typename type_base::type_size type_size;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the random engine:
		typedef Liuze::math::RandomEngine<RandUIntType> type_rand_eng;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef std::iterator<std::input_iterator_tag,type_matrix const> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_RandomUniform();
		Iterator_design_canon_RandomUniform(type_rand_eng & rand_engine);
		Iterator_design_canon_RandomUniform(type_size const & dimension,type_rand_eng & rand_engine);
		Iterator_design_canon_RandomUniform(type_this const & iter);
		Iterator_design_canon_RandomUniform(type_this && iter);
		
		//Destructor:
		virtual ~Iterator_design_canon_RandomUniform()=default;
		
		//operator:
		type_this& operator =(type_this const & iter);
		type_this& operator =(type_this && iter);
		Type_Bool operator ==(type_this const & iter) const;
		virtual type_this& operator ++();
		
		//begin-end:
		virtual void Set_begin();
		void Set_begin(type_size const & dimension);
		
		//Set:
		void Set_rand_engine(type_rand_eng & rand_engine);
	protected:
		std::uniform_real_distribution<type_real> distr_U; //the uniform random variable generator;
		type_rand_eng * engine; //the random engine;
		
		//Cal:
		void Cal_renew();
}; //class Iterator_design_canon_RandomUniform;
template<typename RealType=Type_Real>
class Iterator_design_canon_GoodPoint:
public Iterator_design_canon_sequential<RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_GoodPoint<RealType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<RealType> type_base;
		//Type: the integer for size:
		typedef typename type_base::type_size type_size;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef typename type_base::iterator_category iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_GoodPoint();
		Iterator_design_canon_GoodPoint(type_matrix const & good_point);
		Iterator_design_canon_GoodPoint(type_this const & iter);
		Iterator_design_canon_GoodPoint(type_this && iter);
		
		//Destructor:
		virtual ~Iterator_design_canon_GoodPoint()=default;
		
		//operator:
		type_this& operator =(type_this const & iter);
		type_this& operator =(type_this && iter);
		Type_Bool operator ==(type_this const & iter) const;
		virtual type_this& operator ++();
		
		//begin-end:
		virtual void Set_begin();
		void Set_begin(type_matrix const & good_point);
		
		//static functions:
		//Get: good point:
		//cyclotomic:
		//prime_number>=2*dimension+3:
		static type_matrix Get_GoodPoint_cyclotomic(type_size const & dimension,type_size prime_number=0);
		//power:
		//When dimension=d, prime_number=p, and modify=0,
		//let q=p^(1/(d+1)), then gp=(q^1,q^2,...,q^d);
		//but when modify=1,
		//gp=(a1,b1,a2,b2,...) with (a1,a2,...,b1,b2,...)=(q^1,q^2,...,q^d),
		//see J. G. Liao (1998) on Journal of Computational and Graphical Statistics.
		static type_matrix Get_GoodPoint_power
			(type_size const & dimension,type_size prime_number=2,type_stat modify=0);
	protected:
		type_matrix gp; //the good point;
}; //class Iterator_design_canon_GoodPoint;
template<typename RealType=Type_Real>
class Iterator_design_canon_WrapAroundShifter:
public Iterator_design_canon_sequential<RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_WrapAroundShifter<RealType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<RealType> type_base;
		//Type: the integer for size:
		typedef typename type_base::type_size type_size;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef std::iterator<std::input_iterator_tag,type_matrix const> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		//Type: for storing iterators:
		typedef Iterator_wrapper<value_type,iterator_category,difference_type> type_iter_store;
		
		//static value:
		static type_size constexpr val_size_infinity=0;
		
		//Constructor:
		Iterator_design_canon_WrapAroundShifter();
		template<typename InIterType>
		Iterator_design_canon_WrapAroundShifter
			(Type_ArrayTemp<type_size> const & list_design_size,
			Type_ArrayTemp<InIterType> const & list_design_begin);
		Iterator_design_canon_WrapAroundShifter(type_this const & iter);
		Iterator_design_canon_WrapAroundShifter(type_this && iter);
		
		//Destructor:
		virtual ~Iterator_design_canon_WrapAroundShifter()=default;
		
		//operator:
		type_this& operator =(type_this const & iter);
		type_this& operator =(type_this && iter);
		Type_Bool operator ==(type_this const & iter) const;
		virtual type_this& operator ++();
		
		//begin-end:
		virtual void Set_begin();
		template<typename InIterType>
		void Set_begin
			(Type_ArrayTemp<type_size> const & list_design_size,
			Type_ArrayTemp<InIterType> const & list_design_begin);
	protected:
		Type_ArrayTemp<type_iter_store> list_iter_begin,list_iter;
		Type_ArrayTemp<type_size> list_size,list_i_size;
		Type_ArrayTemp<type_matrix> list_point_buff; //a buffer for partial sums of component designs;
}; //class Iterator_design_canon_WrapAroundShifter;
//}(Iterator for sequential canonical design points)
/****************************************************************************************************/
//good lattice point{
template<typename RealType=Type_Real,typename IntType=Type_UInt>
class Get_design_GoodLatticePoint{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	LZ_DEF_func_check_traits((std::is_integral<IntType>::value));
	public:
		//type: this class:
		typedef Get_design_GoodLatticePoint<RealType,IntType> type_this;
		//type: integer:
		typedef IntType type_int;
		//type: real number:
		typedef RealType type_real;
		//type: point for factorial design:
		typedef Type_MatTemp<type_int> type_pt_fac;
		//type: point for canonical design:
		typedef Type_MatTemp<type_real> type_pt_can;
		//type: integer for size:
		typedef Type_Size type_size;
		//type: state:
		typedef Type_Stat type_stat;
		
		//static value:
		static type_stat constexpr val_stat_null=0; //size is 0;
		static type_stat constexpr val_stat_invalid=1; //size>=2 and vec_gen is not calculated;
		static type_stat constexpr val_stat_valid=2; //size>=2 and vec_gen is valid;
		
		//Constructor:
		Get_design_GoodLatticePoint();
		Get_design_GoodLatticePoint(type_int const & size_design);
		
		//Is:
		Type_Bool Is_null() const;
		Type_Bool Is_invalid() const;
		Type_Bool Is_valid() const;
		
		//Get: value:
		type_int Get_dim() const;
		type_int Get_size() const;
		Type_ArrayTemp<type_int> Get_unit() const;
		Type_ArrayTemp<type_real> Get_design_unary() const;
		type_pt_fac Get_vec_gen() const;
		type_real Get_val_crit() const;
		
		//Set:
		type_stat Set_size(type_int const & size_design);
		
		//Cal: vec_gen and val_crit:
		//find the best GLP over all vec_gen starting with 1:
		template<typename CritType>
		type_stat Cal_best(type_int const & dim_design,CritType & crit_canon);
		//find the best GLP over all power vec_gen:
		template<typename CritType>
		type_stat Cal_best_power(type_int const & dim_design,CritType & crit_canon);
		//find the best GLP over all power vector_gen:
		//(without checking whether the components of the generating vector are different from each other):
		template<typename CritType>
		type_stat Cal_best_power_nocheck(type_int const & dim_design,CritType & crit_canon);
		//find the best leave-one-out GLP over all power vector_gen:
		//(without checking whether the components of the generating vector are different from each other):
		template<typename CritType>
		type_stat Cal_best_power_LeaveOneOut_nocheck(type_int const & dim_design,CritType & crit_canon);
		
		//Get: design:
		//levels are in {0,...,size-1} and the generating vector is vec_gen:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_fac>::iterator>
		OutIterType Get_design_factor(OutIterType result_begin) const;
		//applying the transformation: x |-> (x+(0.5,...,0.5))/size
		//for all points in the corresponding factor good lattice point design:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_can>::iterator>
		OutIterType Get_design_canon(OutIterType result_begin) const;
		//excluding the point (0,...,0) from the GLP(size,vec_gen)
		//and then subtracting (1,...,1) from the rested points:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_fac>::iterator>
		OutIterType Get_design_factor_LeaveOneOut(OutIterType result_begin) const;
		//applying the transformation: x |-> (x+(0.5,...,0.5))/(size-1)
		//for all points in the corresponding leave-one-out factor good lattice point design:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_can>::iterator>
		OutIterType Get_design_canon_LeaveOneOut(OutIterType result_begin) const;
		
		//Convert:
		type_pt_can Convert_to_canon(type_pt_fac const & point) const;
		type_pt_fac Convert_to_factor(type_pt_can const & point) const;
		
		//static function:
		//Is:
		static Type_Bool Is_null(type_stat const & state);
		static Type_Bool Is_invalid(type_stat const & state);
		static Type_Bool Is_valid(type_stat const & state);
		//Get: design:
		//levels are in {0,...,size_design-1}:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_fac>::iterator>
		static OutIterType Get_design_factor
		(type_int const & size_design,type_pt_fac const & vector_gen,OutIterType result_begin);
		//applying the transformation: x |-> (x+(0.5,...,0.5))/size_design
		//for all points in the corresponding factor good lattice point design:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_can>::iterator>
		static OutIterType Get_design_canon
		(type_int const & size_design,type_pt_fac const & vector_gen,OutIterType result_begin);
		//excluding the point (0,...,0) from the GLP(size_design_full,vector_gen)
		//and then subtracting (1,...,1) from the rested points:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_fac>::iterator>
		static OutIterType Get_design_factor_LeaveOneOut
		(type_int const & size_design_full,type_pt_fac const & vector_gen,OutIterType result_begin);
		//applying the transformation: x |-> (x+(0.5,...,0.5))/(size_design_full-1)
		//for all points in the corresponding leave-one-out factor good lattice point design:
		template<typename OutIterType=typename Type_ArrayTemp<type_pt_can>::iterator>
		static OutIterType Get_design_canon_LeaveOneOut
		(type_int const & size_design_full,type_pt_fac const & vector_gen,OutIterType result_begin);
	protected:
		//val:
		type_int size,dim; //size and dim of this design;
		type_int dim_max; //the max dim under this size;
		type_stat stat; //the state: see the static values val_stat_...;
		Type_ArrayTemp<type_int> unit; //{x : 1<x<size and gcd(x,size)=1};
		Type_ArrayTemp<type_real> crd_can; //coordinates for the canonical design;
		
		//val: results:
		type_pt_fac vec_gen; //the generation vector;
		type_real val_crit; //the value of the criterion;
}; //class Get_design_GoodLatticePoint;
//}(good lattice point)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

//implementation:
#include"LZ_H/stats/src/DOE/DOE_design_uniform.cpp"

#endif //#ifdef LZ_DEF_extLIB_Eigen
#endif //#ifndef LZ_DEF_LZ_H_stats_DOE_design_uniform