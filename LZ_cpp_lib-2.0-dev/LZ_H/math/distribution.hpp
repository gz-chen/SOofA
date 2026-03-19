/****************************************************************************************************/
/*
This file defines the base classes of random engine and distribution.
*/
/****************************************************************************************************/

#ifndef LZ_DEF_LZ_H_math_distribution
#define LZ_DEF_LZ_H_math_distribution 202105L

#include"LZ_H/math/function.hpp"
#include<random>
#ifdef LZ_DEF_extLIB_Eigen
	#include<Eigen/Eigenvalues>
#endif

//the following part deals with the std random engines.
namespace Liuze{
namespace math{

/****************************************************************************************************/
//typedef: default types for std random engine{
typedef std::default_random_engine Type_StdRandEngine;
typedef std::uint_fast32_t Type_StdRandEngine_UInt;
typedef unsigned long long Type_StdRandEngine_StepLen;
//}(typedef: default types for std random engine)
/****************************************************************************************************/
//traits: is_std_rand_engine{
template<typename Tp>
struct is_std_rand_engine: public std::false_type{
};
//following are the random engines in C++ std lib random.
template<typename UIntType,UIntType a,UIntType c,UIntType m>
struct is_std_rand_engine<std::linear_congruential_engine<UIntType,a,c,m> >:
public std::true_type{
};
template<typename UIntType,size_t w,size_t n,size_t m,size_t r,UIntType a,size_t u,
	UIntType d,size_t s,UIntType b,size_t t,UIntType c,size_t l,UIntType f>
struct is_std_rand_engine<
	std::mersenne_twister_engine<UIntType,w,n,m,r,a,u,d,s,b,t,c,l,f> >:
public std::true_type{
};
template<typename UIntType,size_t w,size_t s,size_t r>
struct is_std_rand_engine<std::subtract_with_carry_engine<UIntType,w,s,r> >:
public std::true_type{
};
template<typename RandomNumberEngine,size_t p,size_t r>
struct is_std_rand_engine<std::discard_block_engine<RandomNumberEngine,p,r> >:
public std::true_type{
};
template<typename RandomNumberEngine,size_t w,typename UIntType>
struct is_std_rand_engine<std::independent_bits_engine<RandomNumberEngine,w,UIntType> >:
public std::true_type{
};
template<typename RandomNumberEngine,size_t k>
struct is_std_rand_engine<std::shuffle_order_engine<RandomNumberEngine,k> >:
public std::true_type{
};
//}(traits: is_std_rand_engine)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//the following part defines the base class RandomEngine and some basic random engines.
namespace Liuze{
namespace math{

/****************************************************************************************************/
//Random engine base class{
template<typename UIntType=Type_UInt>
class RandomEngine{
	LZ_DEF_func_check_traits(std::is_unsigned_integral<UIntType>::value);
	public:
		//Type: this class:
		typedef RandomEngine<UIntType> type_this;
		//Type: base class:
		typedef void type_base;
		//Type: result type:
		typedef UIntType type_res;
		typedef type_res result_type;
		
		//Constructor:
		RandomEngine(void)=default;
		
		//Destructor:
		virtual ~RandomEngine(void)=default;
		
		//Get: the minimum value:
		virtual result_type min() const =0;
		//Get: the maximum value:
		virtual result_type max() const =0;
		//Set: the seed:
		virtual void seed()=0;
		virtual void seed(result_type const & seed_val)=0;
		//Get: a random result_type value:
		virtual result_type operator ()()=0;
		//Set: advance internal state:
		virtual void discard(Type_StdRandEngine_StepLen const & step)=0;
}; //class RandomEngine;
//}(Random engine base class)
/****************************************************************************************************/
//traits: about class RandomEngine{
//is derived from class RandomEngine:
template<typename Tp>
struct is_derived_RandomEngine:
public std::is_base_of<RandomEngine<typename Tp::result_type>,Tp>::type{
};
//is_std_rand_engine or is_derived_RandomEngine:
template<typename Tp>
struct is_rand_engine:
public std::__or_<is_std_rand_engine<Tp>,is_derived_RandomEngine<Tp> >{
};
//}(traits: about class RandomEngine)
/****************************************************************************************************/
//Random engine class from std random engine{
template<typename StdRandEngineType=Type_StdRandEngine>
class RandomEngine_std: public RandomEngine<typename StdRandEngineType::result_type>{
	LZ_DEF_func_check_traits(is_std_rand_engine<StdRandEngineType>::value);
	public:
		//Type: this class:
		typedef RandomEngine_std<StdRandEngineType> type_this;
		//Type: base class:
		typedef RandomEngine<typename StdRandEngineType::result_type> type_base;
		//Type: result type:
		typedef typename type_base::type_res type_res;
		typedef typename StdRandEngineType::result_type result_type;
		//Type: std random engine:
		typedef StdRandEngineType type_rand_engine;
		
		//public member:
		type_rand_engine rand_engine;
		
		//Constructor:
		RandomEngine_std():
		type_base(),rand_engine(){
		}
		RandomEngine_std(result_type const & seed_val):
		type_base(),rand_engine(seed_val){
		}
		template<typename Sseq>
		RandomEngine_std(Sseq & q):
		type_base(),rand_engine(q){
		}
		RandomEngine_std(type_this const & engine):
		type_base((type_base&)engine),rand_engine(engine.rand_engine){
		}
		RandomEngine_std(type_rand_engine const & engine):
		type_base(),rand_engine(engine){
		}
		
		//Destructor:
		virtual ~RandomEngine_std(void)=default;
		
		//virtual: from type_base:
		virtual result_type min() const final {
			return rand_engine.min();
		}
		virtual result_type max() const final {
			return rand_engine.max();
		}
		virtual void seed() final {
			rand_engine.seed();
		}
		virtual void seed(result_type const & seed_val) final {
			rand_engine.seed(seed_val);
		}
		virtual result_type operator ()() final {
			return rand_engine();
		}
		virtual void discard(Type_StdRandEngine_StepLen const & step) final {
			rand_engine.discard(step);
		}
		Type_Bool operator ==(type_this const & engine) const {
			return rand_engine==engine.rand_engine;
		}
}; //class RandomEngine_std;
//}(Random engine class from std random engine)
/****************************************************************************************************/
//Random engine class from the C-style function rand{
template<typename UIntType=Type_UInt>
class RandomEngine_cstd: public RandomEngine<UIntType>{
	public:
		//Type: this class:
		typedef RandomEngine_cstd<UIntType> type_this;
		//Type: base class:
		typedef RandomEngine<UIntType> type_base;
		//Type: result type:
		typedef typename type_base::type_res type_res;
		typedef typename type_base::result_type result_type;
		
		//Constructor:
		RandomEngine_cstd():
		type_base(){
		}
		RandomEngine_cstd(result_type const & seed_val):
		type_base(){
			srand((unsigned int)seed_val);
		}
		
		//Destructor:
		virtual ~RandomEngine_cstd(void)=default;
		
		//virtual: from type_base:
		virtual result_type min() const final {
			return (result_type)0;
		}
		virtual result_type max() const final {
			return (result_type)RAND_MAX;
		}
		virtual void seed() final {
			srand((unsigned int)0);
		}
		virtual void seed(result_type const & seed_val) final {
			srand((unsigned int)seed_val);
		}
		virtual result_type operator ()() final {
			return (result_type)rand();
		}
		virtual void discard(Type_StdRandEngine_StepLen const & step) final {
			for(Type_StdRandEngine_StepLen i=0;i<step;++i) rand();
		}
		Type_Bool operator ==(type_this const & engine) const {
			return true;
		}
}; //class RandomEngine_cstd;
//}(Random engine class from the C-style function rand)
/****************************************************************************************************/
//typedef: default types from class RandomEngine{
typedef RandomEngine_std<Type_StdRandEngine> Type_RandEngine;
typedef RandomEngine_cstd<Type_UInt> Type_RandEngine_Cstd;
//}(typedef: default types from class RandomEngine)
/****************************************************************************************************/
//Random engine class by reference to other random engine{
template<typename RandEngineType=Type_RandEngine>
class RandomEngine_ref: public RandomEngine<typename Type_RandEngine::result_type>{
	LZ_DEF_func_check_traits(is_rand_engine<RandEngineType>::value);
	public:
		//Type: this class:
		typedef RandomEngine_ref<RandEngineType> type_this;
		//Type: base class:
		typedef RandomEngine<typename Type_RandEngine::result_type> type_base;
		//Type: result type:
		typedef typename type_base::type_res type_res;
		typedef typename type_base::result_type result_type;
		//Type: std random engine:
		typedef RandEngineType type_rand_engine;
		
		//public member:
		type_rand_engine& rand_engine;
		
		//Constructor:
		RandomEngine_ref(type_this const & engine):
		type_base((type_base&)engine),rand_engine(engine.rand_engine){
		}
		RandomEngine_ref(type_rand_engine const & engine):
		type_base(),rand_engine(engine){
		}
		
		//Destructor:
		virtual ~RandomEngine_ref(void)=default;
		
		//virtual: from type_base:
		virtual result_type min() const final {
			return (result_type)(rand_engine.min());
		}
		virtual result_type max() const final {
			return (result_type)(rand_engine.max());
		}
		virtual void seed() final {
			rand_engine.seed();
		}
		virtual void seed(result_type const & seed_val) final {
			rand_engine.seed(seed_val);
		}
		virtual result_type operator ()() final {
			return (result_type)(rand_engine());
		}
		virtual void discard(Type_StdRandEngine_StepLen const & step) final {
			rand_engine.discard(step);
		}
		Type_Bool operator ==(type_this const & engine) const {
			return (&rand_engine)==(&(engine.rand_engine));
		}
}; //class RandomEngine_ref;
//}(Random engine class by reference to other random engine)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//the following part defines the base class Distribution and some derived base classes.
namespace Liuze{
namespace math{

/****************************************************************************************************/
//Distribution base class{
template<typename ResType,typename ProbRealType=Type_Real,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distribution{
	LZ_DEF_func_check_traits(std::is_floating_point<ProbRealType>::value);
	LZ_DEF_func_check_traits(std::is_unsigned_integral<RandUIntType>::value);
	public:
		//Type: this class:
		typedef Distribution<ResType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef void type_base;
		//Type: type of the sample:
		typedef ResType type_res;
		//Type: type of the real number for probability:
		typedef ProbRealType type_real_prob;
		//Type: StatType:
		typedef Type_Stat type_stat;
		//Type: type of the sample sequence:
		typedef Type_ArrayTemp<type_res> type_res_seq;
		//Type: type of the std random engine:
		typedef RandomEngine<RandUIntType> type_rand_engine;
		
		//class: ablilities of this distribution:
		class type_able{
			public:
				//able_rand,able_PDF,able_CDF,able_CF,able_quantile;
				//stored by bits, from low to high:
				//1: true, 0: false.
				unsigned char able_code;
				
				type_able(){able_code=(unsigned char)0;}
				type_able(unsigned char const & code){able_code=code;}
				type_able(type_able const & a){able_code=a.able_code;}
				
				type_able& operator =(type_able const & a){
					able_code=a.able_code;
					return *this;
				}
				type_able& operator =(unsigned char const & code){
					able_code=code;
					return *this;
				}
				
				bool able_rand() const {return (bool)(able_code&(unsigned char)1);}
				bool able_PDF() const {return (bool)(able_code&(unsigned char)(1<<1));}
				bool able_CDF() const {return (bool)(able_code&(unsigned char)(1<<2));}
				bool able_CF() const {return (bool)(able_code&(unsigned char)(1<<3));}
				bool able_quantile() const {return (bool)(able_code&(unsigned char)(1<<4));}
				
				void Set_able_rand(bool const & value){
					able_code= value ? able_code|(unsigned char)1 : able_code&(unsigned char)(~0^1);
				}
				void Set_able_PDF(bool const & value){
					able_code= value ? able_code|(unsigned char)(1<<1) : able_code&(unsigned char)(~0^(1<<1));
				}
				void Set_able_CDF(bool const & value){
					able_code= value ? able_code|(unsigned char)(1<<2) : able_code&(unsigned char)(~0^(1<<2));
				}
				void Set_able_CF(bool const & value){
					able_code= value ? able_code|(unsigned char)(1<<3) : able_code&(unsigned char)(~0^(1<<3));
				}
				void Set_able_quantile(bool const & value){
					able_code= value ? able_code|(unsigned char)(1<<4) : able_code&(unsigned char)(~0^(1<<4));
				}
		};
		
		//state values:
		static type_stat constexpr stat_category_discrete=0;
		static type_stat constexpr stat_category_continuous=1;
		static type_stat constexpr stat_category_singular=2;
		
		//Constructor:
		Distribution(type_stat const & distr_category):
		stat_able((unsigned char)0),stat_category(distr_category){
		}
		Distribution(type_stat const & distr_category,type_able const & able):
		stat_able(able),stat_category(distr_category){
		}
		Distribution(type_this const & distr):
		stat_able(distr.stat_able),stat_category(distr.stat_category){
		}
		//Destructor:
		virtual ~Distribution(void){}
		
		//operator =:
		type_this& operator =(type_this const & distr){
			if(&distr==this) return *this;
			stat_category=distr.stat_category;
			stat_able=distr.stat_able;
			return *this;
		}

		//Cal: generating a random sample:
		virtual type_res rand(type_rand_engine & engine) const =0;
		//Cal: generating a random sample, using operator ():
		virtual type_res operator ()(type_rand_engine & engine) const final {
			return this->rand(engine);
		}
		//Cal: generating a vector of random samples:
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const =0;
		//Cal: calculating the PDF:
		virtual type_real_prob PDF(type_res const & x) const =0;
		
		//Get: type of this distribution:
		type_stat Get_stat_category(void) const {return stat_category;}
		//Get: ables:
		bool Get_able_rand(void) const {return stat_able.able_rand();}
		bool Get_able_PDF(void) const {return stat_able.able_PDF();}
		bool Get_able_CDF(void) const {return stat_able.able_CDF();}
		bool Get_able_CF(void) const {return stat_able.able_CF();}
		bool Get_able_quantile(void) const {return stat_able.able_quantile();}
		type_able Get_able(void) const {return stat_able;}
	protected:
		type_stat stat_category; //0:discrete;1:continuous;2:singular.
		type_able stat_able; //able_rand,able_PDF,able_CDF,able_CF,able_quantile.
}; //class Distribution<ResType>;
//static constexpr in class Distribution:
template<typename ResType,typename ProbRealType,typename RandUIntType>
typename Distribution<ResType,ProbRealType,RandUIntType>::type_stat constexpr
Distribution<ResType,ProbRealType,RandUIntType>::stat_category_discrete;
template<typename ResType,typename ProbRealType,typename RandUIntType>
typename Distribution<ResType,ProbRealType,RandUIntType>::type_stat constexpr
Distribution<ResType,ProbRealType,RandUIntType>::stat_category_continuous;
template<typename ResType,typename ProbRealType,typename RandUIntType>
typename Distribution<ResType,ProbRealType,RandUIntType>::type_stat constexpr
Distribution<ResType,ProbRealType,RandUIntType>::stat_category_singular;
//}(Distribution base class)
/****************************************************************************************************/
//Distribution of the 1-dim real number{
template<typename RealType=Type_Real,
	typename ProbRealType=
	typename std::conditional<std::is_floating_point<RealType>::value,RealType,Type_Real>::type,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distribution_Real: public Distribution<RealType,ProbRealType,RandUIntType>{
	LZ_DEF_func_check_traits(std::is_arithmetic<RealType>::value);
	public:
		//Type: this class:
		typedef Distribution_Real<RealType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef Distribution<RealType,ProbRealType,RandUIntType> type_base;
		//Type: type of the sample:
		typedef typename type_base::type_res type_res;
		//Type: type of the real number for probability:
		typedef typename type_base::type_real_prob type_real_prob;
		//Type: type of the real number:
		typedef
			typename std::conditional<std::is_floating_point<type_res>::value,type_res,type_real_prob>::type
			type_real;
		//Type: type of the complex number:
		typedef Type_CompTemp<type_real> type_comp;
		//Type: StatType:
		typedef typename type_base::type_stat type_stat;
		//Type: type of the sample sequence:
		typedef typename type_base::type_res_seq type_res_seq;
		//Type: type of the std random engine:
		typedef typename type_base::type_rand_engine type_rand_engine;
		//Type: the abilities of this class:
		typedef typename type_base::type_able type_able;
		
		//Constructor:
		Distribution_Real(type_stat const & distr_category):
		type_base(distr_category){
		}
		Distribution_Real(type_stat const & distr_category,type_able const & able):
		type_base(distr_category,able){
		}
		Distribution_Real(type_this const & distr):
		type_base((type_base const &)distr){
		}
		//Destructor:
		virtual ~Distribution_Real(void){}
		
		//operator =:
		type_this& operator =(type_this const & distr){
			if(&distr==this) return *this;
			((type_base*)this)->operator =((type_base const &)distr);
			return *this;
		}
		
		//virtual: from type_base:
		virtual type_res rand(type_rand_engine & engine) const override=0;
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const override=0;
		virtual type_real_prob PDF(type_res const & x) const override=0;	
		
		//Cal: calculating the CDF:
		virtual type_real_prob CDF(type_real const & x) const =0;
		//Cal: calculating the CF:
		virtual type_comp CF(type_real const & x) const =0;
		//Cal: calculating the quantile:
		virtual type_res quantile(type_real_prob const & prob) const =0;
}; //class Distribution_Real<RealType>
//}(Distribution of the 1-dim real number)
/****************************************************************************************************/
//Distribution of the real col-vector{
#ifdef LZ_DEF_extLIB_Eigen
template<typename RealType=Type_Real,
	typename ProbRealType=
	typename std::conditional<std::is_floating_point<RealType>::value,RealType,Type_Real>::type,
	typename RandUIntType=typename Type_RandEngine::result_type>
class Distribution_Real_Vec: public Distribution<Type_MatTemp<RealType>,ProbRealType,RandUIntType>{
	LZ_DEF_func_check_traits(std::is_arithmetic<RealType>::value);
	public:
		//Type: this class:
		typedef Distribution_Real_Vec<RealType,ProbRealType,RandUIntType> type_this;
		//Type: type of the base class:
		typedef
			Distribution<Type_MatTemp<RealType>,ProbRealType,RandUIntType>
			type_base;
		//Type: type of the sample:
		typedef typename type_base::type_res type_res;
		//Type: type of the real number for probability:
		typedef typename type_base::type_real_prob type_real_prob;
		//Type: type of the real number:
		typedef
			typename std::conditional<std::is_floating_point<RealType>::value,RealType,type_real_prob>::type
			type_real;
		//Type: type of the complex number:
		typedef Type_CompTemp<type_real> type_comp;
		//Type: StatType:
		typedef typename type_base::type_stat type_stat;
		//Type: type of the sample sequence:
		typedef typename type_base::type_res_seq type_res_seq;
		//Type: type of the std random engine:
		typedef typename type_base::type_rand_engine type_rand_engine;
		//Type: the abilities of this class:
		typedef typename type_base::type_able type_able;
		
		//Constructor:
		Distribution_Real_Vec
		(Type_UInt const & dimension,type_stat const & distr_category):
		type_base(distr_category),dim(dimension>(Type_UInt)0 ? dimension : (Type_UInt)1){
		}
		Distribution_Real_Vec
		(Type_UInt const & dimension,type_stat const & distr_category,type_able const & able):
		type_base(distr_category,able),
		dim(dimension>(Type_UInt)0 ? dimension : (Type_UInt)1){
		}
		Distribution_Real_Vec(type_this const & distr):
		type_base((type_base const &)distr),dim(distr.dim){
		}
		//Destructor:
		virtual ~Distribution_Real_Vec(void){}
		
		//operator =:
		type_this& operator =(type_this const & distr){
			if(&distr==this) return *this;
			((type_base*)this)->operator =((type_base const &)distr);
			dim=distr.dim;
			return *this;
		}
		
		//virtual: from type_base:
		virtual type_res rand(type_rand_engine & engine) const override=0;
		virtual type_res_seq sample(Type_UInt const num,type_rand_engine & engine) const override=0;
		virtual type_real_prob PDF(type_res const & x) const override=0;
		
		//Is: a matrix of size dim*1:
		Type_Bool Is_type_res(type_res const & x) const {
			return x.rows()==this->dim && x.cols()==(Type_UInt)1;
		}
		
		//Cal: calculating the CDF:
		virtual type_real_prob CDF(Type_MatTemp<type_real> const & x) const =0;
		//Cal: calculating the CF:
		virtual type_comp CF(Type_MatTemp<type_real> const & x) const =0;
		
		//Get: dimension of the sample:
		Type_UInt Get_dim(void) const {
			return dim;
		}
	protected:
		Type_UInt dim;
}; //class Distribution_Real_Vec<RealType>
#endif //#ifdef LZ_DEF_extLIB_Eigen
//}(Distribution of the real col-vector)
/****************************************************************************************************/

//see LZ_H/math/src/distribution/distribution_distr.hpp for the derived distributions.

} //namespace math;
} //namespace Liuze;

#include"LZ_H/math/src/distribution/distribution_distr.hpp"

#endif //#ifndef LZ_DEF_LZ_H_math_distribution