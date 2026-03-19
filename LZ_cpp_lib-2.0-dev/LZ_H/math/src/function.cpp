#if LZ_DEF_LZ_H_math_function!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_function_CPP
#define LZ_DEF_LZ_H_math_src_function_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//general arithmetic{

namespace detail{
	
	//Get zero:
	template<typename ValType>
	class FunTemp_zero{
		public:
			static inline ValType zero(){
				return ValType();
			}
			static inline ValType zero(ValType const &){
				return ValType(0);
			}
	}; //class FunTemp_zero;
	
	//Get one:
	template<typename ValType>
	class FunTemp_one{
		public:
			static inline ValType one(){
				return ValType(1);
			}
			static inline ValType one(ValType const &){
				return ValType(1);
			}
	}; //class FunTemp_one;
	
	//Check whether an object is (approximately) zero:
	template<typename ValType>
	class FunTemp_Is_approx_zero{
		public:
			//for non-floating-point ValType:
			static inline
			Type_Bool Is_approx_zero_alt(ValType const & x,void*){
				return x==zero(x);
			}
			//for floating-point ValType:
			static inline
			Type_Bool Is_approx_zero_alt(ValType const & x,void**){
				return abs(x)<=ValType(LZ_DEF_const_default_precision);
			}
			//for general ValType:
			static inline
			Type_Bool Is_approx_zero(ValType const & x){
				return FunTemp_Is_approx_zero<ValType>::Is_approx_zero_alt
					(x,
					static_cast<typename std::conditional<std::is_floating_point<ValType>::value,
					void*,void>::type *>(NULL));
			}
	}; //class FunTemp_Is_approx_zero;
	
	//Inverse for multiplication:
	template<typename ValType>
	class FunTemp_inverse{
		public:
			static inline ValType inverse(ValType const & x){
				return one(x)/x;
			}
	}; //class FunTemp_inverse;
	
	//Take integer power:
	template<typename ValType,typename IntType>
	class FunTemp_powint_default{
		public:
			//for non-invertible ValType:
			static ValType powint_alt(ValType const & base,IntType const & exp,void*){
				IntType ci0=zero(exp);
				if(exp<ci0) return zero(base);
				if(exp==ci0) return one(base);
				if(exp==one(exp)) return base;
				ValType res=base,b=base;
				IntType e=exp-IntType(1);
				while(e!=ci0){
					while(e%IntType(2)==ci0){
						e/=IntType(2);
						b*=b;
					}
					e-=IntType(1);
					res*=b;
				}
				return res;
			}
			//for invertible ValType:
			static ValType powint_alt(ValType const & base,IntType const & exp,void**){
				IntType ci0=zero(exp);
				if(exp==ci0) return one(base);
				if(exp==one(exp)) return base;
				ValType res=one(base),b=base;
				IntType e=exp;
				if(e<ci0){
					b=inverse(b);
					e=-e;
				}
				while(e!=ci0){
					while(e%IntType(2)==ci0){
						e/=IntType(2);
						b*=b;
					}
					e-=IntType(1);
					res*=b;
				}
				return res;
			}
	}; //class FunTemp_powint_default;
	template<typename ValType,typename IntType,
		typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
	class FunTemp_powint{
		public:
			//for general ValType:
			static inline ValType powint(ValType const & base,IntType const & exp){
				return FunTemp_powint_default<ValType,IntType>::powint_alt
					(base,exp,
					static_cast<typename std::conditional<std::is_floating_point<ValType>::value,
					void*,void>::type *>(NULL));
			}
	}; //class FunTemp_powint;
	
} //namespace detail;

#ifdef LZ_DEF_extLIB_Eigen
namespace detail{
	
	//Get zero:
	template<typename ValType>
	class FunTemp_zero<Type_MatTemp<ValType> >{
		public:
			typedef Type_MatTemp<ValType> type_mat;
			static inline type_mat zero(){
				return type_mat(0,0);
			}
			static inline type_mat zero(type_mat const & X){
				return type_mat::Zero(X.rows(),X.cols());
			}
	}; //class FunTemp_zero;
	
	//Get one:
	template<typename ValType>
	class FunTemp_one<Type_MatTemp<ValType> >{
		public:
			typedef Type_MatTemp<ValType> type_mat;
			static inline type_mat one(){
				return type_mat(0,0);
			}
			static inline type_mat one(type_mat const & X){
				return type_mat::Identity(X.rows(),X.cols());
			}
	}; //class FunTemp_one;
	
	//Check whether an object is (approximately) zero:
	template<typename ValType>
	class FunTemp_Is_approx_zero<Type_MatTemp<ValType> >{
		public:
			static Type_Bool Is_approx_zero(Type_MatTemp<ValType> const & X){
				//`return X.size()==Eigen::Index(0) ? true : Type_Bool(X.isApproxToConstant(zero(X(0,0))));
				if(X.size()==Eigen::Index(0)) return true;
				Eigen::Index ir,ic;
				for(ic=0;ic<X.cols();++ic){
					for(ir=0;ir<X.rows();++ir){
						if(!::Liuze::math::Is_approx_zero(X(ir,ic))) return false;
					}
				}
				return true;
			}
	}; //class FunTemp_Is_approx_zero;
	
	//Inverse for multiplication:
	template<typename ValType>
	class FunTemp_inverse<Type_MatTemp<ValType> >{
		public:
			static inline Type_MatTemp<ValType> inverse(Type_MatTemp<ValType> const & X){
				if(X.size()==Eigen::Index(0) || X.rows()!=X.cols()){
					return Type_MatTemp<ValType>(0,0);
				}
				return X.inverse();
			}
	}; //class FunTemp_inverse;
	
	//Take integer power:
	template<typename ValType,typename IntType>
	class FunTemp_powint<Type_MatTemp<ValType>,IntType>{
		public:
			static inline Type_MatTemp<ValType> powint
			(Type_MatTemp<ValType> const & base,IntType const & exp){
				return FunTemp_powint_default<Type_MatTemp<ValType>,IntType>::powint_alt
					(base,exp,
					static_cast<void*>(NULL));
			}
	}; //class FunTemp_powint;
	
} //namespace detail;
#endif //#ifdef LZ_DEF_extLIB_Eigen

//Get zero:
template<typename ValType>
inline
ValType zero(){
	return detail::FunTemp_zero<ValType>::zero();
}
template<typename ValType>
inline
ValType zero(ValType const & x){
	return detail::FunTemp_zero<ValType>::zero(x);
}

//Get one:
template<typename ValType>
inline
ValType one(){
	return detail::FunTemp_one<ValType>::one();
}
template<typename ValType>
inline
ValType one(ValType const & x){
	return detail::FunTemp_one<ValType>::one(x);
}

//Check whether an object is (approximately) zero:
template<typename ValType>
inline
Type_Bool Is_approx_zero(ValType const & x){
	return detail::FunTemp_Is_approx_zero<ValType>::Is_approx_zero(x);
}

//Inverse for multiplication:
template<typename ValType>
inline
ValType inverse(ValType const & x){
	return detail::FunTemp_inverse<ValType>::inverse(x);
}

//Take integer power:
template<typename ValType,typename IntType>
inline
ValType powint(ValType const & base,IntType const & exp){
	return detail::FunTemp_powint<ValType,IntType>::powint(base,exp);
}

//Division with Remainder:
template<typename FactorType>
Type_Bool div
(FactorType const & numer,FactorType const & denom,FactorType * ptr_quot,FactorType * ptr_rem){
	if(!ptr_quot && !ptr_rem) return false;
	if(denom==zero(denom)) return false;
	if(ptr_quot) *ptr_quot=numer/denom;
	if(ptr_rem) *ptr_rem=numer%denom;
	return true;
}

//Division algorithm (Euclidean algorithm):
//return gcd(fac0,fac1)=coef0*fac0+coef1*fac1:
template<typename FactorType>
FactorType EuclideanDivision
(FactorType const & fac0,FactorType const & fac1,FactorType * coef0,FactorType * coef1){
	FactorType const c0=zero(fac0),c1=one(fac0);
	if(fac1==c0){
		if(coef0!=NULL) *coef0=c1;
		if(coef1!=NULL) *coef1=c0;
		return fac0;
	}
	FactorType coef00=c1,coef01=c0,coef10=c0,coef11=c1,nn0=fac0,nn1=fac1;
	FactorType q,r;
	Liuze::math::div(nn0,nn1,std::addressof(q),std::addressof(r));
	while(r!=c0){
		nn0=nn1;
		nn1=r;
		r=coef01;
		coef01=coef00-q*coef01;
		coef00=r;
		r=coef11;
		coef11=coef10-q*coef11;
		coef10=r;
		Liuze::math::div(nn0,nn1,std::addressof(q),std::addressof(r));
	}
	if(coef0!=NULL) *coef0=coef01;
	if(coef1!=NULL) *coef1=coef11;
	return nn1;
}

//}(general arithmetic)
/****************************************************************************************************/
//functions for the integral number types{

//Take modulo:
template<typename IntType,typename ModIntType>
inline
ModIntType mod(IntType const & number,ModIntType const & modulo){
	typedef IntType type_int;
	typedef ModIntType type_mint;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_integral<type_mint>::value));
	typedef typename std::conditional<
		std::is_signed<type_int>::value || std::is_signed<type_mint>::value,
		typename std::common_type<typename std::make_signed<type_int>::type,
			typename std::make_signed<type_mint>::type>::type,
		typename std::common_type<type_int,type_mint>::type>::type
		type_cint; //the common integral type;
	if(number==(type_int)0 || modulo==(type_mint)0) return (type_mint)number;
	type_cint res=static_cast<type_cint>(number)%static_cast<type_cint>(modulo);
	if(res==(type_cint)0 || number>(type_int)0) return type_mint(res);
	return static_cast<type_mint>(abs(static_cast<type_cint>(modulo))+res);
}
//Take division for integral:
template<typename IntType,typename ModIntType>
inline
IntType divint(IntType const & number,ModIntType const & modulo){
	typedef IntType type_int;
	typedef ModIntType type_mint;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_integral<type_mint>::value));
	typedef typename std::conditional<
		std::is_signed<type_int>::value || std::is_signed<type_mint>::value,
		typename std::common_type<typename std::make_signed<type_int>::type,
			typename std::make_signed<type_mint>::type>::type,
		typename std::common_type<type_int,type_mint>::type>::type
		type_cint; //the common integral type;
	if(number==(type_int)0 || modulo==(type_mint)0) return (type_mint)number;
	type_cint res=static_cast<type_cint>(number)/static_cast<type_cint>(modulo);
	if(static_cast<type_cint>(number)%static_cast<type_cint>(modulo)==(type_cint)0 ||
		number>(type_int)0){
		return static_cast<type_mint>(res);
	}
	return static_cast<type_mint>(res+(modulo>(type_mint)0 ? type_cint(-1) : type_cint(1)));
}

//Is co-prime:
template<typename IntType,typename>
Type_Bool Is_coprime(IntType const & n0,IntType const & n1){
	IntType const c0=zero(n0),c1=one(n0);
	IntType x0=abs(n1),x1=abs(n0);
	if(x0==c1 || x1==c1) return (Type_Bool)true;
	if(x0==c0 || x1==c0) return (Type_Bool)false;
	IntType q=x0/x1;
	IntType r=x0-q*x1;
	while(r!=c0){
		x0=x1;
		x1=r;
		q=x0/x1;
		r=x0-q*x1;
	}
	return Type_Bool(x1==c1);
}

//}(functions for the integral number types)
/****************************************************************************************************/
//functions for matrices{
#ifdef LZ_DEF_extLIB_Eigen

//Vectorisation for a matrix:
template<typename EigenDerived>
Type_MatTemp<typename Eigen::DenseBase<EigenDerived>::Scalar>
Vec
(Eigen::DenseBase<EigenDerived> const & X){
	Type_MatTemp<typename Eigen::DenseBase<EigenDerived>::Scalar> vec(X.size(),1);
	for(Type_Size i=0;i<X.cols();++i) vec.block(i*X.rows(),0,X.rows(),1)=X.col(i);
	return vec;
}

//Converting the `Scalar' type of a matrix:
template<typename NewScalar,typename EigenDerived>
Type_MatTemp<NewScalar>
Matrix_Scalar_cast
(Eigen::DenseBase<EigenDerived> const & X){
	typedef typename Eigen::DenseBase<EigenDerived>::Scalar OldScalar;
	typedef typename Eigen::DenseBase<EigenDerived>::Index type_size;
	type_size ir,ic;
	Type_MatTemp<NewScalar> res(X.rows(),X.cols());
	for(ic=0;ic<X.cols();++ic){
		for(ir=0;ir<X.rows();++ir){
			res(ir,ic)=static_cast<NewScalar>((OldScalar)(X(ir,ic)));
		}
	}
	return res;
}

#endif //#ifdef LZ_DEF_extLIB_Eigen

//}(functions for matrices)
/****************************************************************************************************/
/****************************************************************************************************/
//functions to calculate all kinds of combination numbers{
namespace Comb_num{
	
	//the combination number:
	namespace detail{
		//ResType is an integral type:
		template<typename ResType=Type_Real,typename UIntType=Type_UInt,
			typename=typename std::enable_if<
			std::is_arithmetic<ResType>::value && std::is_integral<UIntType>::value
			>::type>
		ResType comb
		(UIntType const & num_total,UIntType const & num_select,void*){
			if(num_select<UIntType(0) || num_select>num_total) return ResType(0);
			if(num_select==UIntType(0) || num_select==num_total) return ResType(1);
			UIntType n0=(num_select<=num_total/UIntType(2) ? num_select : num_total-num_select);
			ResType res_U=ResType(1),res_L=ResType(1),n1=ResType(n0),m1=ResType(num_total);
			while(n1>ResType(0)){
				res_U*=m1;
				res_L*=n1;
				--m1;
				--n1;
			}
			return res_U/res_L;
		}
		//ResType is not an integral type:
		template<typename ResType=Type_Real,typename UIntType=Type_UInt,
			typename=typename std::enable_if<
			std::is_arithmetic<ResType>::value && std::is_integral<UIntType>::value
			>::type>
		ResType comb
		(UIntType const & num_total,UIntType const & num_select,void**){
			if(num_select<UIntType(0) || num_select>num_total) return ResType(0);
			if(num_select==UIntType(0) || num_select==num_total) return ResType(1);
			UIntType n0=(num_select<=num_total/UIntType(2) ? num_select : num_total-num_select);
			ResType res=ResType(1),n1=ResType(n0),m1=ResType(num_total);
			while(n1>ResType(0)){
				res*=m1/n1;
				--m1;
				--n1;
			}
			return res;
		}
	} //namespace detail;
	template<typename ResType,typename UIntType>
	ResType comb(UIntType const & num_total,UIntType const & num_select){
		return detail::comb<ResType,UIntType>
			(num_total,num_select,
			static_cast<typename std::conditional<std::is_integral<ResType>::value,void,void*>::type *>
			(NULL));
	}
	
	//the permutation or arrange number:
	template<typename ResType,typename UIntType,typename>
	ResType perm(UIntType const & num_total,UIntType const & num_select){
		if(num_select<UIntType(0) || num_select>num_total) return ResType(0);
		UIntType i;
		ResType res=1,m=num_total;
		for(i=UIntType(0);i<num_select;++i){
			res*=m;
			--m;
		}
		return res;
	}
	
} //namespace Comb_num;
//}(functions to calculate all kinds of combination numbers)
/****************************************************************************************************/
/****************************************************************************************************/
//class Comb_enum{
//public:
	//Constructor:
	
	template<typename UIntType,typename TraitsCheck>
	Comb_enum<UIntType,TraitsCheck>::Comb_enum
	():
	val(0),stat(type_this::val_stat_null){
	}
	
	template<typename UIntType,typename TraitsCheck>
	Comb_enum<UIntType,TraitsCheck>::Comb_enum
	(type_this const & comb_enum):
	stat(comb_enum.stat),val(comb_enum.val){
	}
	
	template<typename UIntType,typename TraitsCheck>
	Comb_enum<UIntType,TraitsCheck>::Comb_enum
	(type_this && comb_enum):
	stat(std::move(comb_enum.stat)),val(std::move(comb_enum.val)){
	}
	
	//operator:
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_this &
	Comb_enum<UIntType,TraitsCheck>::operator =
	(type_this const & comb_enum){
		if(&comb_enum==this) return *this;
		stat=comb_enum.stat;
		val=comb_enum.val;
		return *this;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_this &
	Comb_enum<UIntType,TraitsCheck>::operator =
	(type_this && comb_enum){
		stat=std::move(comb_enum.stat);
		val=std::move(comb_enum.val);
		return *this;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Type_Bool
	Comb_enum<UIntType,TraitsCheck>::operator ==
	(type_this const & comb_enum) const {
		return &comb_enum==this ? true :
			(stat!=comb_enum.stat ? false :
			(stat==type_this::val_stat_null ? true : (val==comb_enum.val)));
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Comb_enum<UIntType,TraitsCheck>::operator Type_Bool
	() const {
		return this->Is_stat_deref();
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::value_type &
	Comb_enum<UIntType,TraitsCheck>::operator *
	() const {
		return val;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::pointer
	Comb_enum<UIntType,TraitsCheck>::operator ->
	() const {
		return std::addressof(val);
	}
	
	//Is:
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Type_Bool
	Comb_enum<UIntType,TraitsCheck>::Is_stat_null
	() const {
		return stat==type_this::val_stat_null;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Type_Bool
	Comb_enum<UIntType,TraitsCheck>::Is_stat_deref
	() const {
		return stat==type_this::val_stat_deref;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Type_Bool
	Comb_enum<UIntType,TraitsCheck>::Is_stat_end
	() const {
		return stat==type_this::val_stat_end;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	Type_Bool
	Comb_enum<UIntType,TraitsCheck>::Is_stat_rend
	() const {
		return stat==type_this::val_stat_rend;
	}
	
	//Set: the state:
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_stat
	Comb_enum<UIntType,TraitsCheck>::Set_stat_null
	(){
		stat=type_this::val_stat_null;
		return stat;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_stat
	Comb_enum<UIntType,TraitsCheck>::Set_stat_deref
	(){
		stat=type_this::val_stat_deref;
		return stat;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_stat
	Comb_enum<UIntType,TraitsCheck>::Set_stat_end
	(){
		stat=type_this::val_stat_end;
		return stat;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_stat
	Comb_enum<UIntType,TraitsCheck>::Set_stat_rend
	(){
		stat=type_this::val_stat_rend;
		return stat;
	}
	
	//Get:
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_stat
	Comb_enum<UIntType,TraitsCheck>::Get_stat
	() const {
		return stat;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_val
	Comb_enum<UIntType,TraitsCheck>::Get_val
	() const {
		return val;
	}
	
	template<typename UIntType,typename TraitsCheck>
	inline
	typename Comb_enum<UIntType,TraitsCheck>::type_val
	Comb_enum<UIntType,TraitsCheck>::Get_val_next
	(){
		return *(++(*this));
	}
	
	template<typename UIntType,typename TraitsCheck>
	typename Comb_enum<UIntType,TraitsCheck>::type_val
	Comb_enum<UIntType,TraitsCheck>::Get_val_next
	(type_nat const & num_step){
		if(num_step==(type_nat)0 || stat!=type_this::val_stat_deref) return val;
		type_nat i;
		if(num_step>(type_nat)0){
			for(i=(type_nat)1;i<num_step && stat==type_this::val_stat_deref;i++) ++(*this);
			return *(++(*this));
		} else{
			return val; //to be extended;
		}
	}
	
//}(class Comb_enum)
/****************************************************************************************************/
//class Comb_enum_pow{
//public:
	//Constructor:
	
	template<typename UIntType>
	Comb_enum_pow<UIntType>::Comb_enum_pow
	():
	type_base(),base(0),power(0){
	}
	
	template<typename UIntType>
	Comb_enum_pow<UIntType>::Comb_enum_pow
	(type_nat const & num_base,type_nat const & num_power):
	type_base(),
	base(num_base>=(type_nat)0 ? num_base : (type_nat)0),
	power(num_power>=(type_nat)0 ? num_power : (type_nat)0){
		if(power==0){
			this->stat=this->val_stat_deref;
		} else if(base==0){
			this->Set_stat_null();
		} else {
			this->stat=this->val_stat_deref;
			this->val.assign(power,(type_nat)0);
		}
	}
	
	template<typename UIntType>
	Comb_enum_pow<UIntType>::Comb_enum_pow
	(Comb_enum_pow const & comb_enum_pow):
	type_base((type_base const &)comb_enum_pow),base(comb_enum_pow.base),power(comb_enum_pow.power){
	}
	
	template<typename UIntType>
	Comb_enum_pow<UIntType>::Comb_enum_pow
	(Comb_enum_pow && comb_enum_pow):
	base(std::move(comb_enum_pow.base)),power(std::move(comb_enum_pow.power)),
	type_base(static_cast<type_base &&>(comb_enum_pow)){
	}
	
	//operator and virtual: from type_base:
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_this&
	Comb_enum_pow<UIntType>::operator =
	(type_this const & comb_enum_pow){
		if(&comb_enum_pow==this) return *this;
		((type_base*)this)->operator =((type_base const &)comb_enum_pow);
		base=comb_enum_pow.base;
		power=comb_enum_pow.power;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_this&
	Comb_enum_pow<UIntType>::operator =
	(type_this && comb_enum_pow){
		base=std::move(comb_enum_pow.base);
		power=std::move(comb_enum_pow.power);
		((type_base*)this)->operator =(static_cast<type_base &&>(comb_enum_pow));
		return *this;
	}
	
	template<typename UIntType>
	inline
	Type_Bool
	Comb_enum_pow<UIntType>::operator ==
	(type_this const & comb_enum_pow) const {
		return (*((type_base*)this))!=(type_base const &)comb_enum_pow ? false :
			(base==comb_enum_pow.base && power==comb_enum_pow.power);
	}
	
	template<typename UIntType>
	typename Comb_enum_pow<UIntType>::type_this&
	Comb_enum_pow<UIntType>::operator ++
	(){
		if(this->Is_stat_null() || this->Is_stat_end()){
			return *this;
		}
		if(power==(type_nat)0 || base<=(type_nat)1){
			this->Set_stat_end();
			return *this;
		}
		type_nat i_p=power-(type_nat)1;
		while(this->val[i_p]+(type_nat)1==base){
			if(i_p==(type_nat)0){
				this->Set_stat_end();
				return *this;
			}
			--i_p;
		}
		++(this->val[i_p]);
		while(++i_p<power) this->val[i_p]=(type_nat)0;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_stat
	Comb_enum_pow<UIntType>::Set_begin
	(){
		if(this->Is_stat_null()) return this->stat;
		if(power>(type_nat)0 && base>(type_nat)1) this->val.assign(power,(type_nat)0);
		this->stat=this->val_stat_deref;
		return this->stat;
	}
	
	//virtual: deleted in type_base because of the return type:
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::iterator
	Comb_enum_pow<UIntType>::begin
	() const {
		return iterator(base,power);
	}
	
	template<typename UIntType>
	typename Comb_enum_pow<UIntType>::iterator
	Comb_enum_pow<UIntType>::end
	() const {
		iterator res;
		res.base=base;
		res.power=power;
		if(power==0){
			res.Set_stat_end();
		} else if(base==0){
			res.Set_stat_null();
		} else {
			res.Set_stat_end();
			res.val.assign(power,base-(type_nat)1);
		}
		return res;
	}
	
	//other functions:
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_stat
	Comb_enum_pow<UIntType>::Set_begin
	(type_nat const & num_base,type_nat const & num_power){
		base=num_base>=(type_nat)0 ? num_base : (type_nat)0;
		power=num_power>=(type_nat)0 ? num_power : (type_nat)0;
		if(power==0){
			this->stat=this->val_stat_deref;
			this->val.resize((Type_Size)0);
		} else if(base==0){
			this->Set_stat_null();
			this->val.resize((Type_Size)0);
		} else {
			this->stat=this->val_stat_deref;
			this->val.assign(power,(type_nat)0);
		}
		return this->stat;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_stat
	Comb_enum_pow<UIntType>::Set_begin
	(type_this const & comb_enum_pow){
		return this->Set_begin(comb_enum_pow.base,comb_enum_pow.power);
	}
	
	//Get: the parameters:
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_nat
	Comb_enum_pow<UIntType>::Get_base
	() const {
		return base;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_pow<UIntType>::type_nat
	Comb_enum_pow<UIntType>::Get_power
	() const {
		return power;
	}
	
	//static:
	
	template<typename UIntType>
	template<typename IndexType,typename>
	typename Comb_enum_pow<UIntType>::type_stat
	Comb_enum_pow<UIntType>::value_to_index
	(type_nat const & num_base,type_val const & val_pow,IndexType & index){
		typedef IndexType type_index;
		index=0;
		if(num_base<type_nat(0)) return type_this::val_stat_null;
		type_nat num_power=val_pow.size();
		if(num_power==type_nat(0)) return type_this::val_stat_deref;
		if(num_base==type_nat(0)) return type_this::val_stat_null;
		type_index index_wt=1,index_wt_unit=num_base;
		type_nat i=num_power;
		while(i>type_nat(0)){
			--i;
			if(val_pow[i]<type_nat(0) || val_pow[i]>=num_base){
				index=0;
				return type_this::val_stat_null;
			}
			index+=type_index(val_pow[i])*index_wt;
			index_wt*=index_wt_unit;
		}
		return type_this::val_stat_deref;
	}
	
	template<typename UIntType>
	typename Comb_enum_pow<UIntType>::type_stat
	Comb_enum_pow<UIntType>::index_to_value
	(type_nat const & num_base,type_nat const & num_power,type_nat const & index,
	type_val & val_pow){
		if(num_base<type_nat(0) || num_power<type_nat(0) || index<type_nat(0)){
			val_pow.resize(0);
			return type_this::val_stat_null;
		}
		if(num_power==type_nat(0)){
			val_pow.resize(0);
			return type_this::val_stat_deref;
		}
		if(num_base==type_nat(0)){
			val_pow.resize(0);
			return type_this::val_stat_null;
		}
		val_pow.resize(num_power,type_nat(0));
		type_nat index_rem=index;
		type_nat i_pow=num_power;
		while(i_pow>type_nat(0) && index_rem!=type_nat(0)){
			--i_pow;
			val_pow[i_pow]=index_rem%num_base;
			index_rem-=val_pow[i_pow];
			index_rem/=num_base;
		}
		return type_this::val_stat_deref;
	}
	
//}(class Comb_enum_pow)
/****************************************************************************************************/
//class Comb_enum_cart{
//public:
	//Constructor:
	
	template<typename UIntType>
	Comb_enum_cart<UIntType>::Comb_enum_cart
	():
	type_base(),base(){
	}
	
	template<typename UIntType>
	Comb_enum_cart<UIntType>::Comb_enum_cart
	(type_nat_seq const & num_base):
	type_base(),base(num_base){
		this->stat=this->val_stat_deref;
		for(Type_Size i=0;i<base.size();++i){
			if(base[i]<=(type_nat)0){
				base[i]=(type_nat)0;
				this->stat=this->val_stat_null;
			}
		}
		if(this->Is_stat_deref()) this->val.assign(base.size(),(type_nat)0);
	}
	
	template<typename UIntType>
	Comb_enum_cart<UIntType>::Comb_enum_cart
	(Comb_enum_cart const & comb_enum_cart):
	type_base((type_base const &)comb_enum_cart),base(comb_enum_cart.base){
	}
	
	template<typename UIntType>
	Comb_enum_cart<UIntType>::Comb_enum_cart
	(Comb_enum_cart && comb_enum_cart):
	base(std::move(comb_enum_cart.base)),
	type_base(static_cast<type_base &&>(comb_enum_cart)){
	}
	
	//operator and virtual: from type_base:
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::type_this&
	Comb_enum_cart<UIntType>::operator =
	(type_this const & comb_enum_cart){
		if(&comb_enum_cart==this) return *this;
		((type_base*)this)->operator =((type_base const &)comb_enum_cart);
		base=comb_enum_cart.base;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::type_this&
	Comb_enum_cart<UIntType>::operator =
	(type_this && comb_enum_cart){
		base=std::move(comb_enum_cart.base);
		((type_base*)this)->operator =(static_cast<type_base &&>(comb_enum_cart));
		return *this;
	}
	
	template<typename UIntType>
	inline
	Type_Bool
	Comb_enum_cart<UIntType>::operator ==
	(type_this const & comb_enum_cart) const {
		return ((*((type_base*)this))==(type_base const &)comb_enum_cart) &&
			base==comb_enum_cart.base;
	}
	
	template<typename UIntType>
	typename Comb_enum_cart<UIntType>::type_this&
	Comb_enum_cart<UIntType>::operator ++
	(){
		if(this->Is_stat_null() || this->Is_stat_end()){
			return *this;
		}
		if(base.size()==(Type_Size)0){
			this->stat=this->val_stat_end;
			return *this;
		}
		type_nat i_p=base.size()-(type_nat)1;
		while(this->val[i_p]+(type_nat)1==base[i_p]){
			if(i_p==(type_nat)0){
				this->stat=this->val_stat_end;
				return *this;
			}
			--i_p;
		}
		++(this->val[i_p]);
		while((++i_p)<base.size()) this->val[i_p]=(type_nat)0;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::type_stat
	Comb_enum_cart<UIntType>::Set_begin
	(){
		if(this->Is_stat_null()) return this->stat;
		this->val.assign(base.size(),(type_nat)0);
		this->stat=this->val_stat_deref;
		return this->stat;
	}
	
	//virtual: deleted in type_base because of the return type:
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::iterator
	Comb_enum_cart<UIntType>::begin
	() const {
		return iterator(base);
	}
	
	template<typename UIntType>
	typename Comb_enum_cart<UIntType>::iterator
	Comb_enum_cart<UIntType>::end
	() const {
		iterator res;
		res.base=base;
		if(this->Is_stat_null()){
			return res;
		}
		res.stat=this->val_stat_end;
		res.val.resize(base.size());
		for(Type_Size i=0;i<base.size();++i){
			res.val[i]=base[i]-(type_nat)1;
		}
		return res;
	}
	
	//other functions:
	
	template<typename UIntType>
	typename Comb_enum_cart<UIntType>::type_stat
	Comb_enum_cart<UIntType>::Set_begin
	(type_nat_seq const & num_base){
		base=num_base;
		this->stat=this->val_stat_deref;
		for(Type_Size i=0;i<base.size();++i){
			if(base[i]<=(type_nat)0){
				base[i]=(type_nat)0;
				this->stat=this->val_stat_null;
				this->val.resize((Type_Size)0);
			}
		}
		if(this->Is_stat_deref()) this->val.assign(base.size(),(type_nat)0);
		return this->stat;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::type_stat
	Comb_enum_cart<UIntType>::Set_begin
	(type_this const & comb_enum_cart){
		return this->Set_begin(comb_enum_cart.base);
	}
	
	//Get: the parameters:
	
	template<typename UIntType>
	inline
	typename Comb_enum_cart<UIntType>::type_nat
	Comb_enum_cart<UIntType>::Get_base
	() const {
		return base;
	}
	
//}(class Comb_enum_cart)
/****************************************************************************************************/
//class Comb_enum_comb{
//public:
	//Constructor:
	
	template<typename UIntType>
	Comb_enum_comb<UIntType>::Comb_enum_comb
	():
	type_base(),total(0),select(0){
	}
	
	template<typename UIntType>
	Comb_enum_comb<UIntType>::Comb_enum_comb
	(type_nat const & num_total,type_nat const & num_select):
	type_base(),
	total(num_total>=(type_nat)0 ? num_total : (type_nat)0),
	select(num_select>=(type_nat)0 ? num_select : (type_nat)0){
		if(select>total){
			this->stat=this->val_stat_null;
		} else {
			this->stat=this->val_stat_deref;
			this->val.resize(select);
			std::iota(this->val.begin(),this->val.end(),(type_nat)0);
		}
	}
	
	template<typename UIntType>
	Comb_enum_comb<UIntType>::Comb_enum_comb
	(type_this const & comb_enum_comb):
	type_base((type_base const &)comb_enum_comb),
	total(comb_enum_comb.total),select(comb_enum_comb.select){
	}
	
	template<typename UIntType>
	Comb_enum_comb<UIntType>::Comb_enum_comb
	(type_this && comb_enum_comb):
	total(std::move(comb_enum_comb.total)),select(std::move(comb_enum_comb.select)),
	type_base(static_cast<type_base &&>(comb_enum_comb)){
	}
	
	//operator and virtual: from type_base:
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_this&
	Comb_enum_comb<UIntType>::operator =
	(type_this const & comb_enum_comb){
		if(&comb_enum_comb==this) return *this;
		((type_base*)this)->operator =((type_base const &)comb_enum_comb);
		total=comb_enum_comb.total;
		select=comb_enum_comb.select;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_this&
	Comb_enum_comb<UIntType>::operator =
	(type_this && comb_enum_comb){
		total=std::move(comb_enum_comb.total);
		select=std::move(comb_enum_comb.select);
		((type_base*)this)->operator =(static_cast<type_base &&>(comb_enum_comb));
		return *this;
	}
	
	template<typename UIntType>
	inline
	Type_Bool
	Comb_enum_comb<UIntType>::operator ==
	(type_this const & comb_enum_comb) const {
		return (*((type_base*)this))!=(type_base const &)comb_enum_comb ? false :
			(total==comb_enum_comb.total && select==comb_enum_comb.select);
	}
	
	template<typename UIntType>
	typename Comb_enum_comb<UIntType>::type_this&
	Comb_enum_comb<UIntType>::operator ++
	(){
		if(this->Is_stat_null() || this->Is_stat_end()) return *this;
		if(select==(type_nat)0 || select==total){
			this->stat=this->val_stat_end;
			return *this;
		}
		type_nat i=1;
		while(i<=select && this->val[select-i]==total-i) ++i;
		if(i>select){
			this->stat=this->val_stat_end;
			return *this;
		}
		++(this->val[select-i]);
		while((--i)>=1){
			this->val[select-i]=this->val[select-i-1]+1;
		}
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_stat
	Comb_enum_comb<UIntType>::Set_begin
	(){
		if(this->Is_stat_null()) return this->stat;
		if(select>(type_nat)0 && select<total){
			std::iota(this->val.begin(),this->val.end(),(type_nat)0);
		}
		this->stat=this->val_stat_deref;
		return this->stat;
	}
	
	//virtual: deleted in type_base because of the return type:
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::iterator
	Comb_enum_comb<UIntType>::begin
	() const {
		return iterator(total,select);
	}
	
	template<typename UIntType>
	typename Comb_enum_comb<UIntType>::iterator
	Comb_enum_comb<UIntType>::end
	() const {
		iterator res;
		res.total=total;
		res.select=select;
		if(select>total){
			res.stat=this->val_stat_null;
		} else {
			res.stat=this->val_stat_end;
			res.val.resize(select);
			std::iota(res.val.begin(),res.val.end(),total-select);
		}
		return res;
	}
	
	//other functions:
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_stat
	Comb_enum_comb<UIntType>::Set_begin
	(type_nat const & num_total,type_nat const & num_select){
		total=num_total>=(type_nat)0 ? num_total : (type_nat)0;
		select=num_select>=(type_nat)0 ? num_select : (type_nat)0;
		if(select>total){
			this->stat=this->val_stat_null;
			this->val.resize((Type_Size)0);
		} else {
			this->stat=this->val_stat_deref;
			this->val.resize(select);
			std::iota(this->val.begin(),this->val.end(),(type_nat)0);
		}
		return this->stat;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_stat
	Comb_enum_comb<UIntType>::Set_begin
	(type_this const & comb_enum_comb){
		return this->Set_begin(comb_enum_comb.total,comb_enum_comb.select);
	}
	
	//Get: the parameters:
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_nat
	Comb_enum_comb<UIntType>::Get_total
	() const {
		return total;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_comb<UIntType>::type_nat
	Comb_enum_comb<UIntType>::Get_select
	() const {
		return select;
	}
	
	//static:
	
	template<typename UIntType>
	template<typename IndexType,typename>
	typename Comb_enum_comb<UIntType>::type_stat
	Comb_enum_comb<UIntType>::value_to_index
	(type_nat const & num_total,type_val const & val_comb,IndexType & index){
		typedef IndexType type_index;
		typedef typename std::conditional<std::is_integral<type_index>::value,
			Liuze::Type_Real,type_index>::type
			type_index_real;
		type_nat num_select=val_comb.size();
		if(num_total<type_nat(0) || num_select>num_total){
			index=static_cast<type_index>(0);
			return type_this::val_stat_null;
		}
		type_index_real res=0,res_part;
		type_nat i_ele=0,i_pos=0;
		bool tag_i_ele_first;
		type_nat ti0;
		while(i_pos<num_select){
			tag_i_ele_first=true;
			ti0=num_total-num_select+i_pos+type_nat(1);
			while(i_ele<val_comb[i_pos]){
				if(tag_i_ele_first){
					res_part=Liuze::math::Comb_num::comb<type_index_real,type_nat>
						(num_total-i_ele-type_nat(1),num_select-i_pos-type_nat(1));
					res+=res_part;
					tag_i_ele_first=false;
				} else {
					res_part/=(num_total-i_ele);
					res_part*=(ti0-i_ele);
					res+=res_part;
				}
				++i_ele;
				if(i_ele>=num_total){
					index=static_cast<type_index>(0);
					return type_this::val_stat_null;
				}
			}
			++i_ele;
			++i_pos;
		}
		index=static_cast<type_index>(res);
		return type_this::val_stat_deref;
	}
	
	template<typename UIntType>
	typename Comb_enum_comb<UIntType>::type_stat
	Comb_enum_comb<UIntType>::index_to_value
	(type_nat const & num_total,type_nat const & num_select,type_nat const & index,
	type_val & val_comb){
		typedef Liuze::Type_Real type_real;
		if(num_total<type_nat(0) || num_select<type_nat(0) || num_select>num_total ||
		index<type_nat(0)){
			val_comb.resize(0);
			return type_this::val_stat_null;
		}
		val_comb.resize(num_select);
		if(num_select==type_nat(0)){
			return type_this::val_stat_deref;
		}
		if(num_select==num_total){
			std::iota(val_comb.begin(),val_comb.end(),type_nat(0));
			return type_this::val_stat_deref;
		}
		type_nat i_ele=0,i_pos=0;
		type_real id=static_cast<type_real>(index);
		type_real id_part=Liuze::math::Comb_num::comb<type_real>
			(num_total-type_nat(1),num_select-type_nat(1));
		while(i_pos<num_select){
			if(id<id_part){
				val_comb[i_pos]=i_ele;
				++i_ele;
				++i_pos;
				id_part=id_part/(num_total-i_ele)*(num_select-i_pos);
			} else {
				id-=id_part;
				++i_ele;
				if(i_ele>=num_total){
					val_comb.resize(0);
					return type_this::val_stat_null;
				}
				id_part=id_part/(num_total-i_ele)*(num_total-i_ele-num_select+i_pos+type_nat(1));
			}
		}
		return type_this::val_stat_deref;
	}
	
//}(class Comb_enum_comb)
/****************************************************************************************************/
//class Comb_enum_perm{
//public:
	//Constructor:
	
	template<typename UIntType>
	Comb_enum_perm<UIntType>::Comb_enum_perm
	():
	type_base(),order(0){
		this->stat=this->val_stat_deref;
	}
	
	template<typename UIntType>
	Comb_enum_perm<UIntType>::Comb_enum_perm
	(type_nat const & num_order):
	type_base(),order(num_order>=(type_nat)0 ? num_order : (type_nat)0){
		this->stat=this->val_stat_deref;
		this->val.resize(order);
		std::iota(this->val.begin(),this->val.end(),(type_nat)0);
	}
	
	template<typename UIntType>
	Comb_enum_perm<UIntType>::Comb_enum_perm
	(type_this const & comb_enum_perm):
	type_base((type_base const &)comb_enum_perm),order(comb_enum_perm.order){
	}
	
	template<typename UIntType>
	Comb_enum_perm<UIntType>::Comb_enum_perm
	(type_this && comb_enum_perm):
	order(std::move(comb_enum_perm.order)),
	type_base(static_cast<type_base &&>(comb_enum_perm)){
	}
	
	//operator and virtual: from type_base:
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_this&
	Comb_enum_perm<UIntType>::operator =
	(type_this const & comb_enum_perm){
		if(&comb_enum_perm==this) return *this;
		((type_base*)this)->operator =((type_base const &)comb_enum_perm);
		order=comb_enum_perm.order;
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_this&
	Comb_enum_perm<UIntType>::operator =
	(type_this && comb_enum_perm){
		order=std::move(comb_enum_perm.order);
		((type_base*)this)->operator =(static_cast<type_base &&>(comb_enum_perm));
		return *this;
	}
	
	template<typename UIntType>
	inline
	Type_Bool
	Comb_enum_perm<UIntType>::operator ==
	(type_this const & comb_enum_perm) const {
		return (*((type_base*)this))!=(type_base const &)comb_enum_perm ? false :
			(order==comb_enum_perm.order);
	}
	
	template<typename UIntType>
	typename Comb_enum_perm<UIntType>::type_this&
	Comb_enum_perm<UIntType>::operator ++
	(){
		if(this->Is_stat_null() || this->Is_stat_end()) return *this;
		if(order<=(type_nat)1){
			this->stat=this->val_stat_end;
			return *this;
		}
		type_nat i=order-(type_nat)1;
		while(this->val[i]<this->val[--i]){
			if(i==(type_nat)0){
				this->stat=this->val_stat_end;
				return *this;
			}
		}
		std::reverse(this->val.begin()+(i+(type_nat)1),this->val.end());
		type_nat i0=i+(type_nat)1,i1=order-(type_nat)1,im;
		while(i0!=i1){
			im=(i0+i1)/(type_nat)2;
			if(this->val[im]<this->val[i]){
				i0=im+(type_nat)1;
			} else {
				i1=im;
			}
		}
		std::swap(this->val[i],this->val[i0]);
		return *this;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_stat
	Comb_enum_perm<UIntType>::Set_begin
	(){
		this->stat=this->val_stat_deref;
		std::iota(this->val.begin(),this->val.end(),(type_nat)0);
		return this->stat;
	}
	
	//virtual: deleted in type_base because of the return type:
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::iterator
	Comb_enum_perm<UIntType>::begin
	() const {
		return iterator(order);
	}
	
	template<typename UIntType>
	typename Comb_enum_perm<UIntType>::iterator
	Comb_enum_perm<UIntType>::end
	() const {
		iterator res;
		res.order=order;
		res.stat=this->val_stat_end;
		res.val.resize(order);
		type_nat j=order-(type_nat)1;
		for(type_nat i=0;i<order;++i) res.val[i]=j-i;
		return res;
	}
	
	//other functions:
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_stat
	Comb_enum_perm<UIntType>::Set_begin(type_nat const & num_order){
		order=num_order>=(type_nat)0 ? num_order : (type_nat)0;
		this->stat=this->val_stat_deref;
		this->val.resize(order);
		std::iota(this->val.begin(),this->val.end(),(type_nat)0);
		return this->stat;
	}
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_stat
	Comb_enum_perm<UIntType>::Set_begin
	(type_this const & comb_enum_perm){
		return this->Set_begin(comb_enum_perm.order);
	}
	
	//Get: the parameters:
	
	template<typename UIntType>
	inline
	typename Comb_enum_perm<UIntType>::type_nat
	Comb_enum_perm<UIntType>::Get_order
	() const {
		return order;
	}
	
	//static:
	
	template<typename UIntType>
	template<typename IndexType,typename>
	typename Comb_enum_perm<UIntType>::type_stat
	Comb_enum_perm<UIntType>::value_to_index
	(type_val const & val_perm,IndexType & index){
		typedef IndexType type_index;
		type_nat num_order=val_perm.size();
		index=0;
		if(num_order<=type_nat(1)) return type_this::val_stat_deref;
		type_index index_part=1;
		type_nat i_comp,i_comp_rest,j_comp,coef_index_part;
		for(i_comp_rest=2;i_comp_rest<=num_order;++i_comp_rest){
			i_comp=num_order-i_comp_rest;
			coef_index_part=type_nat(0);
			//consider val_perm as the positions of components in a permutation:
			for(j_comp=i_comp+type_nat(1);j_comp<num_order;++j_comp){
				if(val_perm[j_comp]<val_perm[i_comp]) ++coef_index_part;
			}
			index+=index_part*type_index(coef_index_part);
			index_part*=type_index(i_comp_rest);
		}
		return type_this::val_stat_deref;
	}
	
	template<typename UIntType>
	typename Comb_enum_perm<UIntType>::type_stat
	Comb_enum_perm<UIntType>::index_to_value
	(type_nat const & num_order,type_nat const & index,type_val & val_perm){
		if(num_order<type_nat(0) || index<type_nat(0)){
			val_perm.resize(0);
			return type_this::val_stat_null;
		}
		if(num_order<=type_nat(1)){
			val_perm.assign(num_order,type_nat(0));
			return type_this::val_stat_deref;
		}
		val_perm.resize(num_order);
		type_val val_perm_inv(num_order);
		std::vector<bool> tag_perm(num_order,false);
		type_nat i_comp,coef_id_part,i_pos,i_coef=num_order-type_nat(1);
		type_nat id_part=Liuze::math::Comb_num::perm<type_nat>(i_coef,i_coef);
		type_nat id=index;
		for(i_comp=0;i_comp<num_order;++i_comp){
			if(i_comp>type_nat(0)) id_part/=num_order-i_comp;
			coef_id_part=id/id_part;
			id%=id_part;
			i_coef=i_pos=0;
			while(tag_perm[i_pos] || i_coef<coef_id_part){
				while(tag_perm[i_pos]) ++i_pos;
				while(tag_perm[i_pos]==false && i_coef<coef_id_part){
					++i_coef;
					++i_pos;
				}
			}
			val_perm_inv[i_pos]=i_comp;
			tag_perm[i_pos]=true;
		}
		//the inverse of val_perm_inv:
		for(i_comp=0;i_comp<num_order;++i_comp) val_perm[val_perm_inv[i_comp]]=i_comp;
		return type_this::val_stat_deref;
	}
	
//}(class Comb_enum_perm)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_function_CPP
#endif //#if LZ_DEF_LZ_H_math_function!=202105L