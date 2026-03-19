#if LZ_DEF_LZ_H_math_polynomial!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_polynomial_CPP
#define LZ_DEF_LZ_H_math_src_polynomial_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//class Polynomial{
//public:
	//Constructor:
	
	template<typename CoefType>
	Polynomial<CoefType>::Polynomial
	():
	m_arr_coef(0),m_order(0){
	}
	
	template<typename CoefType>
	template<typename InputIterType,typename>
	Polynomial<CoefType>::Polynomial
	(InputIterType iter_coef_0,InputIterType const & iter_coef_end):
	m_arr_coef(iter_coef_0,iter_coef_end),m_order(0){
		type_size i_ord=m_arr_coef.size();
		while(i_ord>type_size(0)){
			--i_ord;
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])){
				m_order=i_ord+type_size(1);
				break;
			}
		}
	}
	
	template<typename CoefType>
	template<typename InputIterType,typename>
	Polynomial<CoefType>::Polynomial
	(InputIterType iter_coef_0,type_size const & length):
	m_arr_coef(length),m_order(0){
		for(type_size i_ord=0;i_ord<length;++i_ord){
			m_arr_coef[i_ord]=(*iter_coef_0);
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])) m_order=i_ord+type_size(1);
			++iter_coef_0;
		}
	}
	
	template<typename CoefType>
	template<typename CoefVecType,typename>
	Polynomial<CoefType>::Polynomial
	(CoefVecType const & coef_vector):
	m_arr_coef(coef_vector.size()),m_order(0){
		for(type_size i_ord=0;i_ord<m_arr_coef.size();++i_ord){
			m_arr_coef[i_ord]=coef_vector[i_ord];
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])) m_order=i_ord+type_size(1);
		}
	}
	
	template<typename CoefType>
	Polynomial<CoefType>::Polynomial
	(type_this const & poly):
	m_arr_coef(poly.m_arr_coef),m_order(poly.m_order){
	}
	
	template<typename CoefType>
	Polynomial<CoefType>::Polynomial
	(type_this && poly):
	m_arr_coef(std::move(poly.m_arr_coef)),m_order(std::move(poly.m_order)){
	}
	
	//operator:
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator =
	(type_this const & poly){
		if(this==std::addressof(poly)) return *this;
		m_arr_coef=poly.m_arr_coef;
		m_order=poly.m_order;
		return *this;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator =
	(type_this && poly){
		m_arr_coef=std::move(poly.m_arr_coef);
		m_order=std::move(poly.m_order);
		return *this;
	}
	
	template<typename CoefType>
	template<typename CoefVecType,typename>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator =
	(CoefVecType const & coef_vector){
		m_order=0;
		m_arr_coef.resize(coef_vector.size());
		for(type_size i_ord=0;i_ord<m_arr_coef.size();++i_ord){
			m_arr_coef[i_ord]=coef_vector[i_ord];
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])) m_order=i_ord+type_size(1);
		}
		return *this;
	}
	
	template<typename CoefType>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator +=
	(type_this const & poly){
		if(poly.m_order==type_size(0)) return *this;
		type_size ord_max=std::max(m_order,poly.m_order);
		m_arr_coef.resize(ord_max,poly.coef_zero());
		type_size i_ord;
		m_order=0;
		for(i_ord=0;i_ord<ord_max;++i_ord){
			m_arr_coef[i_ord]+=poly.coef(i_ord);
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])) m_order=i_ord+type_size(1);
		}
		return *this;
	}
	
	template<typename CoefType>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator -=
	(type_this const & poly){
		if(poly.m_order==type_size(0)) return *this;
		type_size ord_max=std::max(m_order,poly.m_order);
		m_arr_coef.resize(ord_max,poly.coef_zero());
		type_size i_ord;
		m_order=0;
		for(i_ord=0;i_ord<ord_max;++i_ord){
			m_arr_coef[i_ord]-=poly.coef(i_ord);
			if(!Liuze::math::Is_approx_zero(m_arr_coef[i_ord])) m_order=i_ord+type_size(1);
		}
		return *this;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator *=
	(type_this const & poly){
		return (*this=(*this)*poly);
	}
	
	template<typename CoefType>
	template<typename AuxIntType,typename>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator *=
	(AuxIntType const & num_times){
		if(this->Is_zero()) return *this;
		if(num_times==AuxIntType(0)) return *this=this->zero();
		AuxIntType n_times;
		if(num_times>AuxIntType(0)){
			n_times=num_times;
		} else {
			n_times=-num_times;
			*this=-(*this);
		}
		type_size i_deg;
		for(i_deg=0;i_deg<this->order();++i_deg) this->m_arr_coef[i_deg]*=n_times;
		this->m_order=type_size(0);
		while(i_deg>type_size(0)){
			--i_deg;
			if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_deg])){
				this->m_order=i_deg+type_size(1);
				break;
			}
		}
		return *this;
	}
	
	template<typename CoefType>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator *=
	(type_coef const & val_coef){
		if(Liuze::math::Is_approx_zero(val_coef)) return *this=this->zero();
		type_size i_deg;
		for(i_deg=0;i_deg<this->m_order;++i_deg) (this->m_arr_coef[i_deg])*=val_coef;
		i_deg=this->m_order;
		this->m_order=type_size(0);
		while(i_deg>type_size(0)){
			--i_deg;
			if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_deg])){
				this->m_order=i_deg+type_size(1);
				return *this;
			}
		}
		return *this;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator /=
	(type_this const & poly){
		type_this res;
		type_this::div(*this,poly,std::addressof(res),static_cast<type_this *>(NULL));
		return *this=res;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::operator %=
	(type_this const & poly){
		type_this res;
		type_this::div(*this,poly,static_cast<type_this *>(NULL),std::addressof(res));
		return *this=res;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this
	Polynomial<CoefType>::operator +
	() const {
		return *this;
	}
	
	template<typename CoefType>
	typename Polynomial<CoefType>::type_this
	Polynomial<CoefType>::operator -
	() const {
		type_this res;
		res.m_order=this->m_order;
		res.m_arr_coef.resize(res.order());
		for(type_size i_ord=0;i_ord<res.order();++i_ord) res[i_ord]=-(*this)[i_ord];
		return res;
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator +
	(Polynomial<AuxCoefType> poly0,Polynomial<AuxCoefType> const & poly1){
		return (poly0+=poly1);
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator -
	(Polynomial<AuxCoefType> poly0,Polynomial<AuxCoefType> const & poly1){
		return (poly0-=poly1);
	}
	
	template<typename AuxCoefType>
	Polynomial<AuxCoefType> operator *
	(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1){
		typedef Polynomial<AuxCoefType> type_poly;
		typedef typename type_poly::type_coef type_coef;
		typedef typename type_poly::type_size type_size;
		type_poly res;
		if(poly0.order()==type_size(0) || poly1.order()==type_size(0)) return res;
		res.m_arr_coef.assign(poly0.order()+poly1.order()-type_size(1),poly0.coef_zero());
		type_size i0_ord,i1_ord;
		//calculating `res.m_arr_coef':
		for(i1_ord=0;i1_ord<poly1.order();++i1_ord){
			if(Liuze::math::Is_approx_zero(poly1.coef(i1_ord))) continue;
			for(i0_ord=0;i0_ord<poly0.order();++i0_ord){
				res.m_arr_coef[i0_ord+i1_ord]+=poly0.coef(i0_ord)*poly1.coef(i1_ord);
			}
		}
		//calculating `res.m_order'{
		i0_ord=res.m_arr_coef.size();
		while(i0_ord>type_size(0)){
			--i0_ord;
			if(!Liuze::math::Is_approx_zero(res.m_arr_coef[i0_ord])){
				res.m_order=i0_ord+type_size(1);
				break;
			}
		}
		//}
		return res;
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator *
	(typename Polynomial<AuxCoefType>::type_coef const & val_coef,Polynomial<AuxCoefType> poly){
		return poly*=val_coef;
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator *
	(Polynomial<AuxCoefType> poly,typename Polynomial<AuxCoefType>::type_coef const & val_coef){
		return poly*=val_coef;
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator /
	(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1){
		typedef Polynomial<AuxCoefType> type_poly;
		type_poly res;
		type_poly::div(poly0,poly1,std::addressof(res),static_cast<type_poly *>(NULL));
		return res;
	}
	
	template<typename AuxCoefType>
	inline
	Polynomial<AuxCoefType> operator %
	(Polynomial<AuxCoefType> const & poly0,Polynomial<AuxCoefType> const & poly1){
		typedef Polynomial<AuxCoefType> type_poly;
		type_poly res;
		type_poly::div(poly0,poly1,static_cast<type_poly *>(NULL),std::addressof(res));
		return res;
	}
	
	template<typename CoefType>
	Type_Bool
	Polynomial<CoefType>::operator ==
	(type_this const & poly) const {
		if(m_order!=poly.m_order) return false;
		for(type_size i_ord=0;i_ord<m_order;++i_ord){
			if(m_arr_coef[i_ord]!=poly.m_arr_coef[i_ord]) return false;
		}
		return true;
	}
	
	template<typename CoefType>
	inline
	Type_Bool
	Polynomial<CoefType>::operator !=
	(type_this const & poly) const {
		return !(*this==poly);
	}
	
	template<typename CoefType>
	inline
	Type_Bool
	Polynomial<CoefType>::operator !
	() const {
		return this->Is_zero();
	}
	
	template<typename CoefType>
	inline
	Polynomial<CoefType>::operator Type_Bool
	() const {
		return !(this->Is_zero());
	}
	
	template<typename CoefType>
	template<typename ValType>
	ValType
	Polynomial<CoefType>::operator ()
	(ValType const & x) const {
		typedef ValType type_val;
		if(this->order()==type_size(0)) return type_val(this->coef_zero()*x);
		if(this->order()==type_size(1)) return type_val(this->coef(0)*one(x));
		type_val res=this->coef(0)*one(x);
		type_val x_pow=x;
		for(type_size i_deg=1;i_deg<this->order();++i_deg){
			res+=type_val(this->coef(i_deg)*x_pow);
			x_pow*=x;
		}
		return res;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_coef const &
	Polynomial<CoefType>::operator []
	(type_size const & i_deg) const {
		return m_arr_coef[i_deg];
	}
	
	//Is:
	
	template<typename CoefType>
	inline
	Type_Bool
	Polynomial<CoefType>::Is_zero
	() const {
		return this->order()==type_size(0);
	}
	
	//value:
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_size
	Polynomial<CoefType>::order
	() const {
		return m_order;
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_size
	Polynomial<CoefType>::degree
	() const {
		return m_order>type_size(0) ? m_order-type_size(1) : type_size(0);
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_coef
	Polynomial<CoefType>::coef_zero
	() const {
		return m_order>type_size(0) ?
			Liuze::math::zero(m_arr_coef[m_order-type_size(1)]) :
			Liuze::math::zero<type_coef>();
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_coef
	Polynomial<CoefType>::coef
	(type_size const & i_deg) const {
		return i_deg<this->m_order ? this->m_arr_coef[i_deg] : this->coef_zero();
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this
	Polynomial<CoefType>::zero
	() const {
		return type_this();
	}
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_this
	Polynomial<CoefType>::constant
	(type_coef const & coeff_deg0) const {
		type_this res;
		if(Liuze::math::Is_approx_zero(coeff_deg0)) return res;
		res.m_order=type_size(1);
		res.m_arr_coef.resize(res.m_order);
		res.m_arr_coef[0]=coeff_deg0;
		return res;
	}
	
	//other operator:
	
	template<typename CoefType>
	template<typename InputIterType,typename>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::assign
	(InputIterType iter_coef_0,InputIterType const & iter_coef_end){
		this->m_arr_coef.assign(iter_coef_0,iter_coef_end);
		type_size i_ord=this->m_arr_coef.size();
		this->m_order=0;
		while(i_ord>type_size(0)){
			--i_ord;
			if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_ord])){
				this->m_order=i_ord+type_size(1);
				break;
			}
		}
		return *this;
	}
	
	template<typename CoefType>
	template<typename InputIterType,typename>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::assign
	(InputIterType iter_coef_0,type_size const & length){
		this->m_arr_coef.resize(length);
		this->m_order=0;
		for(type_size i_ord=0;i_ord<length;++i_ord){
			this->m_arr_coef[i_ord]=(*iter_coef_0);
			if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_ord])){
				this->m_order=i_ord+type_size(1);
			}
			++iter_coef_0;
		}
		return *this;
	}
	
	template<typename CoefType>
	template<typename CoefVecType,typename>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::assign
	(CoefVecType const & coef_vector){
		this->m_arr_coef.resize(coef_vector.size());
		this->m_order=0;
		for(type_size i_ord=0;i_ord<this->m_arr_coef.size();++i_ord){
			this->m_arr_coef[i_ord]=coef_vector[i_ord];
			if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_ord])){
				this->m_order=i_ord+type_size(1);
			}
		}
		return *this;
	}
	
	template<typename CoefType>
	template<typename AuxIntType,typename>
	typename Polynomial<CoefType>::type_this
	Polynomial<CoefType>::pow
	(AuxIntType const & num_exp) const {
		return Liuze::math::powint<type_this,AuxIntType>(*this,exp);
	}
	
	//Set:
	
	template<typename CoefType>
	typename Polynomial<CoefType>::type_this &
	Polynomial<CoefType>::Set_coef
	(type_size const & i_deg,type_coef const & val_coef){
		if(i_deg+type_size(1)<this->order()){
			this->m_arr_coef[i_deg]=val_coef;
			return *this;
		}
		if(i_deg+type_size(1)==this->order()){
			this->m_arr_coef[i_deg]=val_coef;
			if(Liuze::math::Is_approx_zero(val_coef)){
				this->m_order=0;
				type_size i_ord=i_deg;
				while(i_ord>type_size(0)){
					--i_ord;
					if(!Liuze::math::Is_approx_zero(this->m_arr_coef[i_ord])){
						this->m_order=i_ord+type_size(1);
						break;
					}
				}
			}
			return *this;
		}
		if(Liuze::math::Is_approx_zero(val_coef)) return *this;
		this->m_order=i_deg+type_size(1);
		this->m_arr_coef.resize(this->m_order,Liuze::math::zero(val_coef));
		this->m_arr_coef[i_deg]=val_coef;
		return *this;
	}
	
	//input and output:
	
	template<typename AuxCoefType>
	std::ostream &
	operator <<
	(std::ostream & out,Polynomial<AuxCoefType> const & poly){
		if(poly.order()==Type_Size(0)) return out;
		Type_Size i_deg=0;
		out<<poly.coef(i_deg);
		for(++i_deg;i_deg<poly.order();++i_deg) out<<" "<<poly.coef(i_deg);
		return out;
	}
	
	//static:
	
	template<typename CoefType>
	Type_Stat
	Polynomial<CoefType>::div
	(type_this const & numer,type_this const & denom,type_this * ptr_quot,type_this * ptr_rem){
		if(!ptr_quot && !ptr_rem) return 2; //no need to calculate;
		if(denom.Is_zero()){
			if(ptr_quot) *ptr_quot=denom;
			if(ptr_rem) *ptr_rem=denom;
			return 1;
		}
		if(numer.order()<denom.order()){
			if(ptr_quot) *ptr_quot=denom.zero();
			if(ptr_rem) *ptr_rem=numer;
			return 0;
		}
		Type_Stat val_ret=0;
		Type_Bool flag_new_quot=false,flag_new_rem=false;
		if(!ptr_quot){
			ptr_quot=new type_this;
			if(!ptr_quot){
				val_ret=11;
				goto GTS_end;
			}
			flag_new_quot=true;
		}
		if(!ptr_rem){
			ptr_rem=new type_this;
			if(!ptr_rem){
				val_ret=11;
				goto GTS_end;
			}
			flag_new_rem=true;
		}
		//dividing:
		{
			(*ptr_rem)=numer;
			(*ptr_quot)=denom.zero();
			type_size deg_denom=denom.order()-type_size(1);
			type_size i_deg_quot=ptr_rem->order()-denom.order();
			ptr_quot->m_order=i_deg_quot+type_size(1);
			ptr_quot->m_arr_coef.resize(ptr_quot->m_order,denom.coef_zero());
			type_size i_deg_rem=ptr_rem->order()-type_size(1),i_deg_denom;
			while(ptr_rem->order()>=denom.order()){
				i_deg_quot=ptr_rem->order()-denom.order();
				(*ptr_quot)[i_deg_quot]=(*ptr_rem)[i_deg_rem]/denom[deg_denom];
				for(i_deg_denom=0;i_deg_denom<deg_denom;++i_deg_denom){
					(*ptr_rem)[i_deg_denom+i_deg_quot]-=(*ptr_quot)[i_deg_quot]*denom[i_deg_denom];
				}
				(*ptr_rem)[i_deg_rem]=denom.coef_zero();
				while(i_deg_rem>type_size(0)){
					--i_deg_rem;
					if(!Liuze::math::Is_approx_zero((*ptr_rem)[i_deg_rem])) break;
				}
				if(Liuze::math::Is_approx_zero((*ptr_rem)[i_deg_rem])){
					(*ptr_rem)=denom.zero();
				} else {
					ptr_rem->m_order=i_deg_rem+type_size(1);
				}
			}
		}
		GTS_end:
		if(flag_new_quot) delete ptr_quot;
		else if(flag_new_rem) delete ptr_quot;
		return val_ret;
	}
	
//protected:
	//operator:
	
	template<typename CoefType>
	inline
	typename Polynomial<CoefType>::type_coef &
	Polynomial<CoefType>::operator []
	(type_size const & i_deg){
		return m_arr_coef[i_deg];
	}
	
//}(class Polynomial)
/****************************************************************************************************/
//operators of class Polynomial{

template<typename IntType,typename CoefType,typename>
inline
Polynomial<CoefType> operator *
(IntType const & num_times,Polynomial<CoefType> poly){
	return poly*=num_times;
}
template<typename IntType,typename CoefType,typename>
inline
Polynomial<CoefType> operator *
(Polynomial<CoefType> poly,IntType const & num_times){
	return poly*=num_times;
}

//}(operators of class Polynomial)
/****************************************************************************************************/
//arithmetic functions{

namespace detail{
	
	template<typename CoefType>
	class FunTemp_zero<Polynomial<CoefType> >{
		public:
			typedef Polynomial<CoefType> type_val;
			static inline type_val zero(){
				return type_val();
			}
			static inline type_val zero(type_val const & x){
				return x.zero();
			}
	};
	
	template<typename CoefType>
	class FunTemp_one<Polynomial<CoefType> >{
		public:
			typedef Polynomial<CoefType> type_val;
			static inline type_val one(){
				return FunTemp_one<type_val>::one(type_val());
			}
			static inline type_val one(type_val const & x){
				return x.constant(one(x.coef_zero()));
			}
	};
	
	template<typename CoefType>
	class FunTemp_Is_approx_zero<Polynomial<CoefType> >{
		public:
			typedef Polynomial<CoefType> type_val;
			static inline Type_Bool Is_approx_zero(type_val const & x){
				return x.Is_zero();
			}
	};
	
	template<typename CoefType,typename ExpIntType>
	class FunTemp_powint<Polynomial<CoefType>,ExpIntType>{
		public:
			typedef Polynomial<CoefType> type_val;
			static inline type_val powint(type_val const & base,ExpIntType const & exp){
				return FunTemp_powint_default<type_val,ExpIntType>::powint_alt
					(base,exp,static_cast<void*>(NULL));
			}
	};
	
} //namespace detail;

template<typename CoefType>
inline
Type_Bool div
(Polynomial<CoefType> const & numer,Polynomial<CoefType> const & denom,
Polynomial<CoefType> * ptr_quot,Polynomial<CoefType> * ptr_rem){
	return Polynomial<CoefType>::div(numer,denom,ptr_quot,ptr_rem)==Type_Stat(0);
}

//}(arithmetic functions)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_polynomial_CPP
#endif //#if LZ_DEF_LZ_H_math_polynomial!=202105L