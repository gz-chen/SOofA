#if LZ_DEF_LZ_H_math_intresidue!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_intresidue_CPP
#define LZ_DEF_LZ_H_math_src_intresidue_CPP 202105L

/****************************************************************************************************/
namespace Liuze{
namespace math{
namespace IntegerResidueClass{

/****************************************************************************************************/
//the ring:

/****************************************************************************************************/
//class Ring_Base{
//public:
	//Constructor:
	
	template<typename IntType,typename TraitsCheck>
	Ring_Base<IntType,TraitsCheck>::Ring_Base
	():
	m_modulo(type_this::val_modulo_universal){
	}
	
	template<typename IntType,typename TraitsCheck>
	template<typename AuxModIntType,typename>
	Ring_Base<IntType,TraitsCheck>::Ring_Base
	(AuxModIntType const & num_modulo):
	m_modulo(abs(num_modulo)){
	}
	
	template<typename IntType,typename TraitsCheck>
	Ring_Base<IntType,TraitsCheck>::Ring_Base
	(type_this const & ring):
	m_modulo(ring.m_modulo){
	}
	
	template<typename IntType,typename TraitsCheck>
	Ring_Base<IntType,TraitsCheck>::Ring_Base
	(type_this && ring):
	m_modulo(std::move(ring.m_modulo)){
	}
	
	//operator:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_this &
	Ring_Base<IntType,TraitsCheck>::operator =
	(type_this const & ring){
		if(std::addressof(ring)==this) return *this;
		m_modulo=ring.m_modulo;
		return *this;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_this &
	Ring_Base<IntType,TraitsCheck>::operator =
	(type_this && ring){
		m_modulo=std::move(ring.m_modulo);
		return *this;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::operator ==
	(type_this const & ring) const {
		return (this->Is_invalid() || ring.Is_invalid()) ?
			(this->Is_invalid() && ring.Is_invalid()) :
			((this->Is_universal() || ring.Is_universal()) ?
			true : (this->m_modulo==ring.m_modulo));
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::operator !=
	(type_this const & ring) const {
		return !(*this==ring);
	}
	
	template<typename IntType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::operator ()
	(AuxIntType const & num_value) const {
		return this->element(num_value);
	}
	
	//Is:
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::Is_valid
	() const {
		return m_modulo!=type_this::val_modulo_invalid;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::Is_invalid
	() const {
		return m_modulo==type_this::val_modulo_invalid;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::Is_universal
	() const {
		return m_modulo==type_this::val_modulo_universal;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::Is_specific
	() const {
		return (this->Is_valid() && !(this->Is_universal()));
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	Type_Bool
	Ring_Base<IntType,TraitsCheck>::Is_finite
	() const {
		return (this->Is_specific() && this->m_modulo!=type_mod(0));
	}
	
	//value:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_val
	Ring_Base<IntType,TraitsCheck>::size
	() const {
		return this->Is_finite() ? type_val(this->m_modulo) : type_val(0);
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_mod
	Ring_Base<IntType,TraitsCheck>::modulo
	() const {
		return this->m_modulo;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_space_act const &
	Ring_Base<IntType,TraitsCheck>::space_act
	() const {
		return *this;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_space_act &
	Ring_Base<IntType,TraitsCheck>::space_act
	(){
		return *this;
	}
	
	//other function:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element_zero
	() const {
		return type_element(Liuze::math::zero<type_val>(),*this);
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element_one
	() const {
		return type_element(Liuze::math::one<type_val>(),*this);
	}
	
	template<typename IntType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element
	(AuxIntType const & num_value) const {
		return type_element(type_val(num_value),*this);
	}
	
	//static:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_this
	Ring_Base<IntType,TraitsCheck>::invalid
	(){
		type_this res;
		res.m_modulo=type_this::val_modulo_invalid;
		return res;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_this
	Ring_Base<IntType,TraitsCheck>::universal
	(){
		return type_this();
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element_invalid
	(){
		return type_element::invalid();
	}
	
	template<typename IntType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element_invalid
	(AuxIntType const & num_value){
		return type_element::invalid(num_value);
	}
	
	template<typename IntType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Ring_Base<IntType,TraitsCheck>::type_element
	Ring_Base<IntType,TraitsCheck>::element_universal
	(AuxIntType const & num_value){
		return type_element::universal(num_value);
	}
	
//}(class Ring_Base)
/****************************************************************************************************/

//(the ring)
/****************************************************************************************************/
/****************************************************************************************************/
//the element:

/****************************************************************************************************/
//class Element{
//public:
	//Constructor:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	Element<IntType,RingType,TraitsCheck>::Element
	():
	m_ring(type_space::universal()),m_value(type_val(0)){
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,RingType,TraitsCheck>::Element
	(AuxIntType const & num_value):
	m_ring(type_space::universal()),m_value(type_val(num_value)){
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,RingType,TraitsCheck>::Element
	(AuxIntType const & num_value,type_space const & ring):
	m_ring(ring),
	m_value(ring.Is_finite() ? type_val(mod(num_value,ring.modulo())) : type_val(num_value)){
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,RingType,TraitsCheck>::Element
	(AuxIntType const & num_value,type_this const & element):
	m_ring(element.m_ring),
	m_value(element.Is_invalid() ? element.m_value :
	(element.space_act().Is_finite() ?
	type_val(mod(num_value,element.space_act().modulo())) : type_val(num_value))){
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	Element<IntType,RingType,TraitsCheck>::Element
	(type_this const & element):
	m_ring(element.m_ring),m_value(element.m_value){
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	Element<IntType,RingType,TraitsCheck>::Element
	(type_this && element):
	m_ring(std::move(element.m_ring)),m_value(std::move(element.m_value)){
	}
	
	//operator:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator =
	(type_this const & x){
		if(std::addressof(x)==this) return *this;
		m_ring=x.m_ring;
		m_value=x.m_value;
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator =
	(type_this && x){
		m_ring=std::move(x.m_ring);
		m_value=std::move(x.m_value);
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator =
	(AuxIntType const & x){
		if(this->Is_invalid()) return *this;
		m_value= this->space_act().Is_finite() ?
			(this->space_act().modulo()==type_mod(1) ?
			type_val(0) : type_val(mod(x,this->space_act().modulo()))) :
			type_val(x);
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator +=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()) return *this=x;
		if(this->space()!=x.space()) return *this=type_this::invalid(type_val(2));
		if(this->Is_universal() && x.Is_specific()){
			this->m_ring=x.m_ring;
		}
		if(this->space_act().modulo()==type_mod(1)){
			this->m_value=type_val(0);
			return *this;
		}
		this->m_value+=x.m_value;
		if(this->space_act().Is_finite()){
			this->m_value=mod(this->m_value,this->space_act().modulo());
		}
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator -=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()) return *this=x;
		if(this->space()!=x.space()) return *this=type_this::invalid(type_val(2));
		if(this->Is_universal() && x.Is_specific()){
			this->m_ring=x.m_ring;
		}
		if(this->space_act().modulo()==type_mod(1)){
			this->m_value=type_val(0);
			return *this;
		}
		this->m_value-=x.m_value;
		if(this->space_act().Is_finite()){
			this->m_value=mod(this->m_value,this->space_act().modulo());
		}
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator *=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()) return *this=x;
		if(this->space()!=x.space()) return *this=type_this::invalid(type_val(2));
		if(this->Is_universal() && x.Is_specific()){
			this->m_ring=x.m_ring;
		}
		if(this->space_act().modulo()==type_mod(1)){
			this->m_value=type_val(0);
			return *this;
		}
		this->m_value*=x.m_value;
		if(this->space_act().Is_finite()){
			this->m_value=mod(this->m_value,this->space_act().modulo());
		}
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator *=
	(AuxIntType const & num_times){
		if(this->Is_invalid()) return *this;
		if(num_times==AuxIntType(0)) return *this=Liuze::math::zero(*this);
		this->m_value*=num_times;
		if(this->space_act().Is_finite()){
			this->m_value=Liuze::math::mod(this->m_value,this->space_act().modulo());
		}
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator /=
	(type_this const & x){
		return (*this)*=(x.inv());
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator ++
	(){
		if(this->Is_invalid()) return *this;
		++this->m_value;
		if(this->space_act().Is_finite() && this->m_value==type_val(this->space_act().modulo())){
			this->m_value=type_val(0);
		}
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this &
	Element<IntType,RingType,TraitsCheck>::operator --
	(){
		if(this->Is_invalid()) return *this;
		if(this->space_act().Is_finite() && this->m_value==type_val(0)){
			this->m_value=type_val(this->space_act().modulo());
		}
		--this->m_value;
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::operator ++
	(int){
		type_this res=*this;
		++(*this);
		return res;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::operator --
	(int){
		type_this res=*this;
		--(*this);
		return res;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::operator +
	() const {
		return *this;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::operator -
	() const {
		if(this->Is_invalid() || this->m_value==type_val(0)) return *this;
		type_this res=*this;
		res.m_value= this->space_act().Is_finite() ?
			type_val(res.space_act().modulo())-res.m_value : -res.m_value;
		return res;
	}
	
	template<typename AuxIntType,typename AuxRingType>
	inline
	Element<AuxIntType,AuxRingType>
	operator +
	(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1){
		return x0+=x1;
	}
	
	template<typename AuxIntType,typename AuxRingType>
	inline
	Element<AuxIntType,AuxRingType>
	operator -
	(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1){
		return x0-=x1;
	}
	
	template<typename AuxIntType,typename AuxRingType>
	inline
	Element<AuxIntType,AuxRingType>
	operator *
	(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1){
		return x0*=x1;
	}
	
	template<typename AuxIntType,typename AuxRingType>
	inline
	Element<AuxIntType,AuxRingType>
	operator /
	(Element<AuxIntType,AuxRingType> x0,Element<AuxIntType,AuxRingType> const & x1){
		return x0/=x1;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::operator ==
	(type_this const & x) const {
		return (this->Is_invalid() || x.Is_invalid()) ? (this->Is_invalid() && x.Is_invalid()) :
			(this->m_value==x.m_value && this->space()==x.space());
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::operator !=
	(type_this const & x) const {
		return !(*this==x);
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::operator !
	() const {
		return !Type_Bool(*this);
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Element<IntType,RingType,TraitsCheck>::operator Type_Bool
	() const {
		return this->Is_valid() && !Liuze::math::Is_approx_zero(this->m_value);
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxType,typename>
	inline
	Element<IntType,RingType,TraitsCheck>::operator AuxType
	() const {
		return this->Is_valid() ? AuxType(this->m_value) : AuxType(type_val(0));
	}
	
	//Is:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_valid
	() const {
		return this->space_act().Is_valid();
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_invalid
	() const {
		return this->space_act().Is_invalid();
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_universal
	() const {
		return this->space_act().Is_universal();
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_specific
	() const {
		return this->space_act().Is_specific();
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_zero
	() const {
		return (this->Is_valid() && this->m_value==type_val(0));
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	Type_Bool
	Element<IntType,RingType,TraitsCheck>::Is_unit
	() const {
		if(this->Is_invalid()) return false;
		if(this->Is_universal() || this->space_act().modulo()==type_mod(0)){
			return abs(this->m_value)==type_val(1);
		}
		type_val c0=0;
		if(this->space_act().modulo()==type_mod(1) || this->m_value==c0) return false;
		type_val n0=this->space_act().modulo(),n1=this->m_value;
		type_val q=n0/n1;
		type_val r=n0-q*n1;
		while(r!=c0){
			n0=n1;
			n1=r;
			q=n0/n1;
			r=n0-q*n1;
		}
		return n1==type_val(1);
	}
	
	//member value:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_val
	Element<IntType,RingType,TraitsCheck>::value
	() const {
		return this->m_value;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_space const &
	Element<IntType,RingType,TraitsCheck>::space
	() const {
		return this->m_ring;
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	inline
	typename Element<IntType,RingType,TraitsCheck>::type_space_act const &
	Element<IntType,RingType,TraitsCheck>::space_act
	() const {
		return this->m_ring.space_act();
	}
	
	//other operator:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::inv
	() const {
		if(this->Is_invalid()) return *this;
		type_val c0=0,c1=1;
		if(this->Is_universal() || this->space_act().modulo()==type_mod(0)){
			return abs(this->m_value)==c1 ? *this : type_this::invalid(type_val(1));
		}
		if(this->m_value==c0) return type_this::invalid(type_val(1));
		type_val n0=this->space_act().modulo(),n1=this->m_value,coef10=c0,coef11=c1;
		type_val q=n0/n1;
		type_val r=n0-q*n1;
		while(r!=c0){
			n0=n1;
			n1=r;
			r=coef11;
			coef11=coef10-q*coef11;
			coef10=r;
			q=n0/n1;
			r=n0-q*n1;
		}
		if(n1==c1) return type_this(coef11,*this);
		else return type_this::invalid(type_val(1));
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::pow
	(AuxIntType const & exp) const {
		if(this->Is_invalid() || (this->Is_specific() && this->space_act().modulo()==type_mod(1))){
			return *this;
		}
		if(exp<AuxIntType(0)){
			type_this res=this->inv();
			if(res.Is_invalid()) return res;
			return Liuze::math::detail::FunTemp_powint_default<type_this,AuxIntType>::powint_alt
				(res,-exp,static_cast<void*>(NULL));
		}
		return Liuze::math::detail::FunTemp_powint_default<type_this,AuxIntType>::powint_alt
			(*this,exp,static_cast<void*>(NULL));
	}
	
	//other function:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,RingType,TraitsCheck>::convert_to_value
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		typedef typename OutputIterType::value_type type_out;
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_out(*iter);
			++result;
		}
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,RingType,TraitsCheck>::convert_from_value
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_this(type_val(*iter),*this);
			++result;
		}
	}
	
	//input and output:
	
	template<typename AuxIntType,typename AuxRingType>
	inline
	std::ostream &
	operator <<
	(std::ostream & out,Element<AuxIntType,AuxRingType> const & x){
		if(x.Is_invalid()) out<<"nan";
		else out<<x.value();
		return out;
	}
	
	//static:
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::invalid
	(){
		return type_this(type_val(0),type_space::invalid());
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::invalid
	(AuxIntType const & num_value){
		return type_this(type_val(num_value),type_space::invalid());
	}
	
	template<typename IntType,typename RingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,RingType,TraitsCheck>::type_this
	Element<IntType,RingType,TraitsCheck>::universal
	(AuxIntType const & num_value){
		return type_this(type_val(num_value),type_space::universal());
	}
	
//}(class Element)
/****************************************************************************************************/
//operators of class Element{

template<typename IntType,typename RingType,typename AuxIntType,typename>
inline
Element<IntType,RingType>
operator *
(Element<IntType,RingType> x,AuxIntType const & num_times){
	return x*=num_times;
}
template<typename IntType,typename RingType,typename AuxIntType,typename>
inline
Element<IntType,RingType>
operator *
(AuxIntType const & num_times,Element<IntType,RingType> x){
	return x*=num_times;
}

//}(operators of class Element)
/****************************************************************************************************/

//(the element)
/****************************************************************************************************/

} //namespace IntegerResidueClass;
} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace Liuze{
namespace math{

//arithmetic functions for class `IntegerResidueClass::Element'{

namespace detail{
	
	template<typename IntType,typename RingType>
	class FunTemp_zero<IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef IntegerResidueClass::Element<IntType,RingType> type_val;
			static inline type_val zero(){
				return type_val(IntType(0));
			}
			static inline type_val zero(type_val const & x){
				return x.space().element_zero();
			}
	};
	
	template<typename IntType,typename RingType>
	class FunTemp_one<IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef IntegerResidueClass::Element<IntType,RingType> type_val;
			static inline type_val one(){
				return type_val(IntType(1));
			}
			static inline type_val one(type_val const & x){
				return x.space().element_one();
			}
	};
	
	template<typename IntType,typename RingType>
	class FunTemp_Is_approx_zero<IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef IntegerResidueClass::Element<IntType,RingType> type_val;
			static inline Type_Bool Is_approx_zero(type_val const & x){
				return x.Is_zero();
			}
	};
	
	template<typename IntType,typename RingType>
	class FunTemp_inverse<IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef IntegerResidueClass::Element<IntType,RingType> type_val;
			static inline type_val inverse(type_val const & x){
				return x.inv();
			}
	};
	
	template<typename IntType,typename RingType,typename ExpIntType>
	class FunTemp_powint<IntegerResidueClass::Element<IntType,RingType>,ExpIntType>{
		public:
			typedef IntegerResidueClass::Element<IntType,RingType> type_val;
			static inline type_val powint(type_val const & base,ExpIntType const & exp){
				return base.pow(exp);
			}
	};
	
} //namespace detail;

//}(arithmetic functions for class `IntegerResidueClass::Element')

} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace std{
	template<typename IntType,typename RingType>
	class numeric_limits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef Liuze::math::IntegerResidueClass::Element<IntType,RingType> type_num;
			
			static constexpr bool is_specialized=true;
			static constexpr type_num min() noexcept {return type_num::invalid();}
			static constexpr type_num max() noexcept {return type_num::invalid();}
			static constexpr type_num lowest() noexcept {return type_num::invalid();}
			static constexpr int digits=numeric_limits<typename type_num::type_val>::digits;
			static constexpr int digits10=numeric_limits<typename type_num::type_val>::digits10;
			static constexpr bool is_signed=true;
			static constexpr bool is_integer=true;
			static constexpr bool is_exact=true;
			static constexpr int radix=numeric_limits<typename type_num::type_val>::radix;
			static constexpr type_num epsilon() noexcept {return type_num::invalid();}
			static constexpr type_num round_error() noexcept {return type_num::invalid();}

			static constexpr int min_exponent=numeric_limits<typename type_num::type_val>::min_exponent;
			static constexpr int min_exponent10=numeric_limits<typename type_num::type_val>::min_exponent10;
			static constexpr int max_exponent=numeric_limits<typename type_num::type_val>::max_exponent;
			static constexpr int max_exponent10=numeric_limits<typename type_num::type_val>::max_exponent10;

			static constexpr bool has_infinity=true;
			static constexpr bool has_quiet_NaN=true;
			static constexpr bool has_signaling_NaN=false;
			static constexpr float_denorm_style has_denorm=denorm_absent;
			static constexpr bool has_denorm_loss=false;
			static constexpr type_num infinity() noexcept {return type_num::invalid(IntType(1));}
			static constexpr type_num quiet_NaN() noexcept {return type_num::invalid(IntType(1));}
			static constexpr type_num signaling_NaN() noexcept {return type_num::invalid(IntType(1));}
			static constexpr type_num denorm_min() noexcept {return type_num::invalid();}

			static constexpr bool is_iec559=false;
			static constexpr bool is_bounded=false;
			static constexpr bool is_modulo=true;

			static constexpr bool traps=true;
			static constexpr bool tinyness_before=false;
			static constexpr float_round_style round_style=round_toward_zero;
	}; //class numeric_limits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >;
} //namespace std;
/****************************************************************************************************/
/****************************************************************************************************/
#ifdef LZ_DEF_extLIB_Eigen
namespace Eigen{
	template<typename IntType,typename RingType>
	class NumTraits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >:
	public GenericNumTraits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >{
		public:
			typedef Liuze::math::IntegerResidueClass::Element<IntType,RingType> type_num;
			typedef type_num Real;
			typedef type_num NonInteger;
			typedef type_num Literal;
			typedef type_num Nested;
			
			enum{
				IsComplex=0,
				IsInteger=1,
				IsSigned=1,
				RequireInitialization=1,
				ReadCost=2*NumTraits<IntType>::ReadCost,
				AddCost=5*ReadCost+NumTraits<IntType>::AddCost+NumTraits<IntType>::MulCost,
				MulCost=5*ReadCost+2*NumTraits<IntType>::MulCost
			};
	}; //class NumTraits<Liuze::math::IntegerResidueClass::Element<IntType,RingType> >;
} //namespace Eigen;
#endif //#ifdef LZ_DEF_extLIB_Eigen
/****************************************************************************************************/

#endif //#ifndef LZ_DEF_LZ_H_math_src_intresidue_CPP
#endif //#if LZ_DEF_LZ_H_math_intresidue!=202105L