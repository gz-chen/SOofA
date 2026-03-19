#if LZ_DEF_LZ_H_math_GaloisField!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_GaloisField_CPP
#define LZ_DEF_LZ_H_math_src_GaloisField_CPP 202105L

/****************************************************************************************************/
namespace Liuze{
namespace math{
namespace GaloisField{

/****************************************************************************************************/
//the field:

/****************************************************************************************************/
//class Field_Base{
//public:
	//Constructor:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	Field_Base<IntType,CoefRingType,TraitsCheck>::Field_Base
	():
	m_ch(0),m_pow(1),m_size(0),m_irr(){
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Field_Base<IntType,CoefRingType,TraitsCheck>::Field_Base
	(AuxIntType const & num_prime,AuxIntType const & num_power,type_poly const & irr_poly):
	m_ch(0),m_pow(1),m_size(0),m_irr(){
		if(num_prime<AuxIntType(0) || num_power<AuxIntType(0) ||
		(num_power>AuxIntType(0) && num_prime==AuxIntType(0))){ //invalid;
			m_pow=2;
		} else if(num_power==AuxIntType(0) || num_prime==AuxIntType(1)){ //trivial;
			if(type_val(irr_poly.degree())==type_val(1)){ //trivial;
				m_ch=1;
				m_pow=1;
				m_size=1;
				m_irr=irr_poly;
			} else { //invalid;
				m_pow=2;
			}
		} else {
			if(AuxIntType(irr_poly.degree())==num_power){
				m_ch=type_val(num_prime);
				m_pow=type_val(num_power);
				m_size=Liuze::math::powint(m_ch,m_pow);
				m_irr=irr_poly;
			} else { //invalid;
				m_pow=2;
			}
		}
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	Field_Base<IntType,CoefRingType,TraitsCheck>::Field_Base
	(type_this const & field):
	m_ch(field.m_ch),m_pow(field.m_pow),m_size(field.m_size),m_irr(field.m_irr){
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	Field_Base<IntType,CoefRingType,TraitsCheck>::Field_Base
	(type_this && field):
	m_ch(std::move(field.m_ch)),m_pow(std::move(field.m_pow)),m_size(std::move(field.m_size)),
	m_irr(std::move(field.m_irr)){
	}
	
	//operator:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_this &
	Field_Base<IntType,CoefRingType,TraitsCheck>::operator =
	(type_this const & field){
		if(this==std::addressof(field)) return *this;
		m_ch=field.m_ch;
		m_pow=field.m_pow;
		m_size=field.m_size;
		m_irr=field.m_irr;
		return *this;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_this &
	Field_Base<IntType,CoefRingType,TraitsCheck>::operator =
	(type_this && field){
		m_ch=std::move(field.m_ch);
		m_pow=std::move(field.m_pow);
		m_size=std::move(field.m_size);
		m_irr=std::move(field.m_irr);
		return *this;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::operator ==
	(type_this const & field) const {
		return (this->Is_invalid() || field.Is_invalid()) ?
			(this->Is_invalid() && field.Is_invalid()) :
			((this->Is_universal() || field.Is_universal()) ? true :
			(this->m_ch==field.m_ch && this->m_pow==field.m_pow && this->m_size==field.m_size &&
			this->m_irr==field.m_irr));
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::operator !=
	(type_this const & field) const {
		return !(*this==field);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::operator ()
	(AuxIntType && num_value) const {
		return this->element(std::forward<AuxIntType>(num_value));
	}
	
	//Is:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::Is_valid
	() const {
		return !(this->Is_invalid());
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::Is_invalid
	() const {
		return (m_ch==type_val(0) && m_pow==type_val(2));
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::Is_universal
	() const {
		return (m_ch==type_val(0) && m_pow==type_val(1));
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	Type_Bool
	Field_Base<IntType,CoefRingType,TraitsCheck>::Is_specific
	() const {
		return (this->Is_valid() && !(this->Is_universal()));
	}
	
	//value:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_val
	Field_Base<IntType,CoefRingType,TraitsCheck>::size
	() const {
		return this->Is_specific() ? m_size : type_val(0);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_val
	Field_Base<IntType,CoefRingType,TraitsCheck>::characteristic
	() const {
		return this->Is_specific() ? m_ch : type_val(0);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_val
	Field_Base<IntType,CoefRingType,TraitsCheck>::power
	() const {
		return this->Is_specific() ? m_pow : type_val(0);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_poly
	Field_Base<IntType,CoefRingType,TraitsCheck>::irrPoly
	() const {
		return this->Is_specific() ? m_irr : type_poly();
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_poly const &
	Field_Base<IntType,CoefRingType,TraitsCheck>::irrPoly_cref
	() const {
		return m_irr;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_space_act const &
	Field_Base<IntType,CoefRingType,TraitsCheck>::space_act
	() const {
		return *this;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_space_act &
	Field_Base<IntType,CoefRingType,TraitsCheck>::space_act
	(){
		return *this;
	}
	
	//other function:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::element_zero
	() const {
		return type_element(Liuze::math::zero<type_val>(),*this);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::element_one
	() const {
		return type_element(Liuze::math::one<type_val>(),*this);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::element
	(AuxIntType && num_value) const {
		return type_element(std::forward<AuxIntType>(num_value),*this);
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_val
	Field_Base<IntType,CoefRingType,TraitsCheck>::poly_to_num
	(type_poly const & poly_value) const {
		if(!(this->Is_specific())) return type_val(0);
		type_val ord=this->power();
		typename Liuze::math::Comb_enum_pow<type_val>::type_val arr_coef_origin(ord);
		type_val i_ord,j_ord=ord;
		for(i_ord=0;i_ord<ord;++i_ord){
			--j_ord;
			arr_coef_origin[j_ord]=type_val(poly_value.coef(i_ord));
		}
		type_val res;
		Liuze::math::Comb_enum_pow<type_val>::value_to_index
			(this->characteristic(),arr_coef_origin,res);
		return res;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_poly
	Field_Base<IntType,CoefRingType,TraitsCheck>::num_to_poly
	(type_val const & num_value) const {
		if(!(this->Is_specific())) return type_poly();
		type_val ord=this->power();
		typename Liuze::math::Comb_enum_pow<type_val>::type_val arr_coef_origin;
		Liuze::math::Comb_enum_pow<type_val>::index_to_value
			(this->characteristic(),ord,num_value,arr_coef_origin);
		typename type_int_coef::type_space ring_coef
			((typename type_int_coef::type_space::type_space_act)(this->characteristic()));
		Type_ArrayTemp<type_int_coef> arr_coef(ord);
		type_val i_ord,j_ord=ord;
		for(i_ord=0;i_ord<ord;++i_ord){
			--j_ord;
			arr_coef[j_ord]=type_int_coef(arr_coef_origin[i_ord],ring_coef);
		}
		return type_poly(arr_coef);
	}
	
	//static:
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_this
	Field_Base<IntType,CoefRingType,TraitsCheck>::invalid
	(){
		type_this res;
		res.m_pow=2;
		return res;
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_this
	Field_Base<IntType,CoefRingType,TraitsCheck>::universal
	(){
		return type_this();
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::element_invalid
	(){
		return type_element::invalid();
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_element
	Field_Base<IntType,CoefRingType,TraitsCheck>::element_universal
	(AuxIntType && num_value){
		return type_element::universal(std::forward<AuxIntType>(num_value));
	}
	
	template<typename IntType,typename CoefRingType,typename TraitsCheck>
	Type_ArrayTemp<typename Field_Base<IntType,CoefRingType,TraitsCheck>::type_poly>
	Field_Base<IntType,CoefRingType,TraitsCheck>::Get_irrPoly_all
	(type_val const & num_prime,type_val const & num_degree_upper){
		typedef Type_ArrayTemp<type_poly> type_res;
		if(num_prime<type_val(2) || num_degree_upper<type_val(1)) return type_res(0);
		if(num_degree_upper==type_val(1)){
			type_res res(num_prime);
			typename type_int_coef::type_space ring_coef
				((typename type_int_coef::type_space::type_space_act)(num_prime));
			Type_ArrayTemp<type_int_coef> arr_coef(2);
			arr_coef[1]=type_int_coef(type_val(1),ring_coef);
			for(type_val i_poly=0;i_poly<num_prime;++i_poly){
				arr_coef[0]=type_int_coef(type_val(i_poly),ring_coef);
				res[i_poly]=type_poly(arr_coef);
			}
			return res;
		}
		type_val i_cand,i_coef,i_deg,i_poly;
		//the number of candidate polynomials (including (0,1)):
		type_val n_cand=Liuze::math::powint(num_prime,num_degree_upper);
		Type_ArrayTemp<type_val> flag_cand(n_cand,type_val(0)); //`0' for irreducible polynomials;
		type_res poly_cand(n_cand); //all candidate polynomials;
		Type_ArrayTemp<type_val> id_start_cand(num_degree_upper);
		Type_ArrayTemp<type_int_coef> arr_coef;
		typename type_int_coef::type_space ring_coef
			((typename type_int_coef::type_space::type_space_act)(num_prime));
		Liuze::math::Comb_enum_pow<type_val> iter_pow;
		//finding all candidate polynomials{
		i_cand=0;
		for(i_deg=1;i_deg<=num_degree_upper;++i_deg){
			//recording the starting index of `i_deg'-degree candidate polynomials:
			id_start_cand[i_deg-type_val(1)]=i_cand;
			arr_coef.resize(i_deg+type_val(1));
			arr_coef[i_deg]=type_int_coef(type_val(1),ring_coef); //monic;
			iter_pow.Set_begin(num_prime,i_deg-type_val(1));
			//0-degree coefficient is non-zero \
				if the degree of the candidate polynomial is greater than 2:
			i_poly= i_deg==type_val(1) ? type_val(0) : type_val(1);
			while(iter_pow.Is_stat_deref()){
				for(i_coef=1;i_coef<i_deg;++i_coef){
					arr_coef[i_deg-i_coef]=type_int_coef((*iter_pow)[i_coef-type_val(1)],ring_coef);
				}
				for(i_coef=i_poly;i_coef<num_prime;++i_coef){
					arr_coef[0]=type_int_coef(i_coef,ring_coef);
					poly_cand[i_cand]=type_poly(arr_coef);
					++i_cand;
				}
				++iter_pow;
			}
		}
		id_start_cand[0]=1;
		//}(finding all candidate polynomials)
		Type_ArrayTemp<type_val> arr_coef_int;
		//`fun_Get_id_cand' calculates the index of a candidate polynoimial in `poly_cand':
		std::function<type_val(type_poly const &)> fun_Get_id_cand=
			[&num_prime,&id_start_cand,&arr_coef_int]
			(type_poly const & poly)->type_val{
				type_val deg_poly_1=poly.degree()-type_val(1);
				arr_coef_int.resize(deg_poly_1);
				type_val i;
				for(i=1;i<=deg_poly_1;++i) arr_coef_int[deg_poly_1-i]=type_val(poly.coef(i));
				Liuze::math::Comb_enum_pow<type_val>::template value_to_index<type_val>
					(num_prime,arr_coef_int,i);
				return id_start_cand[deg_poly_1]
					+i*(num_prime-type_val(1))+type_val(poly.coef(0))-type_val(1);
			};
		//sieve method{
		type_val n_irr=1; //number of irreducible polynomials;
		i_poly=1;
		while(true){
			++n_irr; //`poly_cand[i_poly]' is irreducible;
			if(poly_cand[i_poly].degree()+poly_cand[i_poly].degree()>num_degree_upper){
				break;
			}
			for(i_cand=i_poly;i_cand<n_cand;++i_cand){
				if(flag_cand[i_cand]>type_val(0) && flag_cand[i_cand]<i_poly) continue;
				if(poly_cand[i_poly].degree()+poly_cand[i_cand].degree()>num_degree_upper){
					break;
				}
				flag_cand[fun_Get_id_cand(poly_cand[i_poly]*poly_cand[i_cand])]=i_poly;
			}
			//finding the next irreducible polynomial:
			do{
				++i_poly;
				if(i_poly==n_cand) goto GTS_end;
			}while(flag_cand[i_poly]>type_val(0));
		}
		for(++i_poly;i_poly<n_cand;++i_poly){
			if(flag_cand[i_poly]==type_val(0)) ++n_irr;
		}
		//}(sieve method)
		GTS_end:
		type_res res(n_irr);
		i_poly=0;
		for(i_cand=0;i_poly<n_irr;++i_cand){
			if(flag_cand[i_cand]==type_val(0)){
				res[i_poly]=poly_cand[i_cand];
				++i_poly;
			}
		}
		return res;
	}
	
//}(class Field_Base)
/****************************************************************************************************/
//class Field_tabulated{
//public:
	//Constructor:
	
	template<typename IntType,typename CoefRingType>
	Field_tabulated<IntType,CoefRingType>::Field_tabulated
	():
	type_base(),
	m_tab_add(0),m_tab_neg(0),m_tab_pow(0),m_tab_log(0){
	}
	
	template<typename IntType,typename CoefRingType>
	template<typename AuxIntType,typename>
	Field_tabulated<IntType,CoefRingType>::Field_tabulated
	(AuxIntType const & num_prime,AuxIntType const & num_power,type_poly const & irr_poly):
	type_base(num_prime,num_power,irr_poly),
	m_tab_add(0),m_tab_neg(0),m_tab_pow(0),m_tab_log(0){
		this->Cal_table();
	}
	
	template<typename IntType,typename CoefRingType>
	template<typename FieldType,typename>
	Field_tabulated<IntType,CoefRingType>::Field_tabulated
	(FieldType && field):
	m_tab_add(0),m_tab_neg(0),m_tab_pow(0),m_tab_log(0),
	type_base(field.characteristic(),field.power(),field.irrPoly_cref()){
		this->Cal_table();
	}
	
	template<typename IntType,typename CoefRingType>
	Field_tabulated<IntType,CoefRingType>::Field_tabulated
	(type_this const & field):
	type_base(field),
	m_tab_add(field.m_tab_add),m_tab_neg(field.m_tab_neg),
	m_tab_pow(field.m_tab_pow),m_tab_log(field.m_tab_log){
	}
	
	template<typename IntType,typename CoefRingType>
	Field_tabulated<IntType,CoefRingType>::Field_tabulated
	(type_this && field):
	m_tab_add(std::move(field.m_tab_add)),m_tab_neg(std::move(field.m_tab_neg)),
	m_tab_pow(std::move(field.m_tab_pow)),m_tab_log(std::move(field.m_tab_log)),
	type_base(std::move(field)){
	}
	
	//operator:
	
	template<typename IntType,typename CoefRingType>
	template<typename FieldType,typename>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_this &
	Field_tabulated<IntType,CoefRingType>::operator =
	(FieldType && field){
		static_cast<type_base *>(this)->operator =(std::forward<FieldType>(field));
		this->Cal_table();
		return *this;
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_this &
	Field_tabulated<IntType,CoefRingType>::operator =
	(type_this const & field){
		if(std::addressof(field)==this) return *this;
		static_cast<type_base *>(this)->operator =(static_cast<type_base const &>(field));
		m_tab_add=field.m_tab_add;
		m_tab_neg=field.m_tab_neg;
		m_tab_pow=field.m_tab_pow;
		m_tab_log=field.m_tab_log;
		return *this;
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_this &
	Field_tabulated<IntType,CoefRingType>::operator =
	(type_this && field){
		m_tab_add=std::move(field.m_tab_add);
		m_tab_neg=std::move(field.m_tab_neg);
		m_tab_pow=std::move(field.m_tab_pow);
		m_tab_log=std::move(field.m_tab_log);
		static_cast<type_base *>(this)->operator =(static_cast<type_base &&>(field));
		return *this;
	}
	
	template<typename IntType,typename CoefRingType>
	template<typename AuxIntType,typename>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::operator ()
	(AuxIntType && num_value) const {
		return this->element(std::forward<AuxIntType>(num_value));
	}
	
	//value:
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_space_act const &
	Field_tabulated<IntType,CoefRingType>::space_act
	() const {
		return *this;
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_space_act &
	Field_tabulated<IntType,CoefRingType>::space_act
	(){
		return *this;
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_val const &
	Field_tabulated<IntType,CoefRingType>::table_add
	(type_val const & id) const {
		return this->m_tab_add[id];
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_val const &
	Field_tabulated<IntType,CoefRingType>::table_neg
	(type_val const & id) const {
		return this->m_tab_neg[id];
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_val const &
	Field_tabulated<IntType,CoefRingType>::table_pow
	(type_val const & id) const {
		return this->m_tab_pow[id];
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_val const &
	Field_tabulated<IntType,CoefRingType>::table_log
	(type_val const & id) const {
		return this->m_tab_log[id];
	}
	
	//other function:
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::primitive
	() const {
		if(!(this->Is_specific()) || this->size()==type_val(1)){
			return this->element_invalid();
		}
		if(this->size()==type_val(2)) return this->element_one();
		return this->element(this->m_tab_pow[1]);
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::element_zero
	() const {
		return type_element(Liuze::math::zero<type_val>(),*this);
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::element_one
	() const {
		return type_element(Liuze::math::one<type_val>(),*this);
	}
	
	template<typename IntType,typename CoefRingType>
	template<typename AuxIntType,typename>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::element
	(AuxIntType && num_value) const {
		return type_element(std::forward<AuxIntType>(num_value),*this);
	}
	
	//static:
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_this
	Field_tabulated<IntType,CoefRingType>::invalid
	(){
		return type_this(type_base::invalid());
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_this
	Field_tabulated<IntType,CoefRingType>::universal
	(){
		return type_this(type_base::universal());
	}
	
	template<typename IntType,typename CoefRingType>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::element_invalid
	(){
		return type_element::invalid();
	}
	
	template<typename IntType,typename CoefRingType>
	template<typename AuxIntType,typename>
	inline
	typename Field_tabulated<IntType,CoefRingType>::type_element
	Field_tabulated<IntType,CoefRingType>::element_universal
	(AuxIntType && num_value){
		return type_element::universal(std::forward<AuxIntType>(num_value));
	}
	
//protected:
	
	//Cal:
	
	template<typename IntType,typename CoefRingType>
	void
	Field_tabulated<IntType,CoefRingType>::Cal_table
	(){
		typedef typename type_base::type_element type_base_element;
		if(!(this->Is_specific())){
			m_tab_add.resize(0);
			m_tab_neg.resize(0);
			m_tab_pow.resize(0);
			m_tab_log.resize(0);
			return;
		}
		type_val i_ele,j_ele,i_tab;
		Type_ArrayTemp<type_base_element> list_base_element(this->size());
		for(i_ele=0;i_ele<this->size();++i_ele){
			list_base_element[i_ele]=type_base_element
				(i_ele,static_cast<typename type_base_element::type_space>(*this));
		}
		m_tab_add.resize((this->size()*(this->size()+type_val(1)))/type_val(2));
		m_tab_neg.resize(this->size());
		i_tab=0;
		for(i_ele=0;i_ele<this->size();++i_ele){
			for(j_ele=0;j_ele<=i_ele;++j_ele){
				m_tab_add[i_tab]=type_val(list_base_element[i_ele]+list_base_element[j_ele]);
				++i_tab;
			}
			m_tab_neg[i_ele]=type_val(-list_base_element[i_ele]);
		}
		m_tab_pow.resize(this->size()-type_val(1));
		m_tab_log.resize(this->size()-type_val(1));
		if(this->size()==type_val(1)) return;
		m_tab_pow[0]=1;
		m_tab_log[0]=0;
		if(this->size()==type_val(2)) return;
		type_base_element ele_pow;
		for(i_ele=2;i_ele<this->size();++i_ele){
			ele_pow=list_base_element[1];
			for(j_ele=1;j_ele<m_tab_pow.size();++j_ele){
				ele_pow*=list_base_element[i_ele];
				if(type_val(ele_pow)==type_val(1)) break;
				m_tab_pow[j_ele]=type_val(ele_pow);
			}
			if(j_ele==m_tab_pow.size()) break;
		}
		for(i_ele=1;i_ele<m_tab_pow.size();++i_ele){
			m_tab_log[m_tab_pow[i_ele]-type_val(1)]=i_ele;
		}
		return;
	}
	
//}(class Field_tabulated)
/****************************************************************************************************/
//class Field_const_shared{
//public:
	//Constructor:
	
	template<typename FieldType,typename TraitsCheck>
	Field_const_shared<FieldType,TraitsCheck>::Field_const_shared
	():
	m_ptr_field(type_this::val_ptr_field_universal){
	}
	
	template<typename FieldType,typename TraitsCheck>
	template<typename AuxFieldType,typename AuxAuxFieldType>
	Field_const_shared<FieldType,TraitsCheck>::Field_const_shared
	(AuxFieldType const & field):
	m_ptr_field(field.Is_invalid() ? type_this::val_ptr_field_invalid :
	(field.Is_universal() ? type_this::val_ptr_field_universal :
	std::make_shared<type_space>(AuxAuxFieldType(field)))){
	}
	
	template<typename FieldType,typename TraitsCheck>
	Field_const_shared<FieldType,TraitsCheck>::Field_const_shared
	(type_ptr_space const & ptr_field):
	m_ptr_field((!ptr_field || ptr_field->Is_universal()) ?
	type_this::val_ptr_field_universal :
	(ptr_field->Is_invalid() ? type_this::val_ptr_field_invalid : ptr_field)){
	}
	
	template<typename FieldType,typename TraitsCheck>
	Field_const_shared<FieldType,TraitsCheck>::Field_const_shared
	(type_this const & field):
	m_ptr_field(field.m_ptr_field){
	}
	
	template<typename FieldType,typename TraitsCheck>
	Field_const_shared<FieldType,TraitsCheck>::Field_const_shared
	(type_this && field):
	m_ptr_field(std::move(field.m_ptr_field)){
	}
	
	//operator:
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this &
	Field_const_shared<FieldType,TraitsCheck>::operator =
	(type_this const & field){
		if(std::addressof(field)==this) return *this;
		m_ptr_field=field.m_ptr_field;
		return *this;
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this &
	Field_const_shared<FieldType,TraitsCheck>::operator =
	(type_this && field){
		m_ptr_field=std::move(field.m_ptr_field);
		return *this;
	}
	
	template<typename FieldType,typename TraitsCheck>
	template<typename AuxFieldType,typename AuxAuxFieldType>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this &
	Field_const_shared<FieldType,TraitsCheck>::operator =
	(AuxFieldType const & field){
		m_ptr_field= field.Is_invalid() ? type_this::val_ptr_field_invalid :
			(field.Is_universal() ? type_this::val_ptr_field_universal :
			std::make_shared<type_space>(AuxAuxFieldType(field)));
		return *this;
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this &
	Field_const_shared<FieldType,TraitsCheck>::operator =
	(type_ptr_space const & ptr_field){
		m_ptr_field= (!ptr_field || ptr_field->Is_universal()) ?
			type_this::val_ptr_field_universal :
			(ptr_field->Is_invalid() ? type_this::val_ptr_field_invalid : ptr_field);
		return *this;
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::operator ==
	(type_this const & field) const {
		return this->space()==field.space();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::operator !=
	(type_this const & field) const {
		return !(*this==field);
	}
	
	template<typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::operator ()
	(AuxIntType && num_value) const {
		return this->element(std::forward<AuxIntType>(num_value));
	}
	
	//Is:
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::Is_valid
	() const {
		return this->space().Is_valid();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::Is_invalid
	() const {
		return this->space().Is_invalid();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::Is_universal
	() const {
		return this->space().Is_universal();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Field_const_shared<FieldType,TraitsCheck>::Is_specific
	() const {
		return this->space().Is_specific();
	}
	
	//value:
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_val
	Field_const_shared<FieldType,TraitsCheck>::size
	() const {
		return this->space().size();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_val
	Field_const_shared<FieldType,TraitsCheck>::characteristic
	() const {
		return this->space().characteristic();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_val
	Field_const_shared<FieldType,TraitsCheck>::power
	() const {
		return this->space().power();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_poly
	Field_const_shared<FieldType,TraitsCheck>::irrPoly
	() const {
		return this->space().irrPoly();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_poly const &
	Field_const_shared<FieldType,TraitsCheck>::irrPoly_cref
	() const {
		return this->space().irrPoly_cref();
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_space const &
	Field_const_shared<FieldType,TraitsCheck>::space
	() const {
		return *(this->m_ptr_field);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_space_act const &
	Field_const_shared<FieldType,TraitsCheck>::space_act
	() const {
		return this->space().space_act();
	}
	
	//other function:
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::element_zero
	() const {
		return type_element(Liuze::math::zero<type_val>(),*this);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::element_one
	() const {
		return type_element(Liuze::math::one<type_val>(),*this);
	}
	
	template<typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::element
	(AuxIntType && num_value) const {
		return type_element(std::forward<AuxIntType>(num_value),*this);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_val
	Field_const_shared<FieldType,TraitsCheck>::poly_to_num
	(type_poly const & poly_value) const {
		return this->space().poly_to_num(poly_value);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_poly
	Field_const_shared<FieldType,TraitsCheck>::num_to_poly
	(type_val const & num_value) const {
		return this->space().num_to_poly(num_value);
	}
	
	//static:
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this
	Field_const_shared<FieldType,TraitsCheck>::invalid
	(){
		return type_this(type_this::val_ptr_field_invalid);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_this
	Field_const_shared<FieldType,TraitsCheck>::universal
	(){
		return type_this(type_this::val_ptr_field_universal);
	}
	
	template<typename FieldType,typename TraitsCheck>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::element_invalid
	(){
		return type_element::invalid();
	}
	
	template<typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Field_const_shared<FieldType,TraitsCheck>::type_element
	Field_const_shared<FieldType,TraitsCheck>::element_universal
	(AuxIntType && num_value){
		return type_element::universal(std::forward<AuxIntType>(num_value));
	}
	
//protected:
	//static member value:
	
	template<typename FieldType,typename TraitsCheck>
	typename Field_const_shared<FieldType,TraitsCheck>::type_ptr_space
	Field_const_shared<FieldType,TraitsCheck>::val_ptr_field_invalid
	=std::make_shared<type_space>(type_space::invalid());
	
	template<typename FieldType,typename TraitsCheck>
	typename Field_const_shared<FieldType,TraitsCheck>::type_ptr_space
	Field_const_shared<FieldType,TraitsCheck>::val_ptr_field_universal
	=std::make_shared<type_space>(type_space::universal());
	
//}(class Field_const_shared)
/****************************************************************************************************/

//(the field)
/****************************************************************************************************/
/****************************************************************************************************/
//the element:

/****************************************************************************************************/
//class Element{
//public:
	//Constructor:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	Element<IntType,FieldType,TraitsCheck>::Element
	():
	m_field(type_space::universal()),m_val_num(0),m_val_poly(){
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,TraitsCheck>::Element
	(AuxIntType const & num_value):
	m_field(type_space::universal()),m_val_num(num_value),m_val_poly(){
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,TraitsCheck>::Element
	(AuxIntType const & num_value,type_space const & field):
	m_field(field),m_val_num(num_value),m_val_poly(){
		this->Set_init_val();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,TraitsCheck>::Element
	(AuxIntType const & num_value,type_this const & element):
	m_field(element.m_field),m_val_num(num_value),m_val_poly(){
		this->Set_init_val();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	Element<IntType,FieldType,TraitsCheck>::Element
	(type_this const & element):
	m_field(element.m_field),m_val_num(element.m_val_num),m_val_poly(element.m_val_poly){
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	Element<IntType,FieldType,TraitsCheck>::Element
	(type_this && element):
	m_field(std::move(element.m_field)),
	m_val_num(std::move(element.m_val_num)),m_val_poly(std::move(element.m_val_poly)){
	}
	
	//operator:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator =
	(type_this const & element){
		if(this==std::addressof(element)) return *this;
		m_field=element.m_field;
		m_val_num=element.m_val_num;
		m_val_poly=element.m_val_poly;
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator =
	(type_this && element){
		m_field=std::move(element.m_field);
		m_val_num=std::move(element.m_val_num);
		m_val_poly=std::move(element.m_val_poly);
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator =
	(AuxIntType const & num_value){
		if(this->Is_invalid() || m_val_num==type_val(num_value)) return *this;
		m_val_num=type_val(num_value);
		this->Set_init_val();
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator +=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()){
			this->Set_invalid();
			return *this;
		}
		if(this->Is_universal() && x.Is_universal()){
			if(x.Is_zero()) return *this;
			if(this->Is_zero()) return *this=x;
			this->Set_invalid();
			return *this;
		}
		if(!(x.Is_specific())){
			type_this x_new(x.m_val_num,this->space());
			this->m_val_poly+=x_new.m_val_poly;
		} else {
			if(!(this->Is_specific())){
				this->m_field=x.m_field;
				this->Set_init_val();
			} else if(this->space()!=x.space()){
				this->Set_invalid();
				return *this;
			}
			this->m_val_poly+=x.m_val_poly;
		}
		this->m_val_num=this->space_act().poly_to_num(this->m_val_poly);
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator -=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()){
			this->Set_invalid();
			return *this;
		}
		if(this->Is_universal() && x.Is_universal()){
			if(x.Is_zero()) return *this;
			if(this->Is_zero()) return *this=-x;
			this->Set_invalid();
			return *this;
		}
		if(!(x.Is_specific())){
			type_this x_new(x.m_val_num,this->space());
			this->m_val_poly-=x_new.m_val_poly;
		} else {
			if(!(this->Is_specific())){
				this->m_field=x.m_field;
				this->Set_init_val();
			} else if(this->space()!=x.space()){
				this->Set_invalid();
				return *this;
			}
			this->m_val_poly-=x.m_val_poly;
		}
		this->m_val_num=this->space_act().poly_to_num(this->m_val_poly);
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator *=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()){
			this->Set_invalid();
			return *this;
		}
		if(this->Is_universal() && x.Is_universal()){
			if(x.Is_one() || this->Is_zero()) return *this;
			if(this->Is_one() || x.Is_zero()) return *this=x;
			if((-x).Is_one()) return *this=-(*this);
			if((-(*this)).Is_one()) return *this=-x;
			this->Set_invalid();
			return *this;
		}
		if(!(x.Is_specific())){
			if(this->Is_zero()) return *this;
			type_this x_new(x.m_val_num,this->space());
			this->m_val_poly*=x_new.m_val_poly;
		} else {
			if(!(this->Is_specific())){
				this->m_field=x.m_field;
				this->Set_init_val();
			} else if(this->space()!=x.space()){
				this->Set_invalid();
				return *this;
			}
			if(this->Is_zero()) return *this;
			if(x.Is_zero()) return *this=x;
			this->m_val_poly*=x.m_val_poly;
		}
		this->m_val_poly%=this->space_act().irrPoly_cref();
		this->m_val_num=this->space_act().poly_to_num(this->m_val_poly);
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator *=
	(AuxIntType const num_times){
		if(this->Is_invalid() || this->Is_zero()) return *this;
		if(num_times==AuxIntType(0)) return *this=Liuze::math::zero(*this);
		AuxIntType n_times=num_times;
		if(num_times<AuxIntType(0)){
			*this=-(*this);
			if(this->Is_invalid()) return *this;
			n_times=-num_times;
		}
		if(n_times==AuxIntType(1)) return *this;
		type_this res_unit=*this;
		while(n_times>AuxIntType(1)){
			*this+=res_unit;
			--n_times;
		}
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::operator /=
	(type_this const & x){
		return (*this)*=x.inv();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::operator +
	() const {
		return *this;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::operator -
	() const {
		if(this->Is_invalid()) return *this;
		type_this res;
		res.m_field=this->m_field;
		if(this->Is_specific()){
			res.m_val_poly=-(this->m_val_poly);
			res.m_val_num=res.space_act().poly_to_num(res.m_val_poly);
		} else {
			res.m_val_num=-(this->m_val_num);
		}
		return res;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::operator ==
	(type_this const & element) const {
		return (this->Is_invalid() || element.Is_invalid()) ?
			(this->Is_invalid() && element.Is_invalid()) :
			(this->m_val_num==element.m_val_num ?
			((this->Is_universal() || element.Is_universal()) ? true :
			(this->space()==element.space())) : false);
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::operator !=
	(type_this const & element) const {
		return !(*this==element);
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::operator !
	() const {
		return (this->Is_invalid() || this->m_val_num==type_val(0));
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Element<IntType,FieldType,TraitsCheck>::operator Type_Bool
	() const {
		return (this->Is_valid() && this->m_val_num!=type_val(0));
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxType,typename>
	inline
	Element<IntType,FieldType,TraitsCheck>::operator AuxType
	() const {
		return AuxType(this->value_number());
	}
	
	//Is:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_invalid
	() const {
		return this->space_act().Is_invalid();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_valid
	() const {
		return this->space_act().Is_valid();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_universal
	() const {
		return this->space_act().Is_universal();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_specific
	() const {
		return this->space_act().Is_specific();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_zero
	() const {
		return (this->Is_valid() && this->m_val_num==type_val(0));
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	Type_Bool
	Element<IntType,FieldType,TraitsCheck>::Is_one
	() const {
		return (this->Is_valid() &&
			(this->m_val_num==type_val(1) ||
			(this->Is_specific() && this->space_act().size()==type_val(1))));
	}
	
	//value:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_val
	Element<IntType,FieldType,TraitsCheck>::value_number
	() const {
		return this->Is_valid() ? this->m_val_num : type_val(0);
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_poly
	Element<IntType,FieldType,TraitsCheck>::value_polynomial
	() const {
		return this->Is_specific() ? this->m_val_poly : type_poly();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_space const &
	Element<IntType,FieldType,TraitsCheck>::space
	() const {
		return this->m_field;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_space_act const &
	Element<IntType,FieldType,TraitsCheck>::space_act
	() const {
		return this->space().space_act();
	}
	
	//other operator:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::zero
	() const {
		return type_this(Liuze::math::zero<type_val>(),this->space());
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::one
	() const {
		return type_this(Liuze::math::one<type_val>(),this->space());
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::element
	(AuxIntType const & num_value) const {
		return type_this(type_val(num_value),this->space());
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::inv
	() const {
		if(this->Is_invalid()) return *this;
		if(this->Is_universal()){
			if(this->Is_one() || (-(*this)).Is_one()) return *this;
			return type_this::invalid();
		}
		if(this->space_act().size()==type_val(1)) return *this;
		if(this->Is_zero()) return type_this::invalid();
		type_poly p0=this->space_act().irrPoly_cref();
		type_poly p1=this->m_val_poly;
		type_poly c1=this->m_val_poly.constant
			(type_int_coef(type_val(1),
			typename type_int_coef::type_space(
			typename type_int_coef::type_space::type_space_act(
			this->space_act().characteristic()))));
			//c1=Liuze::math::one(this->m_val_poly);
		type_poly coef10=Liuze::math::zero(this->m_val_poly);
		type_poly coef11=c1;
		type_poly quot,rem;
		type_poly::div(p0,p1,std::addressof(quot),std::addressof(rem));
		while(!(rem.Is_zero())){
			p0=p1;
			p1=rem;
			rem=coef11;
			coef11=coef10-quot*coef11;
			coef10=rem;
			type_poly::div(p0,p1,std::addressof(quot),std::addressof(rem));
		}
		if(p1.order()==type_val(1)){
			type_this res;
			res.m_field=this->m_field;
			res.m_val_poly=Liuze::math::inverse(p1.coef(0))*coef11;
			res.m_val_num=res.space_act().poly_to_num(res.m_val_poly);
			return res;
		}
		return type_this::invalid();
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::pow
	(AuxIntType const & num_exp) const {
		if(this->Is_invalid()) return *this;
		if(num_exp==AuxIntType(0)) return Liuze::math::one(*this);
		if(num_exp==AuxIntType(1)) return *this;
		if(this->Is_zero()){
			return num_exp<AuxIntType(0) ? type_this::invalid() : *this;
		}
		if(this->Is_one()) return *this;
		if((-(*this)).Is_one()){
			return num_exp%AuxIntType(2)==AuxIntType(0) ? -(*this) : *this;
		}
		if(this->Is_universal()) return type_this::invalid();
		//now `*this' is specific and is not 0, -1 or 1, and `num_exp' is not 0 or 1.
		type_this res=Liuze::math::one(*this);
		type_this res_b=*this;
		AuxIntType e=num_exp;
		if(e<AuxIntType(0)){
			res_b=res_b.inv();
			if(res_b.Is_invalid()) return res_b;
			e=-e;
		}
		while(e!=AuxIntType(0)){
			while(e%AuxIntType(2)==AuxIntType(0)){
				e/=AuxIntType(2);
				res_b*=res_b;
			}
			e-=AuxIntType(1);
			res*=res_b;
		}
		return res;
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::sqrt
	() const {
		if(this->Is_invalid()) return *this;
		if(this->Is_zero() || this->Is_one()) return *this;
		if(this->Is_universal()) return type_this::invalid();
		type_this res;
		type_val i_res;
		for(i_res=type_val(2);i_res<this->space_act().size();++i_res){
			res=type_this(i_res,this->space());
			if(res*res==*this) return res;
		}
		return type_this::invalid();
	}
	
	//Set:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this &
	Element<IntType,FieldType,TraitsCheck>::Set_invalid
	(){
		this->m_field=type_space::invalid();
		this->Set_init_val();
		return *this;
	}
	
	//other function:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,TraitsCheck>::convert_to_value_number
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		typedef typename OutputIterType::value_type type_out;
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_out(iter->value_number());
			++result;
		}
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,TraitsCheck>::convert_to_value_polynomial
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		typedef typename OutputIterType::value_type type_out;
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_out(iter->value_polynomial());
			++result;
		}
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,TraitsCheck>::convert_from_value_number
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_this(type_val(*iter),*this);
			++result;
		}
	}
	
	//static:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::invalid
	(){
		return type_this(type_val(0),type_space::invalid());
	}
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,FieldType,TraitsCheck>::type_this
	Element<IntType,FieldType,TraitsCheck>::universal
	(AuxIntType const & num_value){
		return type_this(type_val(num_value),type_space::universal());
	}
	
//private:
	//other function:
	
	template<typename IntType,typename FieldType,typename TraitsCheck>
	void
	Element<IntType,FieldType,TraitsCheck>::Set_init_val
	(){
		if(this->Is_invalid()){
			m_val_num=0;
			m_val_poly=type_poly();
			return;
		}
		if(this->Is_universal()) return;
		//now `*this' is specific.
		if(m_val_num>=type_val(0)){
			m_val_num=Liuze::math::mod(m_val_num,space_act().size());
			m_val_poly=space_act().num_to_poly(m_val_num);
		} else {
			m_val_num=Liuze::math::mod(-m_val_num,space_act().size());
			m_val_poly=-(space_act().num_to_poly(m_val_num));
			m_val_num=space_act().poly_to_num(m_val_poly);
		}
		return;
	}
	
//}(class Element)
/****************************************************************************************************/
//operators of class Element{

template<typename IntType,typename FieldType>
inline
Element<IntType,FieldType> operator +
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1){
	return (x0+=x1);
}

template<typename IntType,typename FieldType>
inline
Element<IntType,FieldType> operator -
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1){
	return (x0-=x1);
}

template<typename IntType,typename FieldType>
inline
Element<IntType,FieldType> operator *
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1){
	return (x0*=x1);
}

template<typename IntType,typename FieldType,typename AuxIntType,typename>
inline
Element<IntType,FieldType> operator *
(Element<IntType,FieldType> x,AuxIntType const & num_times){
	return x*=num_times;
}
template<typename IntType,typename FieldType,typename AuxIntType,typename>
inline
Element<IntType,FieldType> operator *
(AuxIntType const & num_times,Element<IntType,FieldType> x){
	return x*=num_times;
}

template<typename IntType,typename FieldType>
inline
Element<IntType,FieldType> operator /
(Element<IntType,FieldType> x0,Element<IntType,FieldType> const & x1){
	return (x0/=x1);
}

template<typename IntType,typename FieldType>
inline
std::ostream & operator <<
(std::ostream & out,Element<IntType,FieldType> const & element){
	if(element.Is_invalid()) out<<"nan";
	else out<<element.value_number();
	return out;
}

//}(operators of class Element)
/****************************************************************************************************/
//class Element<IntType,FieldType,tag_type_space_tabulated>{
//public:
	//Constructor:
	
	template<typename IntType,typename FieldType>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	():
	m_field(type_space::universal()),m_val_num(0){
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	(AuxIntType const & num_value):
	m_field(type_space::universal()),m_val_num(num_value){
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	(AuxIntType const & num_value,type_space const & field):
	m_field(field),m_val_num(num_value){
		this->Set_init_val();
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	(AuxIntType const & num_value,type_this const & element):
	m_field(element.m_field),m_val_num(num_value){
		this->Set_init_val();
	}
	
	template<typename IntType,typename FieldType>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	(type_this const & element):
	m_field(element.m_field),m_val_num(element.m_val_num){
	}
	
	template<typename IntType,typename FieldType>
	Element<IntType,FieldType,tag_type_space_tabulated>::Element
	(type_this && element):
	m_field(std::move(element.m_field)),m_val_num(std::move(element.m_val_num)){
	}
	
	//operator:
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator =
	(type_this const & element){
		if(this==std::addressof(element)) return *this;
		m_field=element.m_field;
		m_val_num=element.m_val_num;
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator =
	(type_this && element){
		m_field=std::move(element.m_field);
		m_val_num=std::move(element.m_val_num);
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator =
	(AuxIntType const & num_value){
		if(this->Is_invalid() || m_val_num==type_val(num_value)) return *this;
		m_val_num=type_val(num_value);
		this->Set_init_val();
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator +=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()){
			this->Set_invalid();
			return *this;
		}
		if(this->Is_universal() && x.Is_universal()){
			if(x.Is_zero()) return *this;
			if(this->Is_zero()) return *this=x;
			this->Set_invalid();
			return *this;
		}
		type_val val_x;
		if(!(x.Is_specific())){
			val_x=type_this(x.m_val_num,this->space()).m_val_num;
		} else {
			if(!(this->Is_specific())){
				this->m_field=x.m_field;
				this->Set_init_val();
			} else if(this->space()!=x.space()){
				this->Set_invalid();
				return *this;
			}
			val_x=x.m_val_num;
		}
		type_val val_y;
		//make `val_x>=val_y':
		if(this->m_val_num>=val_x){
			val_y=val_x;
			val_x=this->m_val_num;
		} else {
			val_y=this->m_val_num;
		}
		this->m_val_num=this->space_act().table_add(((val_x+type_val(1))*val_x)/type_val(2)+val_y);
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator -=
	(type_this const & x){
		return (*this)+=(-x);
	}
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator *=
	(type_this const & x){
		if(this->Is_invalid()) return *this;
		if(x.Is_invalid()){
			this->Set_invalid();
			return *this;
		}
		if(this->Is_universal() && x.Is_universal()){
			if(x.Is_one() || this->Is_zero()) return *this;
			if(this->Is_one() || x.Is_zero()) return *this=x;
			if((-x).Is_one()) return *this=-(*this);
			if((-(*this)).Is_one()) return *this=-x;
			this->Set_invalid();
			return *this;
		}
		type_val pow_x;
		if(!(x.Is_specific())){
			if(this->Is_zero()) return *this;
			pow_x=type_this(x.m_val_num,this->space()).m_val_num;
			if(pow_x==type_val(0)){
				this->m_val_num=pow_x;
				return *this;
			}
			pow_x=this->space_act().table_log(pow_x-type_val(1));
		} else {
			if(!(this->Is_specific())){
				this->m_field=x.m_field;
				this->Set_init_val();
			} else if(this->space()!=x.space()){
				this->Set_invalid();
				return *this;
			}
			if(this->Is_zero()) return *this;
			if(x.Is_zero()) return *this=x;
			pow_x=this->space_act().table_log(x.m_val_num-type_val(1));
		}
		type_val pow_this=this->space_act().table_log(this->m_val_num-type_val(1));
		this->m_val_num=this->space_act().table_pow
			(Liuze::math::mod(pow_this+pow_x,this->space_act().size()-type_val(1)));
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator *=
	(AuxIntType const num_times){
		if(this->Is_invalid() || this->Is_zero()) return *this;
		if(num_times==AuxIntType(0)) return *this=Liuze::math::zero(*this);
		AuxIntType n_times=num_times;
		if(num_times<AuxIntType(0)){
			*this=-(*this);
			if(this->Is_invalid()) return *this;
			n_times=-num_times;
		}
		if(n_times==AuxIntType(1)) return *this;
		type_this res_unit=*this;
		while(n_times>AuxIntType(1)){
			*this+=res_unit;
			--n_times;
		}
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::operator /=
	(type_this const & x){
		return (*this)*=x.inv();
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::operator +
	() const {
		return *this;
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::operator -
	() const {
		if(this->Is_invalid()) return *this;
		type_this res;
		res.m_field=this->m_field;
		if(this->Is_specific()){
			res.m_val_num=this->space_act().table_neg(this->m_val_num);
		} else {
			res.m_val_num=-(this->m_val_num);
		}
		return res;
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::operator ==
	(type_this const & element) const {
		return (this->Is_invalid() || element.Is_invalid()) ?
			(this->Is_invalid() && element.Is_invalid()) :
			(this->m_val_num==element.m_val_num ?
			((this->Is_universal() || element.Is_universal()) ? true :
			(this->space()==element.space())) : false);
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::operator !=
	(type_this const & element) const {
		return !(*this==element);
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::operator !
	() const {
		return (this->Is_invalid() || this->m_val_num==type_val(0));
	}
	
	template<typename IntType,typename FieldType>
	inline
	Element<IntType,FieldType,tag_type_space_tabulated>::operator Type_Bool
	() const {
		return (this->Is_valid() && this->m_val_num!=type_val(0));
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxType,typename>
	inline
	Element<IntType,FieldType,tag_type_space_tabulated>::operator AuxType
	() const {
		return AuxType(this->value_number());
	}
	
	//Is:
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_invalid
	() const {
		return space_act().Is_invalid();
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_valid
	() const {
		return space_act().Is_valid();
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_universal
	() const {
		return space_act().Is_universal();
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_specific
	() const {
		return space_act().Is_specific();
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_zero
	() const {
		return (this->Is_valid() && this->m_val_num==type_val(0));
	}
	
	template<typename IntType,typename FieldType>
	inline
	Type_Bool
	Element<IntType,FieldType,tag_type_space_tabulated>::Is_one
	() const {
		return (this->Is_valid() &&
			(this->m_val_num==type_val(1) ||
			(this->Is_specific() && this->space_act().size()==type_val(1))));
	}
	
	//value:
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_val
	Element<IntType,FieldType,tag_type_space_tabulated>::value_number
	() const {
		return this->Is_valid() ? this->m_val_num : type_val(0);
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_poly
	Element<IntType,FieldType,tag_type_space_tabulated>::value_polynomial
	() const {
		return this->Is_specific() ?
			this->space_act().num_to_poly(this->m_val_num) : type_poly();
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_space const &
	Element<IntType,FieldType,tag_type_space_tabulated>::space
	() const {
		return this->m_field;
	}
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_space_act const &
	Element<IntType,FieldType,tag_type_space_tabulated>::space_act
	() const {
		return this->space().space_act();
	}
	
	//other operator:
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::zero
	() const {
		return type_this(Liuze::math::zero<type_val>(),this->space());
	}
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::one
	() const {
		return type_this(Liuze::math::one<type_val>(),this->space());
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::element
	(AuxIntType const & num_value) const {
		return type_this(type_val(num_value),this->space());
	}
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::inv
	() const {
		if(this->Is_invalid()) return *this;
		if(this->Is_universal()){
			if(this->Is_one() || (-(*this)).Is_one()) return *this;
			return type_this::invalid();
		}
		if(this->space_act().size()==type_val(1)) return *this;
		if(this->Is_zero()) return type_this::invalid();
		type_val log_this=this->space_act().table_log(this->m_val_num-type_val(1));
		if(log_this>type_val(0)){
			log_this=this->space_act().size()-type_val(1)-log_this;
		}
		return type_this(this->space_act().table_pow(log_this),this->space());
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::pow
	(AuxIntType const & num_exp) const {
		if(this->Is_invalid()) return *this;
		if(num_exp==AuxIntType(0)) return Liuze::math::one(*this);
		if(num_exp==AuxIntType(1)) return *this;
		if(this->Is_zero()){
			return num_exp<AuxIntType(0) ? type_this::invalid() : *this;
		}
		if(this->Is_one()) return *this;
		if((-(*this)).Is_one()){
			return num_exp%AuxIntType(2)==AuxIntType(0) ? -(*this) : *this;
		}
		if(this->Is_universal()) return type_this::invalid();
		//now `*this' is specific and is not 0, -1 or 1, and `num_exp' is not 0 or 1.
		return type_this(
			this->space_act().table_pow(Liuze::math::mod(
			(this->space_act().table_log(this->m_val_num-type_val(1)))*type_val(num_exp),
			this->space_act().size()-type_val(1))),
			this->space());
	}
	
	template<typename IntType,typename FieldType>
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::sqrt
	() const {
		if(this->Is_invalid()) return *this;
		if(this->Is_zero() || this->Is_one()) return *this;
		if(this->Is_universal()) return type_this::invalid();
		type_val log_this=this->space_act().table_log(this->m_val_num-type_val(1));
		if(Liuze::math::mod(log_this,type_val(2))==type_val(0)){
			log_this/=type_val(2);
		} else {
			if(this->space_act().characteristic()==type_val(2)){
				log_this=(log_this+this->space_act().size()-type_val(1))/type_val(2);
			} else {
				return type_this::invalid();
			}
		}
		return type_this(this->space_act().table_pow(log_this),this->space());
	}
	
	//Set:
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this &
	Element<IntType,FieldType,tag_type_space_tabulated>::Set_invalid
	(){
		this->m_field=type_space::invalid();
		this->Set_init_val();
		return *this;
	}
	
	//other function:
	
	template<typename IntType,typename FieldType>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,tag_type_space_tabulated>::convert_to_value_number
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		typedef typename OutputIterType::value_type type_out;
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_out(iter->value_number());
			++result;
		}
	}
	
	template<typename IntType,typename FieldType>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,tag_type_space_tabulated>::convert_to_value_polynomial
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		typedef typename OutputIterType::value_type type_out;
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_out(iter->value_polynomial());
			++result;
		}
	}
	
	template<typename IntType,typename FieldType>
	template<typename InputIterType,typename OutputIterType,typename>
	void
	Element<IntType,FieldType,tag_type_space_tabulated>::convert_from_value_number
	(InputIterType const & first,InputIterType const & last,OutputIterType result) const {
		for(InputIterType iter=first;iter!=last;++iter){
			*result=type_this(type_val(*iter),*this);
			++result;
		}
	}
	
	//static:
	
	template<typename IntType,typename FieldType>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::invalid
	(){
		return type_this(type_val(0),type_space::invalid());
	}
	
	template<typename IntType,typename FieldType>
	template<typename AuxIntType,typename>
	inline
	typename Element<IntType,FieldType,tag_type_space_tabulated>::type_this
	Element<IntType,FieldType,tag_type_space_tabulated>::universal
	(AuxIntType const & num_value){
		return type_this(type_val(num_value),type_space::universal());
	}
	
//private:
	//other function:
	
	template<typename IntType,typename FieldType>
	inline
	void
	Element<IntType,FieldType,tag_type_space_tabulated>::Set_init_val
	(){
		if(this->Is_invalid()){
			m_val_num=0;
			return;
		}
		if(this->Is_universal()) return;
		//now `*this' is specific.
		if(m_val_num>=type_val(0)){
			m_val_num=Liuze::math::mod(m_val_num,space_act().size());
		} else {
			m_val_num=space_act().table_neg(Liuze::math::mod(-m_val_num,space_act().size()));
		}
		return;
	}
	
//}(class Element<IntType,FieldType,tag_type_space_tabulated>)
/****************************************************************************************************/

//(the element)
/****************************************************************************************************/

} //namespace GaloisField;
} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace Liuze{
namespace math{

//arithmetic functions for class `GaloisField::Element'{

namespace detail{
	
	template<typename IntType,typename FieldType>
	class FunTemp_zero<GaloisField::Element<IntType,FieldType> >{
		public:
			typedef GaloisField::Element<IntType,FieldType> type_val;
			static inline type_val zero(){
				return type_val(IntType(0));
			}
			static inline type_val zero(type_val const & x){
				return x.space().element_zero();
			}
	};
	
	template<typename IntType,typename FieldType>
	class FunTemp_one<GaloisField::Element<IntType,FieldType> >{
		public:
			typedef GaloisField::Element<IntType,FieldType> type_val;
			static inline type_val one(){
				return type_val(IntType(1));
			}
			static inline type_val one(type_val const & x){
				return x.space().element_one();
			}
	};
	
	template<typename IntType,typename FieldType>
	class FunTemp_Is_approx_zero<GaloisField::Element<IntType,FieldType> >{
		public:
			typedef GaloisField::Element<IntType,FieldType> type_val;
			static inline Type_Bool Is_approx_zero(type_val const & x){
				return x.Is_zero();
			}
	};
	
	template<typename IntType,typename FieldType>
	class FunTemp_inverse<GaloisField::Element<IntType,FieldType> >{
		public:
			typedef GaloisField::Element<IntType,FieldType> type_val;
			static inline type_val inverse(type_val const & x){
				return x.inv();
			}
	};
	
	template<typename IntType,typename FieldType,typename ExpIntType>
	class FunTemp_powint<GaloisField::Element<IntType,FieldType>,ExpIntType>{
		public:
			typedef GaloisField::Element<IntType,FieldType> type_val;
			static inline type_val powint(type_val const & base,ExpIntType const & exp){
				return base.pow(exp);
			}
	};
	
} //namespace detail;

//}(arithmetic functions for class `GaloisField::Element')

} //namespace math;
} //namespace Liuze;
/****************************************************************************************************/

/****************************************************************************************************/
namespace std{
	template<typename IntType,typename FieldType>
	class numeric_limits<Liuze::math::GaloisField::Element<IntType,FieldType> >{
		public:
			typedef Liuze::math::GaloisField::Element<IntType,FieldType> type_num;
			
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
			static constexpr type_num infinity() noexcept {return type_num::invalid();}
			static constexpr type_num quiet_NaN() noexcept {return type_num::invalid();}
			static constexpr type_num signaling_NaN() noexcept {return type_num::invalid();}
			static constexpr type_num denorm_min() noexcept {return type_num::invalid();}

			static constexpr bool is_iec559=false;
			static constexpr bool is_bounded=false;
			static constexpr bool is_modulo=true;

			static constexpr bool traps=true;
			static constexpr bool tinyness_before=false;
			static constexpr float_round_style round_style=round_toward_zero;
	}; //class numeric_limits<Liuze::math::GaloisField::Element<IntType,FieldType> >;
} //namespace std;
/****************************************************************************************************/
/****************************************************************************************************/
#ifdef LZ_DEF_extLIB_Eigen
namespace Eigen{
	template<typename IntType,typename FieldType,typename FieldTypeTag>
	class NumTraits<Liuze::math::GaloisField::Element<IntType,FieldType,FieldTypeTag> >:
	public GenericNumTraits<Liuze::math::GaloisField::Element<IntType,FieldType,FieldTypeTag> >{
		public:
			typedef Liuze::math::GaloisField::Element<IntType,FieldType,FieldTypeTag> type_num;
			typedef typename type_num::type_int_coef type_num_int_coef;
			typedef type_num Real;
			typedef type_num NonInteger;
			typedef type_num Literal;
			typedef type_num Nested;
			
			enum{
				IsComplex=0,
				IsInteger=1,
				IsSigned=1,
				RequireInitialization=1,
				ReadCost=10*NumTraits<type_num_int_coef>::ReadCost,
				AddCost=10*NumTraits<type_num_int_coef>::AddCost,
				MulCost=200*NumTraits<type_num_int_coef>::MulCost+200*AddCost
			};
	}; //class NumTraits<Liuze::math::GaloisField::Element<IntType,FieldType,FieldTypeTag> >;
	template<typename IntType,typename FieldType>
	class NumTraits<Liuze::math::GaloisField::Element<
		IntType,FieldType,Liuze::math::GaloisField::tag_type_space_tabulated> >:
	public GenericNumTraits<Liuze::math::GaloisField::Element<
	IntType,FieldType,Liuze::math::GaloisField::tag_type_space_tabulated> >{
		public:
			typedef
				Liuze::math::GaloisField::Element<
					IntType,FieldType,Liuze::math::GaloisField::tag_type_space_tabulated>
				type_num;
			typedef typename type_num::type_int_coef type_num_int_coef;
			typedef type_num Real;
			typedef type_num NonInteger;
			typedef type_num Literal;
			typedef type_num Nested;
			
			enum{
				IsComplex=0,
				IsInteger=1,
				IsSigned=1,
				RequireInitialization=1,
				ReadCost=NumTraits<type_num_int_coef>::ReadCost,
				AddCost=5*
					(NumTraits<type_num_int_coef>::AddCost+NumTraits<type_num_int_coef>::MulCost),
				MulCost=AddCost
			};
	};
	//class NumTraits<Liuze::math::GaloisField::Element<\
		IntType,FieldType,Liuze::math::GaloisField::tag_type_space_tabulated> >;
} //namespace Eigen;
#endif //#ifdef LZ_DEF_extLIB_Eigen
/****************************************************************************************************/

#endif //#ifndef LZ_DEF_LZ_H_math_src_GaloisField_CPP
#endif //#if LZ_DEF_LZ_H_math_GaloisField!=202105L