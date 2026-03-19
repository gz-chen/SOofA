#if LZ_DEF_LZ_H_stats_DOE_design_uniform!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_design_uniform_CPP
#define LZ_DEF_LZ_H_stats_src_DOE_DOE_design_uniform_CPP 202105L

namespace Liuze{
namespace stats{

/****************************************************************************************************/
//class Iterator_design_canon_sequential{
//public:
	//Constructor:
	
	template<typename RealType>
	Iterator_design_canon_sequential<RealType>::Iterator_design_canon_sequential
	():
	dim(0),stat(this->val_stat_null),point(){
	}
	
	template<typename RealType>
	Iterator_design_canon_sequential<RealType>::Iterator_design_canon_sequential
	(type_size const & dimension){
		if(dimension<=(type_size)0){
			stat=val_stat_null;
			dim=(type_size)0;
			point.resize(0,0);
		} else {
			stat=val_stat_deref;
			dim=dimension;
			point.resize(dim,1);
		}
	}
	
	template<typename RealType>
	Iterator_design_canon_sequential<RealType>::Iterator_design_canon_sequential
	(type_this const & iter):
	dim(iter.dim),point(iter.point),stat(iter.stat){
	}
	
	template<typename RealType>
	Iterator_design_canon_sequential<RealType>::Iterator_design_canon_sequential
	(type_this && iter):
	dim(iter.dim),point(iter.point),stat(iter.stat){
	}
	
	//operator:
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_sequential<RealType>::value_type &
	Iterator_design_canon_sequential<RealType>::operator *
	() const {
		return point;
	}
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_sequential<RealType>::value_type *
	Iterator_design_canon_sequential<RealType>::operator ->
	() const {
		return std::addressof(point);
	}
	
	//Is:
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_sequential<RealType>::Is_deref
	() const {
		return stat==val_stat_deref;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_sequential<RealType>::Is_end
	() const {
		return stat==val_stat_end;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_sequential<RealType>::Is_null
	() const {
		return stat==val_stat_null;
	}
	
	//Get:
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_sequential<RealType>::type_size
	Iterator_design_canon_sequential<RealType>::Get_dim
	() const {
		return dim;
	}
	
//protected:
	
	//operator:
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_sequential<RealType>::type_this &
	Iterator_design_canon_sequential<RealType>::operator =
	(type_this const & iter){
		if(&iter==this) return *this;
		dim=iter.dim;
		point=iter.point;
		stat=iter.stat;
		return *this;
	}
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_sequential<RealType>::type_this &
	Iterator_design_canon_sequential<RealType>::operator =
	(type_this && iter){
		dim=iter.dim;
		point=iter.point;
		stat=iter.stat;
		return *this;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_sequential<RealType>::operator ==
	(type_this const & iter) const {
		return stat!=iter.stat ? false :
			(stat==val_stat_null ? true :
			(dim!=iter.dim ? false : (stat==val_stat_end ? true : point==iter.point)));
	}
	
//}(class Iterator_design_canon_sequential)
/****************************************************************************************************/
//class Iterator_design_canon_RandomUniform{
//public:
	//Constructor:
	
	template<typename RealType,typename RandUIntType>
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Iterator_design_canon_RandomUniform
	():
	type_base(),engine(NULL),distr_U((type_real)0,(type_real)1){
	}
	
	template<typename RealType,typename RandUIntType>
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Iterator_design_canon_RandomUniform
	(type_rand_eng & rand_engine):
	type_base(),engine(std::addressof(rand_engine)),distr_U((type_real)0,(type_real)1){
	}
	
	template<typename RealType,typename RandUIntType>
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Iterator_design_canon_RandomUniform
	(type_size const & dimension,type_rand_eng & rand_engine):
	type_base(dimension),engine(std::addressof(rand_engine)),distr_U((type_real)0,(type_real)1){
		this->Cal_renew();
	}
	
	template<typename RealType,typename RandUIntType>
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Iterator_design_canon_RandomUniform
	(type_this const & iter):
	engine(iter.engine),distr_U((type_real)0,(type_real)1),type_base((type_base const &)iter){
	}
	
	template<typename RealType,typename RandUIntType>
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Iterator_design_canon_RandomUniform
	(type_this && iter):
	engine(iter.engine),distr_U((type_real)0,(type_real)1),type_base((type_base &&)iter){
	}
	
	//operator:
	
	template<typename RealType,typename RandUIntType>
	inline
	typename Iterator_design_canon_RandomUniform<RealType,RandUIntType>::type_this &
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::operator =
	(type_this const & iter){
		if(&iter==this) return *this;
		this->type_base::operator =((type_base const &)iter);
		engine=iter.engine;
		return *this;
	}
	
	template<typename RealType,typename RandUIntType>
	inline
	typename Iterator_design_canon_RandomUniform<RealType,RandUIntType>::type_this &
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::operator =
	(type_this && iter){
		engine=iter.engine;
		this->type_base::operator =((type_base &&)iter);
		return *this;
	}
	
	template<typename RealType,typename RandUIntType>
	inline
	Type_Bool
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::operator ==
	(type_this const & iter) const {
		return this->type_base::operator ==((type_base const &)iter) && engine==iter.engine;
	}
	
	template<typename RealType,typename RandUIntType>
	inline
	typename Iterator_design_canon_RandomUniform<RealType,RandUIntType>::type_this &
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::operator ++
	(){
		this->Cal_renew();
		return *this;
	}
	
	//begin-end:
	
	template<typename RealType,typename RandUIntType>
	inline
	void
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Set_begin
	(){
		if(!(this->Is_null())) this->stat=this->val_stat_deref;
		this->Cal_renew();
	}
	
	template<typename RealType,typename RandUIntType>
	void
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Set_begin
	(type_size const & dimension){
		if(engine==NULL || dimension<=(type_size)0){
			this->stat=this->val_stat_null;
			this->point.resize(0,0);
		}
		this->stat=this->val_stat_deref;
		if(dimension!=this->dim){
			this->dim=dimension;
			this->point.resize(this->dim,1);
		}
		this->Cal_renew();
	}
	
	//Set:
	
	template<typename RealType,typename RandUIntType>
	inline
	void
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Set_rand_engine
	(type_rand_eng & rand_engine){
		engine=std::addressof(rand_engine);
	}
	
//protected:
	
	template<typename RealType,typename RandUIntType>
	void
	Iterator_design_canon_RandomUniform<RealType,RandUIntType>::Cal_renew
	(){
		if(this->Is_deref()){
			for(type_size i_dim=0;i_dim<this->dim;++i_dim) this->point(i_dim,0)=distr_U(*engine);
		}
	}
	
//}(class Iterator_design_canon_RandomUniform)
/****************************************************************************************************/
//class Iterator_design_canon_GoodPoint{
//public:
	//Constructor:
	
	template<typename RealType>
	Iterator_design_canon_GoodPoint<RealType>::Iterator_design_canon_GoodPoint
	():
	type_base(),gp(){
	}
	
	template<typename RealType>
	Iterator_design_canon_GoodPoint<RealType>::Iterator_design_canon_GoodPoint
	(type_matrix const & good_point):
	type_base(),gp(){
		if(good_point.rows()>(type_size)0 && good_point.cols()>(type_size)0){
			this->dim=good_point.rows();
			this->stat=type_this::val_stat_deref;
			this->gp=good_point.col(0);
			this->point=this->gp;
		}
	}
	
	template<typename RealType>
	Iterator_design_canon_GoodPoint<RealType>::Iterator_design_canon_GoodPoint
	(type_this const & iter):
	gp(iter.gp),type_base((type_base const &)iter){
	}
	
	template<typename RealType>
	Iterator_design_canon_GoodPoint<RealType>::Iterator_design_canon_GoodPoint
	(type_this && iter):
	gp(iter.gp),type_base((type_base &&)iter){
	}
	
	//operator:
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_GoodPoint<RealType>::type_this &
	Iterator_design_canon_GoodPoint<RealType>::operator =
	(type_this const & iter){
		if(&iter==this) return *this;
		this->type_base::operator =((type_base const &)iter);
		gp=iter.gp;
		return *this;
	}
	
	template<typename RealType>
	inline
	typename Iterator_design_canon_GoodPoint<RealType>::type_this &
	Iterator_design_canon_GoodPoint<RealType>::operator =
	(type_this && iter){
		gp=iter.gp;
		this->type_base::operator =((type_base &&)iter);
		return *this;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_GoodPoint<RealType>::operator ==
	(type_this const & iter) const {
		return this->type_base::operator ==((type_base const &)iter) ?
			(this->Is_deref() ? gp==iter.gp : true) :
			false;
	}
	
	template<typename RealType>
	typename Iterator_design_canon_GoodPoint<RealType>::type_this &
	Iterator_design_canon_GoodPoint<RealType>::operator ++
	(){
		if(!(this->Is_deref())) return *this;
		for(type_size i_dim=0;i_dim<this->dim;++i_dim){
			this->point(i_dim,0)+=gp(i_dim,0);
			if(this->point(i_dim,0)>=(type_real)1) this->point(i_dim,0)-=(type_real)1;
		}
		return *this;
	}
	
	//begin-end:
	
	template<typename RealType>
	inline
	void
	Iterator_design_canon_GoodPoint<RealType>::Set_begin
	(){
		if(this->Is_null()) return;
		this->stat=type_this::val_stat_deref;
		this->point=gp;
	}
	
	template<typename RealType>
	void
	Iterator_design_canon_GoodPoint<RealType>::Set_begin
	(type_matrix const & good_point){
		if(good_point.rows()<=(type_size)0 || good_point.cols()<=(type_size)0){
			this->dim=(type_size)0;
			this->stat=type_this::val_stat_null;
			this->point.resize(0,0);
			this->gp.resize(0,0);
			return;
		}
		this->dim=good_point.rows();
		this->stat=type_this::val_stat_deref;
		this->gp=good_point.col(0);
		this->point=this->gp;
	}
	
//public static functions:
	
	//Get: good point:
	
	template<typename RealType>
	typename Iterator_design_canon_GoodPoint<RealType>::type_matrix
	Iterator_design_canon_GoodPoint<RealType>::Get_GoodPoint_cyclotomic
	(type_size const & dimension,type_size prime_number){
		if(dimension<=(type_size)0) return type_matrix(0,0);
		type_matrix res(dimension,1);
		type_size prime=(type_size)2*dimension+(type_size)3;
		prime= prime_number>=prime ? prime_number : prime;
		type_real angle0=(type_real)LZ_DEF_const_math_2pi/(type_real)prime;
		for(type_size i_dim=0;i_dim<dimension;++i_dim){
			res(i_dim,0)=(type_real)2*cos(angle0*(type_real)(i_dim+1));
			res(i_dim,0)-=type_real(floor(res(i_dim,0)));
		}
		return res;
	}
	
	template<typename RealType>
	typename Iterator_design_canon_GoodPoint<RealType>::type_matrix
	Iterator_design_canon_GoodPoint<RealType>::Get_GoodPoint_power
	(type_size const & dimension,type_size prime_number,type_stat modify){
		if(dimension<=(type_size)0) return type_matrix(0,0);
		type_matrix res(dimension,1);
		type_size prime= prime_number>=(type_size)2 ? prime_number : (type_size)2;
		type_size i_dim;
		res(0,0)=type_real(pow((type_real)prime,(type_real)1/(type_real(dimension)+type_real(1))));
		if(modify==(type_stat)1){
			type_real tr=res(0,0);
			for(type_size i_dim_0=2;i_dim_0>=1;--i_dim_0){
				for(i_dim=i_dim_0;i_dim<dimension;i_dim+=2) res(i_dim,0)=(tr*=res(0,0));
			}
		} else {
			for(i_dim=1;i_dim<dimension;++i_dim) res(i_dim,0)=res(i_dim-1,0)*res(0,0);
		}
		for(i_dim=0;i_dim<dimension;++i_dim) res(i_dim,0)-=type_real(floor(res(i_dim,0)));
		return res;
	}
	
//}(class Iterator_design_canon_GoodPoint)
/****************************************************************************************************/
//class Iterator_design_canon_WrapAroundShifter{
//public:
	//Constructor:
	
	template<typename RealType>
	Iterator_design_canon_WrapAroundShifter<RealType>::Iterator_design_canon_WrapAroundShifter
	():
	type_base(),list_size(0),list_i_size(0),list_iter_begin(0),list_iter(0),list_point_buff(0){
	}
	
	template<typename RealType>
	template<typename InIterType>
	Iterator_design_canon_WrapAroundShifter<RealType>::Iterator_design_canon_WrapAroundShifter
	(Type_ArrayTemp<type_size> const & list_design_size,
	Type_ArrayTemp<InIterType> const & list_design_begin):
	type_base(){
		this->Set_begin<InIterType>(list_design_size,list_design_begin);
	}
	
	template<typename RealType>
	Iterator_design_canon_WrapAroundShifter<RealType>::Iterator_design_canon_WrapAroundShifter
	(type_this const & iter):
	type_base(static_cast<type_base const &>(iter)),
	list_size(iter.list_size),list_i_size(iter.list_i_size),
	list_iter_begin(iter.list_iter_begin),list_iter(iter.list_iter),
	list_point_buff(iter.list_point_buff){
	}
	
	template<typename RealType>
	Iterator_design_canon_WrapAroundShifter<RealType>::Iterator_design_canon_WrapAroundShifter
	(type_this && iter):
	type_base(static_cast<type_base &&>(iter)),
	list_size(std::forward<Type_ArrayTemp<type_size> >(iter.list_size)),
	list_i_size(std::forward<Type_ArrayTemp<type_size> >(iter.list_i_size)),
	list_iter_begin(std::forward<Type_ArrayTemp<type_iter_store> >(iter.list_iter_begin)),
	list_iter(std::forward<Type_ArrayTemp<type_iter_store> >(iter.list_iter)),
	list_point_buff(std::forward<Type_ArrayTemp<type_matrix> >(iter.list_point_buff)){
	}
	
	//operator:
	
	template<typename RealType>
	typename Iterator_design_canon_WrapAroundShifter<RealType>::type_this &
	Iterator_design_canon_WrapAroundShifter<RealType>::operator =
	(type_this const & iter){
		if(&iter==this) return *this;
		this->type_base::operator =((type_base const &)iter);
		list_size=iter.list_size;
		list_i_size=iter.list_i_size;
		list_iter_begin=iter.list_iter_begin;
		list_iter=iter.list_iter;
		list_point_buff=iter.list_point_buff;
		return *this;
	}
	
	template<typename RealType>
	typename Iterator_design_canon_WrapAroundShifter<RealType>::type_this &
	Iterator_design_canon_WrapAroundShifter<RealType>::operator =
	(type_this && iter){
		list_size=std::forward<Type_ArrayTemp<type_size> >(iter.list_size);
		list_i_size=std::forward<Type_ArrayTemp<type_size> >(iter.list_i_size);
		list_iter_begin=std::forward<Type_ArrayTemp<type_iter_store> >(iter.list_iter_begin);
		list_iter=std::forward<Type_ArrayTemp<type_iter_store> >(iter.list_iter);
		list_point_buff=std::forward<Type_ArrayTemp<type_matrix> >(iter.list_point_buff);
		this->type_base::operator =((type_base &&)iter);
		return *this;
	}
	
	template<typename RealType>
	inline
	Type_Bool
	Iterator_design_canon_WrapAroundShifter<RealType>::operator ==
	(type_this const & iter) const {
		return this->type_base::operator ==((type_base const &)iter) ?
			(this->Is_deref() ?
			(list_size==iter.list_size && list_i_size=iter.list_i_size &&
			list_iter_begin=iter.list_iter_begin && list_iter=iter.list_iter) :
			true) : false;
	}
	
	template<typename RealType>
	typename Iterator_design_canon_WrapAroundShifter<RealType>::type_this &
	Iterator_design_canon_WrapAroundShifter<RealType>::operator ++
	(){
		if(!(this->Is_deref())) return *this;
		type_size i=list_size.size()-(type_size)1,i_dim;
		while(list_size[i]!=type_this::val_size_infinity && list_i_size[i]+(type_size)1==list_size[i]){
			if(i==(type_size)0){
				this->stat=type_this::val_stat_end;
				return *this;
			}
			--i;
		}
		++(list_iter[i]);
		if(list_size[i]!=type_this::val_size_infinity) ++(list_i_size[i]);
		this->point= i==(type_size)0 ? type_matrix::Zero(this->dim,1) : list_point_buff[i-(type_size)1];
		while(true){
			this->point+=(*(list_iter[i]));
			list_point_buff[i]=this->point;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				if(list_point_buff[i](i_dim,0)>=(type_real)1){
					list_point_buff[i](i_dim,0)-=(type_real)1;
					this->point(i_dim,0)=list_point_buff[i](i_dim,0);
				}
			}
			++i;
			if(i==list_size.size()) break;
			list_i_size[i]=(type_size)0;
			list_iter[i]=list_iter_begin[i];
		}
		return *this;
	}
	
	//begin-end:
	
	template<typename RealType>
	void
	Iterator_design_canon_WrapAroundShifter<RealType>::Set_begin
	(){
		if(this->Is_null()) return;
		this->stat=type_this::val_stat_deref;
		list_i_size.assign(list_size.size(),(type_size)0);
		list_iter=list_iter_begin;
		this->point=type_matrix::Zero(this->dim,1);
		type_size i,i_dim;
		for(i=0;i<list_iter.size();++i){
			this->point+=(*(list_iter[i]));
			this->list_point_buff[i]=this->point;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				if(this->list_point_buff[i](i_dim,0)>=(type_real)1){
					this->list_point_buff[i](i_dim,0)-=(type_real)1;
					this->point(i_dim,0)=this->list_point_buff[i](i_dim,0);
				}
			}
		}
	}
	
	template<typename RealType>
	template<typename InIterType>
	void
	Iterator_design_canon_WrapAroundShifter<RealType>::Set_begin
	(Type_ArrayTemp<type_size> const & list_design_size,
	Type_ArrayTemp<InIterType> const & list_design_begin){
		LZ_DEF_func_check_traits((std::is_convertible<InIterType,type_iter_store>::value));
		type_size n_iter=std::min(list_design_size.size(),list_design_begin.size()),i,i_dim;
		if(n_iter>(type_size)0){
			i_dim=list_design_begin[0]->rows();
			if(i_dim<=(type_size)0) n_iter=(type_size)0;
			else for(i=1;i<n_iter;++i){
				if(list_design_begin[i]->rows()!=i_dim){
					n_iter=(type_size)0;
					break;
				}
			}
		}
		if(n_iter>(type_size)0){
			this->dim=i_dim;
			this->stat=type_this::val_stat_deref;
			this->point=type_matrix::Zero(this->dim,1);
			this->list_size.resize(n_iter);
			this->list_i_size.assign(n_iter,(type_size)0);
			this->list_iter_begin.resize(n_iter);
			this->list_iter.resize(n_iter);
			this->list_point_buff.resize(n_iter);
			for(i=0;i<n_iter;++i){
				list_size[i]=list_design_size[i];
				list_iter[i]=list_iter_begin[i]=list_design_begin[i];
				this->point+=(*(list_design_begin[i]));
				this->list_point_buff[i]=this->point;
				for(i_dim=0;i_dim<this->dim;++i_dim){
					if(this->list_point_buff[i](i_dim,0)>=(type_real)1){
						this->list_point_buff[i](i_dim,0)-=(type_real)1;
						this->point(i_dim,0)=this->list_point_buff[i](i_dim,0);
					}
				}
			}
		} else {
			this->dim=(type_size)0;
			this->stat=type_this::val_stat_null;
			this->point.resize(0,0);
			this->list_size.resize(0);
			this->list_i_size.resize(0);
			this->list_iter_begin.resize(0);
			this->list_iter.resize(0);
			this->list_point_buff.resize(0);
		}
	}
	
//}(class Iterator_design_canon_WrapAroundShifter)
/****************************************************************************************************/
/****************************************************************************************************/
//class Get_design_GoodLatticePoint{
//public:
	//Constructor:
	
	template<typename RealType,typename IntType>
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_GoodLatticePoint
	():
	size(0),dim_max(0),stat(type_this::val_stat_null),unit(0),crd_can(0){
	}
	
	template<typename RealType,typename IntType>
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_GoodLatticePoint
	(type_int const & size_design):
	size(size_design>=(type_int)2 ? size_design : (type_int)0){
		if(size==(type_int)0){
			dim_max=(type_int)0;
			stat=type_this::val_stat_null;
		} else {
			stat=type_this::val_stat_invalid;
			//the following calculates dim_max and unit.
			dim_max=(type_int)1;
			Type_ArrayTemp<type_int> unit_init(size-1);
			unit_init[0]=(type_int)1;
			type_int i;
			for(i=(type_int)2;i<size;++i){
				if(math::Is_coprime(size,i)) unit_init[dim_max++]=i;
			}
			unit.assign(unit_init.begin(),unit_init.begin()+dim_max);
			//the following calculates crd_can.
			crd_can.resize(size);
			type_real tr0=size,tr1=(type_real)1/((type_real)((type_int)2*size));
			for(i=(type_int)0;i<size;++i) crd_can[i]=(type_real)i/tr0+tr1;
		}
	}
	
	//Is:
	
	template<typename RealType,typename IntType>
	inline
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_null
	() const {
		return this->stat==type_this::val_stat_null;
	}
	
	template<typename RealType,typename IntType>
	inline
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_invalid
	() const {
		return this->stat==type_this::val_stat_invalid;
	}
	
	template<typename RealType,typename IntType>
	inline
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_valid
	() const {
		return this->stat==type_this::val_stat_valid;
	}
	
	//Get: value:
	
	template<typename RealType,typename IntType>
	inline
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_int
	Get_design_GoodLatticePoint<RealType,IntType>::Get_dim
	() const {
		return this->dim;
	}
	
	template<typename RealType,typename IntType>
	inline
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_int
	Get_design_GoodLatticePoint<RealType,IntType>::Get_size
	() const {
		return this->size;
	}
	
	template<typename RealType,typename IntType>
	inline
	Type_ArrayTemp<typename Get_design_GoodLatticePoint<RealType,IntType>::type_int>
	Get_design_GoodLatticePoint<RealType,IntType>::Get_unit
	() const {
		return this->unit;
	}
	
	template<typename RealType,typename IntType>
	inline
	Type_ArrayTemp<typename Get_design_GoodLatticePoint<RealType,IntType>::type_real>
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_unary
	() const {
		return this->crd_can;
	}
	
	template<typename RealType,typename IntType>
	inline
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_pt_fac
	Get_design_GoodLatticePoint<RealType,IntType>::Get_vec_gen
	() const {
		return this->vec_gen;
	}
	
	template<typename RealType,typename IntType>
	inline
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_real
	Get_design_GoodLatticePoint<RealType,IntType>::Get_val_crit
	() const {
		return this->val_crit;
	}
	
	//Set:
	
	template<typename RealType,typename IntType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_stat
	Get_design_GoodLatticePoint<RealType,IntType>::Set_size
	(type_int const & size_design){
		if(size_design==size || (size_design<(type_int)2 && size==(type_int)0)) return stat;
		if(size_design<(type_int)2){
			size=(type_int)0;
			dim_max=(type_int)0;
			stat=type_this::val_stat_null;
		} else {
			size=size_design;
			stat=type_this::val_stat_invalid;
			//the following calculates dim_max and unit.
			dim_max=(type_int)1;
			Type_ArrayTemp<type_int> unit_init(size-1);
			unit_init[0]=(type_int)1;
			type_int i;
			for(i=(type_int)2;i<size;++i){
				if(math::Is_coprime(size,i)) unit_init[dim_max++]=i;
			}
			unit.assign(unit_init.begin(),unit_init.begin()+dim_max);
			//the following calculates crd_can.
			crd_can.resize(size);
			type_real tr0=size,tr1=(type_real)1/((type_real)((type_int)2*size));
			for(i=(type_int)0;i<size;++i) crd_can[i]=(type_real)i/tr0+tr1;
		}
		return stat;
	}
	
	//Cal: vec_gen and val_crit:
	
	template<typename RealType,typename IntType>
	template<typename CritType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_stat
	Get_design_GoodLatticePoint<RealType,IntType>::Cal_best
	(type_int const & dim_design,CritType & crit_canon){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<CritType>::type(
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				Iterator_constval<type_real>,
				Iterator_constval<type_real>
			)>::type,
			type_real>::value));
		if(size==(type_int)0) return stat;
		if(dim_design<=(type_int)0 || dim_design>dim_max){
			stat=type_this::val_stat_invalid;
			return stat;
		}
		//the following begins to search for the best design.
		dim=dim_design;
		vec_gen.resize(dim,1);
		Type_Bool tag_first_des=true;
		type_real v_crit; //inner value of criterion;
		Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
		math::Comb_enum_comb<type_int> iter_comb(dim_max,dim);
		Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
		type_int i_dim,i_size;
		for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
		while(iter_comb.Is_stat_deref() && (*iter_comb)[0]==(type_int)0){
			for(i_size=0;i_size<size;++i_size){
				for(i_dim=0;i_dim<dim;++i_dim){
					des[i_size](i_dim,0)=crd_can[(unit[(*iter_comb)[i_dim]]*i_size)%size];
				}
			}
			v_crit=crit_canon(des.begin(),des.end(),iter_wt,iter_wt.end());
			if(v_crit<val_crit || tag_first_des){
				tag_first_des=false;
				val_crit=v_crit;
				for(i_dim=0;i_dim<dim;++i_dim) vec_gen(i_dim,0)=unit[(*iter_comb)[i_dim]];
			}
			++iter_comb;
		}
		stat=type_this::val_stat_valid;
		return stat;
	}
	
	template<typename RealType,typename IntType>
	template<typename CritType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_stat
	Get_design_GoodLatticePoint<RealType,IntType>::Cal_best_power
	(type_int const & dim_design,CritType & crit_canon){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<CritType>::type(
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				Iterator_constval<type_real>,
				Iterator_constval<type_real>
			)>::type,
			type_real>::value));
		if(size==(type_int)0) return stat;
		if(dim_design<=(type_int)0 || dim_design>dim_max){
			stat=type_this::val_stat_invalid;
			return stat;
		}
		//the following begins to search for the best design.
		dim=dim_design;
		type_int i_dim,i_size,i_unit;
		Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
		for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
		Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
		Type_Bool tag_first_des=true;
		type_real v_crit; //inner value of criterion;
		type_pt_fac v_gen(dim,1); //inner vec_gen;
		Type_ArrayTemp<type_int> val_v_gen(dim); //the sorted entries of v_gen;
		val_v_gen[0]=v_gen(0,0)=(type_int)1;
		for(i_unit=(dim_max>(type_int)1 ? (type_int)1 : (type_int)0);i_unit<dim_max;++i_unit){
			//first, check if unit[i_unit] can generate a feasible generating vector.
			for(i_dim=1;i_dim<dim;++i_dim){
				val_v_gen[i_dim]=v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
			}
			std::sort(val_v_gen.begin(),val_v_gen.end());
			for(i_dim=1;i_dim<dim;++i_dim){
				if(val_v_gen[i_dim]==val_v_gen[i_dim-1]) i_dim=dim;
			}
			if(i_dim>dim) continue;
			//second, check the criterion value of this generating vector.
			for(i_size=0;i_size<size;++i_size){
				for(i_dim=0;i_dim<dim;++i_dim){
					des[i_size](i_dim,0)=crd_can[(v_gen(i_dim,0)*i_size)%size];
				}
			}
			v_crit=crit_canon(des.begin(),des.end(),iter_wt,iter_wt.end());
			if(v_crit<val_crit || tag_first_des){
				tag_first_des=false;
				val_crit=v_crit;
				vec_gen=v_gen;
			}
		}
		stat= tag_first_des ? type_this::val_stat_invalid : type_this::val_stat_valid;
		return stat;
	}
	
	template<typename RealType,typename IntType>
	template<typename CritType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_stat
	Get_design_GoodLatticePoint<RealType,IntType>::Cal_best_power_nocheck
	(type_int const & dim_design,CritType & crit_canon){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<CritType>::type(
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				Iterator_constval<type_real>,
				Iterator_constval<type_real>
			)>::type,
			type_real>::value));
		if(size==(type_int)0) return stat;
		if(dim_design<=(type_int)0){
			stat=type_this::val_stat_invalid;
			return stat;
		}
		//the following begins to search for the best design.
		dim=dim_design;
		type_int i_dim,i_size,i_unit;
		Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
		for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
		Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
		Type_Bool tag_first_des=true;
		type_real v_crit; //inner value of criterion;
		type_pt_fac v_gen(dim,1); //inner vec_gen;
		v_gen(0,0)=(type_int)1;
		for(i_unit=(dim_max>(type_int)1 ? (type_int)1 : (type_int)0);i_unit<dim_max;++i_unit){
			for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
			for(i_size=0;i_size<size;++i_size){
				for(i_dim=0;i_dim<dim;++i_dim){
					des[i_size](i_dim,0)=crd_can[(v_gen(i_dim,0)*i_size)%size];
				}
			}
			v_crit=crit_canon(des.begin(),des.end(),iter_wt,iter_wt.end());
			if(v_crit<val_crit || tag_first_des){
				tag_first_des=false;
				val_crit=v_crit;
				vec_gen=v_gen;
			}
		}
		stat=type_this::val_stat_valid;
		return stat;
	}
	
	template<typename RealType,typename IntType>
	template<typename CritType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_stat
	Get_design_GoodLatticePoint<RealType,IntType>::Cal_best_power_LeaveOneOut_nocheck
	(type_int const & dim_design,CritType & crit_canon){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<CritType>::type(
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				typename Type_ArrayTemp<type_pt_can>::const_iterator,
				Iterator_constval<type_real>,
				Iterator_constval<type_real>
			)>::type,
			type_real>::value));
		if(size==(type_int)0) return stat;
		if(dim_design<=(type_int)0){
			stat=type_this::val_stat_invalid;
			return stat;
		}
		//the following begins to search for the best design.
		dim=dim_design;
		type_int i_dim,i_size,i_unit;
		type_real v_crit; //inner value of criterion;
		type_pt_fac v_gen(dim,1); //inner vec_gen;
		v_gen(0,0)=(type_int)1;
		type_int size_loo=size-(type_int)1; //size of the leave-one-out ("loo") GLP;
		Type_ArrayTemp<type_real> crd_can_loo(size_loo);
		Type_ArrayTemp<type_pt_can> des(size_loo); //inner canonical design;
		v_crit=(type_real)size_loo;
		for(i_size=0;i_size<size_loo;++i_size){
			des[i_size].resize(dim,1);
			crd_can_loo[i_size]=((type_real)i_size+(type_real)0.5)/v_crit;
		}
		Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size_loo,size_loo);
		Type_Bool tag_first_des=true;
		for(i_unit=(dim_max>(type_int)1 ? (type_int)1 : (type_int)0);i_unit<dim_max;++i_unit){
			for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
			for(i_size=0;i_size<size_loo;++i_size){
				for(i_dim=0;i_dim<dim;++i_dim){
					des[i_size](i_dim,0)=
						crd_can_loo[(v_gen(i_dim,0)*(i_size+(type_int)1))%size-(type_int)1];
				}
			}
			v_crit=crit_canon(des.begin(),des.end(),iter_wt,iter_wt.end());
			if(v_crit<val_crit || tag_first_des){
				tag_first_des=false;
				val_crit=v_crit;
				vec_gen=v_gen;
			}
		}
		stat=type_this::val_stat_valid;
		return stat;
	}
	
	//Get: design:
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	inline
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_factor
	(OutIterType result_begin) const {
		LZ_DEF_func_check_traits((std::is_same<type_pt_fac,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		if(!(this->Is_valid())) return result_begin;
		return type_this::Get_design_factor(size,vec_gen,result_begin);
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	inline
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_canon
	(OutIterType result_begin) const {
		LZ_DEF_func_check_traits((std::is_same<type_pt_can,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		if(!(this->Is_valid())) return result_begin;
		return type_this::Get_design_canon(size,vec_gen,result_begin);
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	inline
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_factor_LeaveOneOut
	(OutIterType result_begin) const {
		LZ_DEF_func_check_traits((std::is_same<type_pt_fac,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		if(!(this->Is_valid())) return result_begin;
		return type_this::Get_design_factor_LeaveOneOut(size,vec_gen,result_begin);
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	inline
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_canon_LeaveOneOut
	(OutIterType result_begin) const {
		LZ_DEF_func_check_traits((std::is_same<type_pt_can,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		if(!(this->Is_valid())) return result_begin;
		return type_this::Get_design_canon_LeaveOneOut(size,vec_gen,result_begin);
	}
	
	//Convert:
	
	template<typename RealType,typename IntType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_pt_can
	Get_design_GoodLatticePoint<RealType,IntType>::Convert_to_canon
	(type_pt_fac const & point) const {
		if(size==(type_int)0) return type_pt_can::Zero(point.rows(),point.cols());
		type_pt_can res(point.rows(),point.cols());
		Type_Size ir,ic;
		for(ic=(Type_Size)0;ic<res.cols();++ic){
			for(ir=(Type_Size)0;ir<res.rows();++ir){
				res(ir,ic)=crd_can[Liuze::math::mod(point(ir,ic),size)];
			}
		}
		return res;
	}
	
	template<typename RealType,typename IntType>
	typename Get_design_GoodLatticePoint<RealType,IntType>::type_pt_fac
	Get_design_GoodLatticePoint<RealType,IntType>::Convert_to_factor
	(type_pt_can const & point) const {
		if(size==(type_int)0) return type_pt_fac::Zero(point.rows(),point.cols());
		type_pt_fac res(point.rows(),point.cols());
		Type_Size ir,ic;
		for(ic=(Type_Size)0;ic<res.cols();++ic){
			for(ir=(Type_Size)0;ir<res.rows();++ir){
				res(ir,ic)=type_int(floor((point(ir,ic)-type_real(floor(point(ir,ic))))*type_real(size)));
			}
		}
		return res;
	}
	
//public static function:
	
	//Is:
	
	template<typename RealType,typename IntType>
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_null
	(type_stat const & state){
		return state==type_this::val_stat_null;
	}
	
	template<typename RealType,typename IntType>
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_invalid
	(type_stat const & state){
		return state==type_this::val_stat_invalid;
	}
	
	template<typename RealType,typename IntType>
	Type_Bool
	Get_design_GoodLatticePoint<RealType,IntType>::Is_valid
	(type_stat const & state){
		return state==type_this::val_stat_valid;
	}
	
	//Get: design:
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_factor
	(type_int const & size_design,type_pt_fac const & vector_gen,OutIterType result_begin){
		LZ_DEF_func_check_traits((std::is_same<type_pt_fac,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		type_int const dim=vector_gen.rows();
		if(size_design<=(type_int)0 || dim<=(type_int)0 || vector_gen.cols()<=(type_int)0){
			return result_begin;
		}
		//checking ends.
		type_pt_fac point=type_pt_fac::Zero(dim,1);
		type_int i_dim,i_level=0;
		while(true){
			*result_begin=point;
			++result_begin;
			++i_level;
			if(i_level==size_design) return result_begin;
			for(i_dim=0;i_dim<dim;++i_dim){
				point(i_dim,0)=(point(i_dim,0)+vector_gen(i_dim,0))%size_design;
			}
		}
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_canon
	(type_int const & size_design,type_pt_fac const & vector_gen,OutIterType result_begin){
		LZ_DEF_func_check_traits((std::is_same<type_pt_can,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		type_int const dim=vector_gen.rows();
		if(size_design<=(type_int)0 || dim<=(type_int)0 || vector_gen.cols()<=(type_int)0){
			return result_begin;
		}
		//checking ends.
		type_pt_can point(dim,1);
		Type_ArrayTemp<type_real> crd(size_design);
		Type_ArrayTemp<type_int> point_id(dim,(type_int)0); //the point in the factor design;
		type_int i_dim,i_level;
		//calculating the 1-dim coordinates:
		for(i_level=0;i_level<size_design;++i_level){
			crd[i_level]=(type_real(i_level)+type_real(0.5))/type_real(size_design);
		}
		//transformation from factor design to canonical design:
		auto fun_Get_point=[&i_dim,&dim,&point,&crd,&point_id]()->type_pt_can&{
				for(i_dim=(type_int)0;i_dim<dim;++i_dim){
					point(i_dim,0)=crd[point_id[i_dim]];
				}
				return point;
			};
		i_level=(type_int)0;
		while(true){
			*result_begin=fun_Get_point();
			++result_begin;
			++i_level;
			if(i_level==size_design) return result_begin;
			for(i_dim=0;i_dim<dim;++i_dim){
				point_id[i_dim]=(point_id[i_dim]+vector_gen(i_dim,0))%size_design;
			}
		}
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_factor_LeaveOneOut
	(type_int const & size_design_full,type_pt_fac const & vector_gen,OutIterType result_begin){
		LZ_DEF_func_check_traits((std::is_same<type_pt_fac,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		type_int const dim=vector_gen.rows();
		if(size_design_full<=(type_int)1 || dim<=(type_int)0 || vector_gen.cols()<=(type_int)0){
			return result_begin;
		}
		//checking ends.
		type_pt_fac point=type_pt_fac::Zero(dim,1),point_shift=type_pt_fac::Ones(dim,1);
		type_int i_dim,i_level;
		for(i_level=(type_int)1;i_level<size_design_full;++i_level){
			for(i_dim=0;i_dim<dim;++i_dim){
				point(i_dim,0)=(point(i_dim,0)+vector_gen(i_dim,0))%size_design_full;
			}
			*result_begin=point-point_shift;
			++result_begin;
		}
		return result_begin;
	}
	
	template<typename RealType,typename IntType>
	template<typename OutIterType>
	OutIterType
	Get_design_GoodLatticePoint<RealType,IntType>::Get_design_canon_LeaveOneOut
	(type_int const & size_design_full,type_pt_fac const & vector_gen,OutIterType result_begin){
		LZ_DEF_func_check_traits((std::is_same<type_pt_can,
			typename std::iterator_traits<OutIterType>::value_type>::value));
		type_int const dim=vector_gen.rows();
		if(size_design_full<=(type_int)1 || dim<=(type_int)0 || vector_gen.cols()<=(type_int)0){
			return result_begin;
		}
		//checking ends.
		type_int size_design=size_design_full-(type_int)1;
		type_pt_can point(dim,1);
		Type_ArrayTemp<type_real> crd(size_design);
		Type_ArrayTemp<type_int> point_id(dim,(type_int)0); //the point in the factor design;
		type_int i_dim,i_level;
		//calculating the 1-dim coordinates:
		for(i_level=0;i_level<size_design;++i_level){
			crd[i_level]=(type_real(i_level)+type_real(0.5))/type_real(size_design);
		}
		//transformation from the full factor design to leave-one-out canonical design:
		auto fun_Get_point=[&i_dim,&dim,&point,&crd,&point_id]()->type_pt_can&{
				for(i_dim=(type_int)0;i_dim<dim;++i_dim){
					point(i_dim,0)=crd[point_id[i_dim]-(type_int)1];
				}
				return point;
			};
		for(i_level=(type_int)0;i_level<size_design;++i_level){
			for(i_dim=0;i_dim<dim;++i_dim){
				point_id[i_dim]=(point_id[i_dim]+vector_gen(i_dim,0))%size_design_full;
			}
			*result_begin=fun_Get_point();
			++result_begin;
		}
		return result_begin;
	}
	
//}(class Get_design_GoodLatticePoint)
/****************************************************************************************************/

} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_design_uniform_CPP
#endif //#if LZ_DEF_LZ_H_stats_DOE_design_uniform!=202105L