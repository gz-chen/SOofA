#if LZ_DEF_LZ_H_basic_basic!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_basic_src_basic_utility
#define LZ_DEF_LZ_H_basic_src_basic_utility 202105L

namespace Liuze{
/****************************************************************************************************/
//default operator relying on related operator{
template<typename Tp>
inline Type_Bool operator !=(Tp const & x0,Tp const & x1){
	return !(x0==x1);
}
template<typename Tp>
inline Tp operator ++(Tp & x,int){
	Tp res=x;
	++x;
	return res;
}
template<typename Tp>
inline Tp operator --(Tp & x,int){
	Tp res=x;
	--x;
	return res;
}
template<typename Tp>
inline Tp & operator +=(Tp & x0,Tp const & x1){
	return x0=x0+x1;
}
template<typename Tp>
inline Tp & operator -=(Tp & x0,Tp const & x1){
	return x0=x0-x1;
}
template<typename Tp>
inline Tp & operator *=(Tp & x0,Tp const & x1){
	return x0=x0*x1;
}
template<typename Tp>
inline Tp & operator /=(Tp & x0,Tp const & x1){
	return x0=x0/x1;
}
template<typename Tp>
inline Tp & operator %=(Tp & x0,Tp const & x1){
	return x0=x0%x1;
}
template<typename Tp>
inline Tp & operator &=(Tp & x0,Tp const & x1){
	return x0=x0&x1;
}
template<typename Tp>
inline Tp & operator |=(Tp & x0,Tp const & x1){
	return x0=x0|x1;
}
template<typename Tp>
inline Tp & operator ^=(Tp & x0,Tp const & x1){
	return x0=x0^x1;
}
template<typename Tp>
inline Tp & operator <<=(Tp & x0,Tp const & x1){
	return x0=x0<<x1;
}
template<typename Tp>
inline Tp & operator >>=(Tp & x0,Tp const & x1){
	return x0=x0>>x1;
}
//}(default operator relying on related operator)
/****************************************************************************************************/
//The iterator for a finite or infinite sequence of a common constant value{
template<typename Tp,typename UIntType=Type_Size>
class Iterator_constval{
	LZ_DEF_func_check_traits((std::is_unsigned_integral<UIntType>::value));
	public:
		//Type: this class:
		typedef Iterator_constval<Tp> type_this;
		//Type: the value type:
		typedef Tp type_val;
		//Type: the size:
		typedef UIntType type_size;
		//Type: its iterator:
		typedef type_this iterator;
		//Type: for iterator:
		typedef std::iterator<std::forward_iterator_tag,type_val> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//values:
		static Type_Stat constexpr val_stat_deref=0; //de-reference-able;
		static Type_Stat constexpr val_stat_end=1; //at after the last entry;
		static Type_Stat constexpr val_stat_rend=2; //at before the first entry;
		static Type_Stat constexpr val_stat_deref_inf=10; //infinite length;
		static Type_Stat constexpr val_stat_end_inf=11; //virtual end of infinite length;
		static Type_Stat constexpr val_stat_rend_inf=12; //virtual reversed-end of infinite length;
		
		//Constructor:
		Iterator_constval():
		val(),stat(val_stat_end_inf){
		}
		Iterator_constval(type_this const & iter):
		val(iter.val),stat(iter.stat),n_size(iter.n_size),i_size(iter.i_size){
		}
		Iterator_constval(type_this && iter):
		val(std::forward<type_val>(iter.val)),stat(iter.stat),
		n_size(iter.n_size),i_size(iter.i_size){
		}
		explicit Iterator_constval(type_val const & value):
		val(value),stat(val_stat_deref_inf){
		}
		explicit Iterator_constval(type_val && value):
		val(std::forward<type_val>(value)),stat(val_stat_deref_inf){
		}
		Iterator_constval(type_val const & value,type_size const & num_size):
		val(value),n_size(num_size),i_size(0),
		stat(num_size>(type_size)0 ? val_stat_deref : val_stat_end){
		}
		Iterator_constval(type_val && value,type_size const & num_size):
		val(std::forward<type_val>(value)),n_size(num_size),i_size(0),
		stat(num_size>(type_size)0 ? val_stat_deref : val_stat_end){
		}
		
		//Destructor:
		~Iterator_constval()=default;
		
		//operator:
		type_this& operator =(type_this const & iter){
			if(&iter==this) return *this;
			val=iter.val;
			n_size=iter.n_size;
			i_size=iter.i_size;
			stat=iter.stat;
			return *this;
		}
		type_this& operator =(type_this && iter){
			val=std::forward<type_val>(iter.val);
			n_size=iter.n_size;
			i_size=iter.i_size;
			stat=iter.stat;
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return (stat==iter.stat && val==iter.val) ?
				(this->Is_finite() ? (n_size==iter.n_size && i_size==iter.i_size) : true) :
				false;
		}
		Type_Bool operator !=(type_this const & iter) const {
			return !((*this)==iter);
		}
		type_val& operator *(){
			return val;
		}
		type_val const & operator *() const {
			return val;
		}
		type_val* operator ->(){
			return std::addressof(val);
		}
		type_val const * operator ->() const {
			return std::addressof(val);
		}
		type_this& operator ++(){
			if(this->Is_finite()){
				switch(stat){
					case val_stat_deref:
						++i_size;
						if(i_size==n_size) stat=val_stat_end;
						break;
					case val_stat_end:
						break;
					case val_stat_rend:
						stat=val_stat_rend;
						break;
				}
			}
			return *this;
		}
		type_this operator ++(int){
			type_this res=*this;
			++(*this);
			return res;
		}
		
		//Is:
		Type_Bool Is_finite() const {
			return (stat==val_stat_deref || stat==val_stat_end || stat==val_stat_rend);
		}
		Type_Bool Is_infinite() const {
			return !(this->Is_finite());
		}
		Type_Bool Is_deref() const {
			return (stat==val_stat_deref || stat==val_stat_deref_inf);
		}
		Type_Bool Is_end() const {
			return (stat==val_stat_end || stat==val_stat_end_inf);
		}
		Type_Bool Is_rend() const {
			return (stat==val_stat_rend || stat==val_stat_rend_inf);
		}
		
		//begin and end:
		type_this begin() const {
			return this->Is_finite() ? type_this(val,n_size) : type_this(val);
		}
		type_this end() const {
			type_this res(*this);
			res.stat=val_stat_end;
			if(this->Is_finite()) res.i_size=n_size;
			return res;
		}
		
		//Get:
		type_size Get_size() const {
			return n_size;
		}
		type_size size() const {
			return n_size;
		}
		type_val Get_value() const {
			return val;
		}
		
		//Set:
		void resize(){
			stat=val_stat_deref_inf;
		}
		void resize(type_size const & num_size){
			n_size=num_size;
			i_size=(type_size)0;
			stat= num_size==(type_size)0 ? val_stat_end : val_stat_deref;
		}
		void assign(type_val const & value){
			val=value;
			stat=val_stat_deref_inf;
		}
		void assign(type_val && value){
			val=std::forward<type_val>(value);
			stat=val_stat_deref_inf;
		}
		void assign(type_val const & value,type_size const & num_size){
			val=value;
			n_size=num_size;
			i_size=(type_size)0;
			stat= num_size>(type_size)0 ? val_stat_deref : val_stat_end;
		}
		void assign(type_val && value,type_size const & num_size){
			val=std::forward<type_val>(value);
			n_size=num_size;
			i_size=(type_size)0;
			stat= num_size>(type_size)0 ? val_stat_deref : val_stat_end;
		}
	protected:
		type_val val; //the value;
		Type_Stat stat; //the state;
		type_size n_size,i_size; //size and id of size;
}; //class Iterator_constval;
//}(The iterator for a finite or infinite sequence of a common constant value)
/****************************************************************************************************/
} //namespace Liuze;

//A multi-way time recorder:
#include"LZ_H/basic/src/basic_utility_time_recorder.hpp"
//A wrapper for iterators which have specific std::iterator_traits<>::value_type:
#include"LZ_H/basic/src/basic_utility_iterator_wrapper.hpp"

#endif //#ifndef LZ_DEF_LZ_H_basic_src_basic_utility
#endif //#if LZ_DEF_LZ_H_basic_basic!=202105L