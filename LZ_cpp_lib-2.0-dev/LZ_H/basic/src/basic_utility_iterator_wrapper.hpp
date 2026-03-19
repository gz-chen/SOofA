#if LZ_DEF_LZ_H_basic_basic!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_basic_src_basic_utility_iterator_wrapper
#define LZ_DEF_LZ_H_basic_src_basic_utility_iterator_wrapper 202105L

namespace Liuze{
/****************************************************************************************************/
//A wrapper for iterators which have specific std::iterator_traits<>::value_type{
#define LZ_DEF_LZ_H_basic_CODE_has_member_function(func_name,operation_pre,operation_post) \
	template<typename Tp,typename RetType> class has_member_function_##func_name{public: \
		template<typename Tp1,typename RetType1> static constexpr \
		typename std::enable_if<std::is_convertible< \
		decltype(operation_pre(std::declval<Tp1>())operation_post),RetType1>::value,bool>::type \
		test(int){return true;} \
		template<typename Tp1,typename RetType1> static constexpr bool test(double){return false;} \
		static bool constexpr value=test<Tp,RetType>(int(0)); \
	};
#define LZ_DEF_LZ_H_basic_CODE_INNER_has_member_function_true(class_type,ret_type,func_name) \
	typename std::enable_if<has_member_function_##func_name<class_type,ret_type>::value>::type* assval=NULL
#define LZ_DEF_LZ_H_basic_CODE_INNER_has_member_function_false(class_type,ret_type,func_name) \
	typename std::enable_if<!(has_member_function_##func_name<class_type,ret_type>::value)>::type* assval=NULL
#define LZ_DEF_LZ_H_basic_CODE_INNER_var_list(args...) args
#define LZ_DEF_LZ_H_basic_CODE_INNER_var_list_head(args...) args,
#define LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(args...) ,args
#define LZ_DEF_LZ_H_basic_CODE_INNER_define_invoke_member_function(func_name,var_list,ret_true,ret_false) \
	template<typename Tp,typename RetType> static inline RetType invoke_member_function_##func_name \
	(Tp & x,var_list LZ_DEF_LZ_H_basic_CODE_INNER_has_member_function_false(Tp,RetType,func_name)) \
	{return ret_false;} \
	template<typename Tp,typename RetType> static inline RetType invoke_member_function_##func_name \
	(Tp & x,var_list LZ_DEF_LZ_H_basic_CODE_INNER_has_member_function_true(Tp,RetType,func_name)) \
	{return ret_true;}
#define LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke\
(func_name,operation_pre,operation_post,var_list,ret_true,ret_false) \
	LZ_DEF_LZ_H_basic_CODE_has_member_function(func_name,operation_pre,operation_post) \
	LZ_DEF_LZ_H_basic_CODE_INNER_define_invoke_member_function \
	(func_name,LZ_DEF_LZ_H_basic_CODE_INNER_var_list(var_list),ret_true,ret_false)
#define LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function\
(func_name,class_type,ret_type,class_var,var_list) \
	type_this::invoke_member_function_##func_name<class_type,ret_type>(class_var var_list)
	
template<typename ValueType,typename IterCategory=std::input_iterator_tag,
	typename DifferenceType=Type_Int>
class Iterator_wrapper{
	LZ_DEF_func_check_traits((std::is_integral<DifferenceType>::value &&
		std::is_signed<DifferenceType>::value));
	public:
		//this class:
		typedef Iterator_wrapper<ValueType,IterCategory,DifferenceType> type_this;
		//typedefs for std::iterator_traits.
		typedef IterCategory iterator_category;
		typedef ValueType value_type;
		typedef DifferenceType difference_type;
		typedef ValueType* pointer;
		typedef ValueType& reference;
		//other typedefs:
		typedef typename std::add_const<value_type>::type value_type_const;
	protected:
		//member function invokers:
		//for input iterator:
		template<typename Tp,typename RetType>
		static inline RetType invoke_member_function_operator_rightarrow
		(Tp & x,
		typename std::enable_if<!(std::is_pointer<Tp>::value)>::type* assval=NULL){
			return x.operator ->();
		}
		template<typename Tp,typename RetType>
		static inline RetType invoke_member_function_operator_rightarrow
		(Tp & x,
		typename std::enable_if<std::is_pointer<Tp>::value>::type* assval=NULL){
			return x;
		}
		//for output iterator:
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_dereference,*,,
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(value_type & false_result,),*x,false_result)
		//for forward iterator:
		//for bidirectional iterator:
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_minusminus_prefix,--,,,--x,x)
		//for random_access iterator:
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_lt,,<(std::declval<Tp1>()),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(Tp const & y,),x<y,false)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_le,,<=(std::declval<Tp1>()),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(Tp const & y,),x<=y,false)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_pluseq_difftype,,+=difference_type(0),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),x+=step,x)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_minus_iterator,,-(std::declval<Tp1>()),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(Tp const & y,),x-y,difference_type(0))
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_index,,[difference_type(0)],
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),x[step],*x)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_plus_difftype,,+difference_type(0),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),x+step,x)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_difftype_plus,difference_type(0)+,,
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),step+x,x)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_minus_difftype,,-difference_type(0),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),x-step,x)
		LZ_DEF_LZ_H_basic_CODE_INNER_member_function_def_has_invoke
			(operator_minuseq_difftype,,-=difference_type(0),
			LZ_DEF_LZ_H_basic_CODE_INNER_var_list(difference_type const & step,),x-=step,x)
		
		class type_store_base{
			public:
				//Constructor:
				type_store_base()=default;
				//Destructor:
				virtual ~type_store_base()=default;
				
				//operator:
				//for input iterator:
				virtual Type_Bool operator ==(type_store_base const & iter_other) const=0;
				virtual value_type_const& operator *() const=0;
				virtual value_type_const* operator ->() const=0;
				virtual type_store_base& operator ++()=0;
				//for output iterator:
				virtual value_type& operator *()=0;
				virtual pointer operator ->()=0;
				//for forward iterator:
				//for bidirectional iterator:
				virtual type_store_base& operator --()=0;
				//for random_access iterator:
				virtual Type_Bool operator <(type_store_base const & iter_other) const=0;
				virtual Type_Bool operator >(type_store_base const & iter_other) const=0;
				virtual Type_Bool operator <=(type_store_base const & iter_other) const=0;
				virtual Type_Bool operator >=(type_store_base const & iter_other) const=0;
				virtual type_store_base& operator +=(difference_type const & step)=0;
				virtual difference_type operator -(type_store_base const & iter_other) const=0;
				virtual value_type_const& operator [](difference_type const & step) const=0;
				virtual value_type& operator [](difference_type const & step)=0;
				
				//other function:
				virtual type_store_base* clone()=0;
		};
		template<typename IterType>
		class type_store:public type_store_base{
			LZ_DEF_func_check_traits((std::is_convertible<
				typename std::remove_cv<typename std::iterator_traits<
				typename std::decay<IterType>::type>::value_type>::type,
				value_type>::value));
			public:
				typedef IterType type_iter;
				//Constructor:
				template<typename IterTypeType>
				type_store(IterTypeType && iter_source):iter(std::forward<IterTypeType>(iter_source)){}
				//Destructor:
				virtual ~type_store()=default;
				
				//operator:
				//for input iterator:
				virtual Type_Bool operator ==(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						iter==static_cast<type_store<type_iter> const &>(iter_other).iter : false;
				}
				virtual value_type_const& operator *() const {
					return LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_dereference,type_iter const,value_type_const&,
						const_cast<type_iter const &>(iter),
						LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(type_this::return_false));
				}
				virtual value_type_const* operator ->() const {
					return const_cast<value_type_const*>(
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_rightarrow,
						type_iter const,typename std::iterator_traits<type_iter>::value_type const *,iter,));
				}
				virtual type_store<type_iter>& operator ++(){
					++iter;
					return *this;
				}
				//for output iterator:
				virtual value_type& operator *(){
					return LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_dereference,type_iter,value_type&,
						iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(type_this::return_false));
				}
				virtual pointer operator ->(){
					return const_cast<pointer>(
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_rightarrow,
						type_iter,typename std::iterator_traits<type_iter>::pointer,iter,));
				}
				//for forward iterator:
				//for bidirectional iterator:
				virtual type_store<type_iter>& operator --(){
					LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_minusminus_prefix,type_iter,type_iter&,iter,);
					return *this;
				}
				//for random_access iterator:
				virtual Type_Bool operator <(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_lt,type_iter const,Type_Bool,
						iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(
						static_cast<type_store<type_iter> const &>(iter_other).iter)) :
						false;
				}
				virtual Type_Bool operator >(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_lt,type_iter const,Type_Bool,
						static_cast<type_store<type_iter> const &>(iter_other).iter,
						LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(iter)) :
						false;
				}
				virtual Type_Bool operator <=(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						(has_member_function_operator_le<type_iter const,Type_Bool>::value ?
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_le,type_iter const,Type_Bool,
						iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(
						static_cast<type_store<type_iter> const &>(iter_other).iter)) :
						(this->operator <(iter_other) || this->operator ==(iter_other))
						) : false;
				}
				virtual Type_Bool operator >=(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						(has_member_function_operator_le<type_iter const,Type_Bool>::value ?
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_le,type_iter const,Type_Bool,
						static_cast<type_store<type_iter> const &>(iter_other).iter,
						LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(iter)) :
						(this->operator >(iter_other) || this->operator ==(iter_other))
						) : false;
				}
				virtual type_store<type_iter>& operator +=(difference_type const & step){
					if(has_member_function_operator_pluseq_difftype<type_iter,type_iter&>::value){
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
							(operator_pluseq_difftype,type_iter,type_iter&,
							iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(step));
					} else {
						if(step>(difference_type)0){
							for(difference_type i=0;i<step;++i) ++iter;
						} else if(step<(difference_type)0){
							for(difference_type i=0;i>step;--i) --(*this);
						}
					}
					return *this;
				}
				virtual difference_type operator -(type_store_base const & iter_other) const {
					return typeid(iter_other)==typeid(*this) ?
						LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
						(operator_minus_iterator,type_iter const,difference_type,
						iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(
						static_cast<type_store<type_iter> const &>(iter_other).iter)) :
						difference_type(0);
				}
				virtual value_type_const& operator [](difference_type const & step) const {
					if(has_member_function_operator_index<type_iter const,value_type_const&>::value){
						return LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
							(operator_index,type_iter const,value_type_const&,
							const_cast<type_iter const &>(iter),
							LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(step));
					} else {
						type_store<type_iter> iter_new(iter);
						iter_new+=step;
						return *const_cast<type_store<type_iter> const &>(iter_new);
					}
				}
				virtual value_type& operator [](difference_type const & step){
					if(has_member_function_operator_index<type_iter,value_type&>::value){
						return LZ_DEF_LZ_H_basic_CODE_INNER_call_invoke_member_function
							(operator_index,type_iter,value_type&,
							iter,LZ_DEF_LZ_H_basic_CODE_INNER_var_list_tail(step));
					} else {
						type_store<type_iter> iter_new(iter);
						iter_new+=step;
						return *iter_new;
					}
				}
				
				//other function:
				virtual type_store_base* clone(){
					return static_cast<type_store_base*>(new type_store<type_iter>(iter));
				}
			protected:
				type_iter iter;
		};
	public:
		//Constructor:
		Iterator_wrapper():pt_iter(NULL){}
		template<typename IterType,
			typename=typename std::enable_if<!(std::is_same<type_this,
			typename std::decay<IterType>::type>::value)>::type>
		Iterator_wrapper(IterType && iter):
			pt_iter(static_cast<type_store_base*>
			(new type_store<typename std::decay<IterType>::type>(std::forward<IterType>(iter)))){}
		Iterator_wrapper(type_this const & iter):
			pt_iter(iter.pt_iter ? iter.pt_iter->clone() : NULL){}
		Iterator_wrapper(type_this && iter):pt_iter(iter.pt_iter){
			iter.pt_iter=NULL;
		}
		
		//Destructor:
		~Iterator_wrapper(){
			if(pt_iter) delete pt_iter;
		}
		
		//operator:
		//for input iterator:
		template<typename IterType,
			typename=typename std::enable_if<!(std::is_same<type_this,
			typename std::decay<IterType>::type>::value)>::type>
		type_this& operator =(IterType && iter){
			if(std::addressof(iter)==this) return *this;
			if(pt_iter) delete pt_iter;
			pt_iter=static_cast<type_store_base*>
				(new type_store<typename std::decay<IterType>::type>(std::forward<IterType>(iter)));
			return *this;
		}
		type_this& operator =(type_this const & iter){
			if(this==&iter) return *this;
			if(pt_iter) delete pt_iter;
			pt_iter= iter.pt_iter ? iter.pt_iter->clone() : NULL;
			return *this;
		}
		type_this& operator =(type_this && iter){
			if(pt_iter) delete pt_iter;
			pt_iter=iter.pt_iter;
			iter.pt_iter=NULL;
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return (*pt_iter)==(*(iter.pt_iter));
		}
		value_type_const& operator *() const {
			return const_cast<value_type_const&>(*static_cast<type_store_base const &>(*pt_iter));
		}
		value_type_const* operator ->() const {
			return const_cast<value_type_const*>
				(static_cast<type_store_base const *>(pt_iter)->operator ->());
		}
		type_this& operator ++(){
			++(*pt_iter);
			return *this;
		}
		template<typename ValueType_1,typename IterCategory_1,typename DifferenceType_1>
		friend void std::swap
		(Iterator_wrapper<ValueType_1,IterCategory_1,DifferenceType_1> & iter0,
		Iterator_wrapper<ValueType_1,IterCategory_1,DifferenceType_1> & iter1);
		//for output iterator:
		value_type& operator *(){
			return *(*pt_iter);
		}
		pointer operator ->(){
			return pt_iter->operator ->();
		}
		//for forward iterator:
		//for bidirectional iterator:
		type_this& operator --(){
			--(*pt_iter);
			return *this;
		}
		//for random_access iterator:
		Type_Bool operator <(type_this const & iter) const {
			return (*pt_iter)<(*(iter.pt_iter));
		}
		Type_Bool operator >(type_this const & iter) const {
			return (*pt_iter)>(*(iter.pt_iter));
		}
		Type_Bool operator <=(type_this const & iter) const {
			return (*pt_iter)<=(*(iter.pt_iter));
		}
		Type_Bool operator >=(type_this const & iter) const {
			return (*pt_iter)>=(*(iter.pt_iter));
		}
		type_this& operator +=(difference_type const & step){
			(*pt_iter)+=step;
			return *this;
		}
		type_this& operator -=(difference_type const & step){
			return (*this)+=(-step);
		}
		friend type_this operator +(type_this iter,difference_type const & step){
			return iter+=step;
		}
		friend type_this operator +(difference_type const & step,type_this iter){
			return iter+=step;
		}
		friend type_this operator -(type_this iter,difference_type const & step){
			return iter-=step;
		}
		friend difference_type operator -(type_this const & iter0,type_this const & iter1){
			return (*(iter0.pt_iter))-(*(iter1.pt_iter));
		}
		value_type_const& operator [](difference_type const & step) const {
			return const_cast<value_type_const&>(static_cast<type_store_base const &>(*pt_iter)[step]);
		}
		value_type& operator [](difference_type const & step){
			return (*pt_iter)[step];
		}
	protected:
		type_store_base* pt_iter; //the pointer to the stored iterator;
	private:
		static value_type return_false;
}; //class Iterator_wrapper;
//static values in class Iterator_wrapper{
template<typename ValueType,typename IterCategory,typename DifferenceType>
typename Iterator_wrapper<ValueType,IterCategory,DifferenceType>::value_type
Iterator_wrapper<ValueType,IterCategory,DifferenceType>::return_false;
//}(static values in class Iterator_wrapper)
//}(A wrapper for iterators which have specific std::iterator_traits<>::value_type)
/****************************************************************************************************/
} //namespace Liuze;

namespace std{
template<typename ValueType,typename IterCategory,typename DifferenceType>
void swap
(Liuze::Iterator_wrapper<ValueType,IterCategory,DifferenceType> & iter0,
Liuze::Iterator_wrapper<ValueType,IterCategory,DifferenceType> & iter1){
	swap(iter0.pt_iter,iter1.pt_iter);
}
} //namespace std;

#endif //#ifndef LZ_DEF_LZ_H_basic_src_basic_utility_iterator_wrapper
#endif //#if LZ_DEF_LZ_H_basic_basic!=202105L