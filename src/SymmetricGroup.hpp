#ifndef LZ_DEF_hpp_SymmetricGroup
#define LZ_DEF_hpp_SymmetricGroup

//`#ifndef LZ_DEF_LZ_H_math_minimisation
//`#define LZ_DEF_LZ_H_math_minimisation 202105L

#include<LZ_H/math/basic.hpp>
//`#include"LZ_H/math/differentiation.hpp"

//`#ifdef LZ_DEF_extLIB_Eigen

namespace Liuze{
namespace math{

/****************************************************************************************************/
//symmetric group on finite set {0,...,n-1}{
template<typename IntType=Type_UInt,typename RealType=Type_Real,
	typename=typename std::enable_if<std::is_integral<IntType>::value &&
	std::is_arithmetic<RealType>::value>::type>
class FiniteSymmetricGroup{
	public:
		typedef FiniteSymmetricGroup<IntType,RealType,void> type_this;
		typedef IntType type_int;
		typedef RealType type_real;
		
		//Constructor:
		FiniteSymmetricGroup():ord_base(0){
		}
		explicit FiniteSymmetricGroup(type_int const & base_order):
		ord_base(base_order>=type_int(0) ? base_order : type_int(0)){
		}
		FiniteSymmetricGroup(type_this const & group):
		ord_base(group.ord_base){
		}
		
		//Destructor:
		~FiniteSymmetricGroup(void)=default;
		
		//operator:
		type_this & operator =(type_this const & group){
			ord_base=group.ord_base;
			return *this;
		}
		type_this & operator =(type_int const & base_order){
			if(base_order>=type_int(0)) ord_base=base_order;
			return *this;
		}
		
		template<typename iter_permutation>
		static inline type_int act
		(iter_permutation const & permutation,type_int const & x){
			return permutation[x];
		}
		template<typename vec_permutation=Type_ArrayTemp<type_int>,typename iter_permutation>
		static vec_permutation inverse
		(iter_permutation const & permutation,type_int const & base_order){
			if(base_order<=type_int(1)) return permutation;
			vec_permutation res(base_order);
			for(type_int i=0;i<base_order;++i) res[permutation[i]]=i;
			return res;
		}
		template<typename vec_permutation>
		static inline vec_permutation inverse
		(vec_permutation const & permutation){
			return type_this::template inverse<vec_permutation,vec_permutation>
				(permutation,permutation.size());
		}
		template<typename vec_permutation=Type_ArrayTemp<type_int>,typename iter_permutation>
		static vec_permutation compose
		(iter_permutation const & perm_left,iter_permutation const & perm_right,
		type_int const & base_order){
			if(base_order<=type_int(1)) return perm_left;
			vec_permutation res(base_order);
			for(type_int i=0;i<base_order;++i) res[i]=perm_left[perm_right[i]];
			return res;
		}
		template<typename vec_permutation>
		static inline vec_permutation compose
		(vec_permutation const & perm_left,vec_permutation const & perm_right){
			return type_this::template compose<vec_permutation,vec_permutation>
				(perm_left,perm_right,perm_right.size());
		}
		
		template<typename Scalar=typename type_this::type_real,typename iter_permutation,
			typename=typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>
		static Type_MatTemp<Scalar> matrix_permutation_row
		(iter_permutation const & permutation,type_int const & base_order){
			typedef Scalar type_scalar;
			typedef Type_MatTemp<Scalar> type_mat_perm;
			if(base_order<=type_int(0)) return type_mat_perm(0,0);
			type_mat_perm res=type_mat_perm::Zero(base_order,base_order);
			type_int i_comp;
			for(i_comp=0;i_comp<base_order;++i_comp){
				res(i_comp,permutation[i_comp])=type_scalar(1);
			}
			return res;
		}
		template<typename Scalar=typename type_this::type_real,typename vec_permutation,
			typename=typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>
		static inline Type_MatTemp<Scalar> matrix_permutation_row
		(vec_permutation const & permutation){
			return type_this::template matrix_permutation_row<Scalar,vec_permutation>
				(permutation,permutation.size());
		}
		template<typename Scalar=typename type_this::type_real,typename iter_permutation,
			typename=typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>
		static inline Type_MatTemp<Scalar> matrix_permutation_col
		(iter_permutation const & permutation,type_int const & base_order){
			return type_this::template matrix_permutation_row<Scalar,iter_permutation>
				(permutation,base_order).transpose();
		}
		template<typename Scalar=typename type_this::type_real,typename vec_permutation,
			typename=typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>
		static inline Type_MatTemp<Scalar> matrix_permutation_col
		(vec_permutation const & permutation){
			return type_this::template matrix_permutation_row<Scalar,vec_permutation>
				(permutation,permutation.size()).transpose();
		}
		
		template<typename iter_permutation>
		static type_real weight_Cayley
		(iter_permutation const & permutation,type_int const & base_order){
			if(base_order<=type_int(1)) return type_real(0);
			Type_ArrayTemp<bool> flag(base_order,false);
			type_real res=0;
			type_int i0,i,s0=1,s;
			for(i0=0;i0<base_order && s0<base_order;++i0){
				if(flag[i0]) continue;
				s=0;
				i=i0;
				do{
					flag[i]=true;
					i=permutation[i];
					++s;
				}while(i!=i0);
				s0+=s;
				if(s>type_int(1)) res+=type_real(s-type_int(1));
			}
			return res;
		}
		template<typename iter_permutation>
		static inline type_real distance_Cayley
		(iter_permutation const & perm0,iter_permutation const & perm1,
		type_int const & base_order){
			return type_this::weight_Cayley(
				type_this::compose(type_this::inverse(perm0,base_order),perm1,base_order),
				base_order);
		}
		template<typename vec_permutation>
		static inline type_real distance_Cayley
		(vec_permutation const & perm0,vec_permutation const & perm1){
			return type_this::weight_Cayley(
				type_this::compose(type_this::inverse(perm0),perm1),
				perm1.size());
		}
		template<typename iter_permutation>
		static type_real weight_Hamming
		(iter_permutation const & permutation,type_int const & base_order){
			if(base_order<=type_int(1)) return type_real(0);
			type_real res=0;
			for(type_int i=0;i<base_order;++i){
				if(permutation[i]!=i) ++res;
			}
			return res;
		}
		template<typename iter_permutation>
		static inline type_real distance_Hamming
		(iter_permutation const & perm0,iter_permutation const & perm1,
		type_int const & base_order){
			return type_this::weight_Hamming(
				type_this::compose(type_this::inverse(perm0,base_order),perm1,base_order),
				base_order);
		}
		template<typename vec_permutation>
		static inline type_real distance_Hamming
		(vec_permutation const & perm0,vec_permutation const & perm1){
			return type_this::weight_Hamming(
				type_this::compose(type_this::inverse(perm0),perm1),
				perm1.size());
		}
	protected:
		type_int ord_base; //it is the symmetric group on {0,...,ord_base-1};
}; //class FiniteSymmetricGroup;
//}(symmetric group on finite set {0,...,n-1})
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

//implementation:
//`#include"LZ_H/math/src/minimisation.cpp"

//`#endif //#ifdef LZ_DEF_extLIB_Eigen
//`#endif //#ifndef LZ_DEF_LZ_H_math_minimisation

#endif //#ifndef LZ_DEF_hpp_SymmetricGroup