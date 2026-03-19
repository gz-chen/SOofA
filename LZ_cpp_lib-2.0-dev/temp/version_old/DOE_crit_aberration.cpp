#if LZ_DEF_LZ_H_stats_DOE_crit_aberration!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_aberration_CPP
#define LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_aberration_CPP 202105L

namespace Liuze{
namespace stats{
namespace crit{

/****************************************************************************************************/
//class WordlengthPattern{
//public:
	//Constructor:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	WordlengthPattern<RealType,IntType,TraitsCheck>::WordlengthPattern
	(type_int const & num_level):
	m_num_level(num_level>=type_int(1) ? num_level : type_int(2)){
		typedef Type_MatTemp<type_real> type_matrix;
		m_tab_powval.Set_enlarge(m_num_level,m_num_level);
		m_tab_poly.resize(m_num_level);
		m_tab_polyval.resize(m_num_level);
		type_int i_poly,i_level,i_pow,j_poly;
		m_tab_poly[0].resize(1);
		m_tab_poly[0][0]=type_real(1);
		m_tab_polyval[0].assign(m_num_level,type_real(1));
		type_matrix poly_LE_coef,poly_LE_val,poly_LE_sol;
		type_real sqnorm;
		for(i_poly=1;i_poly<m_num_level;++i_poly){
			poly_LE_coef=type_matrix::Zero(i_poly,i_poly);
			poly_LE_val=type_matrix::Zero(i_poly,1);
			for(j_poly=0;j_poly<i_poly;++j_poly){
				for(i_level=0;i_level<m_num_level;++i_level){
					for(i_pow=0;i_pow<i_poly;++i_pow){
						poly_LE_coef(j_poly,i_pow)+=
							m_tab_powval[i_level][i_pow]*m_tab_polyval[j_poly][i_level];
					}
					poly_LE_val(j_poly,0)-=
						m_tab_powval[i_level][i_poly]*m_tab_polyval[j_poly][i_level];
				}
			}
			poly_LE_sol=poly_LE_coef.colPivHouseholderQr().solve(poly_LE_val);
			m_tab_poly[i_poly].resize(i_poly+type_int(1));
			for(i_pow=0;i_pow<i_poly;++i_pow) m_tab_poly[i_poly][i_pow]=poly_LE_sol(i_pow,0);
			m_tab_poly[i_poly][i_poly]=type_real(1);
			m_tab_polyval[i_poly].assign(m_num_level,m_tab_poly[i_poly][0]);
			sqnorm=type_real(0);
			for(i_level=0;i_level<m_num_level;++i_level){
				for(i_pow=1;i_pow<=i_poly;++i_pow){
					m_tab_polyval[i_poly][i_level]+=
						m_tab_poly[i_poly][i_pow]*m_tab_powval[i_level][i_pow];
				}
				sqnorm+=m_tab_polyval[i_poly][i_level]*m_tab_polyval[i_poly][i_level];
			}
			m_tab_poly[i_poly][i_poly]=type_real(sqrt(type_real(m_num_level)/sqnorm));
			for(i_pow=0;i_pow<i_poly;++i_pow)
				m_tab_poly[i_poly][i_pow]*=m_tab_poly[i_poly][i_poly];
			for(i_level=0;i_level<m_num_level;++i_level)
				m_tab_polyval[i_poly][i_level]*=m_tab_poly[i_poly][i_poly];
		}
	}
	
	//other function:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,typename>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_res
	WordlengthPattern<RealType,IntType,TraitsCheck>::beta
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
		typedef std::complex<type_real> type_comp;
		typedef typename std::make_signed<type_int>::type type_sint;
		type_int dim=des_begin->rows();
		if(dim==type_int(0) || m_num_level<type_int(2)) return type_res(1,type_real(1));
		type_int maxlen=dim*(m_num_level-type_int(1));
		type_int len_res=maxlen+type_int(1);
		type_res res(len_res);
		res[0]=type_real(1);
		type_int i_root,i_len;
		type_real tr0=type_real(LZ_DEF_const_math_2pi)/type_real(maxlen);
		Type_ArrayTemp<type_comp> list_root(maxlen);
		for(i_root=0;i_root<maxlen;++i_root){
			list_root[i_root]=type_comp
				(type_real(cos(type_real(i_root)*tr0)),type_real(sin(type_real(i_root)*tr0)));
		}
		Type_ArrayTemp<type_comp> enum_in(m_num_level-type_int(1));
		Type_ArrayTemp<type_comp> enum_val(maxlen);
		for(i_len=0;i_len<maxlen;++i_len){
			for(i_root=1;i_root<m_num_level;++i_root){
				enum_in[i_root-type_int(1)]=list_root[((i_len+type_int(1))*i_root)%maxlen];
			}
			enum_val[i_len]=
				this->enumerator(des_begin,des_end,weight_begin,weight_end,enum_in.begin())
				-type_comp(1);
		}
		tr0=type_real(maxlen);
		type_comp res0;
		for(i_len=1;i_len<=maxlen;++i_len){
			res0=0;
			for(i_root=0;i_root<maxlen;++i_root){
				res0+=enum_val[i_root]*
					list_root[Liuze::math::mod(-type_sint(i_len*(i_root+type_int(1))),type_sint(maxlen))];
			}
			res[i_len]=res0.real()/tr0;
		}
		return res;
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,
		typename ForwardIterType_enumval,typename>
	typename std::decay<typename std::iterator_traits<ForwardIterType_enumval>::value_type>::type
	WordlengthPattern<RealType,IntType,TraitsCheck>::enumerator
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	ForwardIterType_enumval const & enumval_begin) const {
		typedef
			typename std::decay<typename std::iterator_traits<ForwardIterType_enumval>::value_type>::type
			type_val;
		type_int dim=des_begin->rows();
		if(dim==type_int(0)) return Liuze::math::one(*enumval_begin);
		type_int i_level,i0_level,i1_level,i_dim;
		type_val cv1=Liuze::math::one(*enumval_begin);
		ForwardIterType_enumval iter_enumval;
		Type_ArrayTemp<type_val> ContrastSimil(m_num_level*m_num_level);
		for(i0_level=0;i0_level<m_num_level;++i0_level){
			for(i1_level=0;i1_level<i0_level;++i1_level){
				ContrastSimil[i0_level*m_num_level+i1_level]=
					ContrastSimil[i1_level*m_num_level+i0_level];
			}
			for(;i1_level<m_num_level;++i1_level){
				i_dim=i0_level*m_num_level+i1_level;
				iter_enumval=enumval_begin;
				ContrastSimil[i_dim]=cv1;
				for(i_level=1;i_level<m_num_level;++i_level){
					ContrastSimil[i_dim]+=type_val(
						m_tab_polyval[i_level][i0_level]*m_tab_polyval[i_level][i1_level]*
						(*iter_enumval));
					++iter_enumval;
				}
			}
		}
		type_val res_1;
		type_val res_ndiag=Liuze::math::zero(*enumval_begin),res_diag=res_ndiag;
		type_real wt_1;
		type_real wtsum_ndiag=0,wtsum_diag=0;
		ForwardIterType_des iter0_des=des_begin,iter1_des;
		ForwardIterType_weight iter0_wt=weight_begin,iter1_wt;
		while(iter0_des!=des_end && iter0_wt!=weight_end){
			iter1_des=des_begin;
			iter1_wt=weight_begin;
			while(iter1_des!=iter0_des){
				res_1=cv1;
				for(i_dim=0;i_dim<dim;++i_dim){
					res_1*=ContrastSimil[(*iter0_des)(i_dim,0)*m_num_level+(*iter1_des)(i_dim,0)];
				}
				wt_1=(*iter0_wt)*(*iter1_wt);
				res_ndiag+=wt_1*res_1;
				wtsum_ndiag+=wt_1;
				++iter1_des;
				++iter1_wt;
			}
			{
				res_1=cv1;
				for(i_dim=0;i_dim<dim;++i_dim){
					res_1*=ContrastSimil[(*iter0_des)(i_dim,0)*m_num_level+(*iter0_des)(i_dim,0)];
				}
				wt_1=(*iter0_wt)*(*iter0_wt);
				res_diag+=wt_1*res_1;
				wtsum_diag+=wt_1;
			}
			++iter0_des;
			++iter0_wt;
		}
		return type_val((type_real(1)/(wtsum_ndiag+wtsum_ndiag+wtsum_diag))*
			(res_ndiag+res_ndiag+res_diag));
	}
	
//protected:
	//types:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_powval{
		public:
			//types:
			typedef type_tab_powval type_this;
			
			//Constructor:
			type_tab_powval(void):
			m_val(0){
			}
			type_tab_powval(type_int const & num_base,type_int const & num_power):
			m_val(num_base>type_int(0) ? num_base : type_int(0),
			Type_ArrayTemp<type_real>(num_power>type_int(0) ? num_power : type_int(0))){
				if(num_base>type_int(0) && num_power>type_int(0)){
					type_int i_base,i_pow;
					for(i_base=0;i_base<num_base;++i_base){
						m_val[i_base][0]=type_real(1);
						for(i_pow=1;i_pow<num_power;++i_pow){
							m_val[i_base][i_pow]=m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
						}
					}
				}
			}
			
			//Destructor:
			~type_tab_powval(void)=default;
			
			//operator:
			type_real operator ()(type_int const & base,type_int const & power) const {
				return (base>=type_int(0) && base<Get_num_base() &&
					power>=type_int(0) && power<Get_num_power()) ?
					m_val[base][power] : std::numeric_limits<type_real>::quiet_NaN();
			}
			Type_ArrayTemp<type_real> const & operator [](type_int const & base) const {
				return m_val[base];
			}
			
			//other function:
			type_int Get_num_base(void) const {
				return type_int(m_val.size());
			}
			type_int Get_num_power(void) const {
				return Get_num_base()==type_int(0) ? type_int(0) : type_int(m_val[0].size());
			}
			void Set_enlarge(type_int const & num_base,type_int const & num_power){
				type_int n_base=Get_num_base(),n_pow=Get_num_power();
				if(num_base<=n_base && num_power<=n_pow) return;
				type_int nn_base=std::max(n_base,num_base);
				type_int nn_pow=std::max(n_pow,num_power);
				type_int i_base,i_pow;
				if(nn_pow>n_pow){
					for(i_base=0;i_base<n_base;++i_base){
						m_val[i_base].resize(nn_pow);
						i_pow= n_pow==type_int(0) ?
							(m_val[i_base][0]=type_real(1),type_int(1)) : n_pow;
						for(;i_pow<nn_pow;++i_pow){
							m_val[i_base][i_pow]=m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
						}
					}
				}
				if(nn_base>n_base){
					m_val.resize(nn_base);
					for(i_base=n_base;i_base<nn_base;++i_base){
						m_val[i_base].resize(nn_pow);
						m_val[i_base][0]=type_real(1);
						for(i_pow=1;i_pow<nn_pow;++i_pow){
							m_val[i_base][i_pow]=m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
						}
					}
				}
				return;
			}
			type_real Cal(type_int const & base,type_int const & power){
				if(base<type_int(0) || power<type_int(0)){
					return std::numeric_limits<type_real>::quiet_NaN();
				}
				if(base<Get_num_base() && power<Get_num_power()) return m_val[base][power];
				this->Set_enlarge(base+type_int(1),power+type_int(1));
				return m_val[base][power];
			}
		protected:
			Type_ArrayTemp<Type_ArrayTemp<type_real> > m_val;
	}; //class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_powval
	
	//static member values:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_powval
	WordlengthPattern<RealType,IntType,TraitsCheck>::m_tab_powval;
	
//}(class WordlengthPattern)
/****************************************************************************************************/
} //namespace crit;
} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_aberration_CPP
#endif //#if LZ_DEF_LZ_H_stats_DOE_crit_aberration!=202105L