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
	//other function:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,typename>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_res
	WordlengthPattern<RealType,IntType,TraitsCheck>::Alphabeta
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	Type_ArrayTemp<type_int> const & num_level_qual,
	Type_ArrayTemp<type_int> const & num_level_quan) const {
		typedef std::complex<type_real> type_comp;
		typedef typename std::make_signed<type_int>::type type_sint;
		type_int dim_qual=num_level_qual.size();
		type_int dim_quan=num_level_quan.size();
		type_int i_dim;
		type_int maxlen=0;
		for(i_dim=0;i_dim<dim_qual;++i_dim){
			if(num_level_qual[i_dim]>type_int(1)) ++maxlen;
		}
		for(i_dim=0;i_dim<dim_quan;++i_dim){
			maxlen+=num_level_quan[i_dim]-type_int(1);
		}
		if(maxlen==type_int(0)) return type_res(1,type_real(1));
		type_int len_res=maxlen+type_int(1);
		type_res res(len_res);
		res[0]=type_real(1);
		type_int i_len;
		type_real tr0=type_real(LZ_DEF_const_math_2pi)/type_real(maxlen);
		Type_ArrayTemp<type_comp> list_root(maxlen);
		for(i_len=0;i_len<maxlen;++i_len){
			list_root[i_len]=type_comp
				(type_real(cos(type_real(i_len)*tr0)),type_real(sin(type_real(i_len)*tr0)));
		}
		Type_ArrayTemp<type_comp> enum_val(maxlen);
		for(i_len=0;i_len<maxlen;++i_len){
			enum_val[i_len]=
				this->enumerator_Alphabeta(des_begin,des_end,weight_begin,weight_end,
					num_level_qual,num_level_quan,list_root[(i_len+type_int(1))%maxlen])
				-type_comp(1);
		}
		tr0=type_real(maxlen);
		type_comp res0;
		for(i_len=1;i_len<=maxlen;++i_len){
			res0=0;
			for(i_dim=0;i_dim<maxlen;++i_dim){
				res0+=enum_val[i_dim]*
					list_root[Liuze::math::mod(-type_sint(i_len*(i_dim+type_int(1))),type_sint(maxlen))];
			}
			res[i_len]=res0.real()/tr0;
		}
		return res;
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,
		typename ValType,typename>
	ValType
	WordlengthPattern<RealType,IntType,TraitsCheck>::enumerator_Alphabeta
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	Type_ArrayTemp<type_int> const & num_level_qual,Type_ArrayTemp<type_int> const & num_level_quan,
	ValType const & val) const {
		typedef ValType type_val;
		type_int dim_qual=num_level_qual.size();
		type_int dim_quan=num_level_quan.size();
		type_int dim=dim_qual+dim_quan;
		if(dim==type_int(0)) return Liuze::math::one(val);
		type_int i_dim,i_level,n_level,i_n_level,n_n_level;
		type_val cv1=Liuze::math::one(val);
		//contrast similarity for qualitative factors{
		Type_ArrayTemp<type_int> id_fac_qual(dim_qual);
		Type_ArrayTemp<type_val> ContrastSimil_qual(0);
		if(dim_qual>type_int(0)){
			for(i_dim=0;i_dim<dim_qual;++i_dim) id_fac_qual[i_dim]=i_dim;
			std::sort(id_fac_qual.begin(),id_fac_qual.end(),
				[&num_level_qual](type_int const & i0,type_int const & i1)->bool
				{return num_level_qual[i0]<num_level_qual[i1];});
			n_n_level=n_level=0;
			for(i_dim=0;i_dim<dim_qual;++i_dim){
				if(num_level_qual[id_fac_qual[i_dim]]>n_level){
					++n_n_level;
					n_level=num_level_qual[id_fac_qual[i_dim]];
				}
			}
			if(n_n_level==type_int(0)) return std::numeric_limits<type_val>::quiet_NaN();
			ContrastSimil_qual.resize(n_n_level*type_int(2),cv1);
			i_dim=i_n_level=n_level=0;
			while(i_n_level<n_n_level){
				while(num_level_qual[id_fac_qual[i_dim]]==n_level) ++i_dim;
				n_level=num_level_qual[id_fac_qual[i_dim]];
				ContrastSimil_qual[i_n_level*type_int(2)]-=val;
				ContrastSimil_qual[i_n_level*type_int(2)+type_int(1)]+=
					type_val(n_level-type_int(1))*val;
				++i_n_level;
			}
		}
		//}
		//contrast similarity for quantitative factors{
		Type_ArrayTemp<type_int> id_fac_quan(dim_quan);
		Type_ArrayTemp<Type_ArrayTemp<type_val> > ContrastSimil_quan(0);
		if(dim_quan>type_int(0)){
			for(i_dim=0;i_dim<dim_quan;++i_dim) id_fac_quan[i_dim]=i_dim;
			std::sort(id_fac_quan.begin(),id_fac_quan.end(),
				[&num_level_quan](type_int const & i0,type_int const & i1)->bool
				{return num_level_quan[i0]<num_level_quan[i1];});
			n_n_level=n_level=0;
			for(i_dim=0;i_dim<dim_quan;++i_dim){
				if(num_level_quan[id_fac_quan[i_dim]]>n_level){
					++n_n_level;
					n_level=num_level_quan[id_fac_quan[i_dim]];
				}
			}
			if(n_n_level==type_int(0)) return std::numeric_limits<type_val>::quiet_NaN();
			Type_ArrayTemp<type_val> list_val(n_level); //`n_level' is greater than 0;
			list_val[0]=cv1;
			for(i_level=1;i_level<n_level;++i_level){
				list_val[i_level]=list_val[i_level-type_int(1)]*val;
			}
			ContrastSimil_quan.resize(n_n_level);
			type_int i0_level,i1_level,i_simil;
			typename type_tab_polyval::type_mat_real const * ptr_tab_polyval;
			i_dim=i_n_level=n_level=0;
			while(i_n_level<n_n_level){
				while(num_level_quan[id_fac_quan[i_dim]]==n_level) ++i_dim;
				n_level=num_level_quan[id_fac_quan[i_dim]];
				ContrastSimil_quan[i_n_level].resize(n_level*n_level);
				type_this::m_tab_polyval.Cal_polyval(n_level);
				ptr_tab_polyval=std::addressof(type_this::m_tab_polyval[n_level]);
				for(i0_level=0;i0_level<n_level;++i0_level){
					for(i1_level=0;i1_level<i0_level;++i1_level){
						ContrastSimil_quan[i_n_level][i0_level*n_level+i1_level]=
							ContrastSimil_quan[i_n_level][i1_level*n_level+i0_level];
					}
					for(;i1_level<n_level;++i1_level){
						i_simil=i0_level*n_level+i1_level;
						ContrastSimil_quan[i_n_level][i_simil]=cv1;
						for(i_level=1;i_level<n_level;++i_level){
							ContrastSimil_quan[i_n_level][i_simil]+=type_val(
								(*ptr_tab_polyval)(i_level,i0_level)*
								(*ptr_tab_polyval)(i_level,i1_level)*
								list_val[i_level]);
						}
					}
				}
				++i_n_level;
			}
		}
		//}
		type_val res_1;
		type_val res_ndiag=Liuze::math::zero(val),res_diag=res_ndiag;
		type_real wtsum=0;
		ForwardIterType_des iter0_des=des_begin,iter1_des;
		ForwardIterType_weight iter0_wt=weight_begin,iter1_wt;
		while(iter0_des!=des_end && iter0_wt!=weight_end){
			iter1_des=des_begin;
			iter1_wt=weight_begin;
			while(true){
				res_1=cv1;
				i_n_level=n_level=0;
				for(i_dim=0;i_dim<dim_qual;++i_dim){
					if(num_level_qual[id_fac_qual[i_dim]]>n_level){
						n_level=num_level_qual[id_fac_qual[i_dim]];
						++i_n_level;
					}
					res_1*=ContrastSimil_qual[i_n_level*type_int(2)-
						((*iter0_des)(id_fac_qual[i_dim],0)==(*iter1_des)(id_fac_qual[i_dim],0) ?
						type_int(1) : type_int(2))];
				}
				i_n_level=n_level=0;
				for(i_dim=0;i_dim<dim_quan;++i_dim){
					if(num_level_quan[id_fac_quan[i_dim]]>n_level){
						n_level=num_level_quan[id_fac_quan[i_dim]];
						++i_n_level;
					}
					res_1*=ContrastSimil_quan[i_n_level-type_int(1)]
						[(*iter0_des)(dim_qual+id_fac_quan[i_dim],0)*n_level
						+(*iter1_des)(dim_qual+id_fac_quan[i_dim],0)];
				}
				if(iter1_des!=iter0_des){
					res_ndiag+=(*iter0_wt)*(*iter1_wt)*res_1;
				} else {
					res_diag+=(*iter0_wt)*(*iter1_wt)*res_1;
					break;
				}
				++iter1_des;
				++iter1_wt;
			}
			wtsum+=(*iter0_wt);
			++iter0_des;
			++iter0_wt;
		}
		return type_val((type_real(1)/(wtsum*wtsum))*(res_ndiag+res_ndiag+res_diag));
	}
	
	//static function:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_polyval::type_mat_real
	WordlengthPattern<RealType,IntType,TraitsCheck>::get_contrast
	(type_int const & num_level){
		if(num_level<type_int(1)){
			return typename type_tab_polyval::type_mat_real(0,0);
		}
		type_this::m_tab_polyval.Cal_polyval(num_level);
		return type_this::m_tab_polyval[num_level];
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_haarval::type_mat_real
	WordlengthPattern<RealType,IntType,TraitsCheck>::get_contrast_Haar_value
	(type_int const & num_level,type_int const & base){
		if(num_level<type_int(1) || base<type_int(2)){
			return typename type_tab_haarval::type_mat_real(0,0);
		}
		type_this::m_tab_haarval.Cal_val(num_level,base);
		return std::get<0>(type_this::m_tab_haarval.get(num_level,base));
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_haarval::type_arr_int
	WordlengthPattern<RealType,IntType,TraitsCheck>::get_contrast_Haar_order
	(type_int const & num_level,type_int const & base){
		if(num_level<type_int(1) || base<type_int(2)){
			return typename type_tab_haarval::type_arr_int(0);
		}
		type_this::m_tab_haarval.Cal_val(num_level,base);
		return std::get<1>(type_this::m_tab_haarval.get(num_level,base));
	}
	
//protected:
	//types:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_polyval{
		public:
			//types:
			typedef Type_MatTemp<type_real> type_mat_real;
		protected:
			//types:
			typedef std::map<type_int,type_mat_real*> type_map_addr_valtab;
			class type_tab_powval;
		public:
			//Constructor:
			type_tab_polyval(void):
			m_tab_powval(),m_addr_tab_polyval(){
			}
			type_tab_polyval(type_tab_polyval const &)=delete;
			type_tab_polyval(type_tab_polyval &&)=delete;
			
			//Destructor:
			~type_tab_polyval(void){
				for(typename type_map_addr_valtab::iterator iter=m_addr_tab_polyval.begin();
				iter!=m_addr_tab_polyval.end();++iter){
					if(iter->second) delete iter->second;
				}
			}
			
			//operator:
			type_tab_polyval & operator =(type_tab_polyval const &)=delete;
			type_tab_polyval & operator =(type_tab_polyval &&)=delete;
			
			//other function:
			type_real operator ()
			(type_int const & num_level,type_int const & degree,type_int const & level){
				if(num_level<type_int(1)) return std::numeric_limits<type_real>::quiet_NaN();
				if(degree<type_int(0) || degree>=num_level || level<type_int(0) || level>=num_level){
					return std::numeric_limits<type_real>::quiet_NaN();
				}
				if(m_addr_tab_polyval.find(num_level)==m_addr_tab_polyval.end()){
					this->Cal_polyval(num_level);
				}
				return (*(m_addr_tab_polyval.at(num_level)))(degree,level);
			}
			type_mat_real const & operator [](type_int const & num_level) const {
				return *(m_addr_tab_polyval.at(num_level));
			}
			void Cal_polyval(type_int const & num_level){
				if(num_level<type_int(1)) return;
				if(m_addr_tab_polyval.find(num_level)!=m_addr_tab_polyval.end()) return;
				m_addr_tab_polyval.insert(std::pair<type_int const,type_mat_real*>
					(num_level,new type_mat_real(num_level,num_level)));
				type_mat_real & m_tab_polyval=*m_addr_tab_polyval[num_level];
				m_tab_powval.Set_enlarge(num_level,num_level);
				type_int i_poly,i_level,i_pow,j_poly;
				m_tab_polyval.row(0)=type_mat_real::Ones(1,num_level);
				type_mat_real poly_LE_coef,poly_LE_val,poly_LE_sol;
				type_real sqnorm;
				for(i_poly=1;i_poly<num_level;++i_poly){
					poly_LE_coef=type_mat_real::Zero(i_poly,i_poly);
					poly_LE_val=type_mat_real::Zero(i_poly,1);
					for(j_poly=0;j_poly<i_poly;++j_poly){
						for(i_level=0;i_level<num_level;++i_level){
							for(i_pow=0;i_pow<i_poly;++i_pow){
								poly_LE_coef(j_poly,i_pow)+=
									m_tab_powval[i_level][i_pow]*m_tab_polyval(j_poly,i_level);
							}
							poly_LE_val(j_poly,0)-=
								m_tab_powval[i_level][i_poly]*m_tab_polyval(j_poly,i_level);
						}
					}
					poly_LE_sol=poly_LE_coef.colPivHouseholderQr().solve(poly_LE_val);
					//`poly_LE_sol' stores the lower `i_poly' coefficients of\
					//the monic `i_poly'-th degree polynomial;
					m_tab_polyval.row(i_poly)=type_mat_real::Constant(1,num_level,poly_LE_sol(0,0));
					sqnorm=type_real(0);
					for(i_level=0;i_level<num_level;++i_level){
						for(i_pow=1;i_pow<i_poly;++i_pow){
							m_tab_polyval(i_poly,i_level)+=
								poly_LE_sol(i_pow,0)*m_tab_powval[i_level][i_pow];
						}
						m_tab_polyval(i_poly,i_level)+=m_tab_powval[i_level][i_poly];
						sqnorm+=m_tab_polyval(i_poly,i_level)*m_tab_polyval(i_poly,i_level);
					}
					sqnorm=type_real(sqrt(type_real(num_level)/sqnorm));
					//now `sqnorm' is the `i_poly'-th coefficient of the `i_poly'-th degree polynomial;
					m_tab_polyval.row(i_poly)*=sqnorm;
				}
			}
		protected:
			//types:
			class type_tab_powval{
				public:
					//types:
					typedef Type_ArrayTemp<Type_ArrayTemp<type_real> > type_arr2d_real;
					
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
									m_val[i_base][i_pow]=
										m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
								}
							}
						}
					}
					type_tab_powval(type_tab_powval const &)=default;
					type_tab_powval(type_tab_powval &&)=default;
					
					//Destructor:
					~type_tab_powval(void)=default;
					
					//operator:
					type_tab_powval & operator =(type_tab_powval const &)=default;
					type_tab_powval & operator =(type_tab_powval &&)=default;
					inline
					type_real operator ()(type_int const & base,type_int const & power) const {
						return (base>=type_int(0) && base<Get_num_base() &&
							power>=type_int(0) && power<Get_num_power()) ?
							m_val[base][power] : std::numeric_limits<type_real>::quiet_NaN();
					}
					inline
					Type_ArrayTemp<type_real> const & operator [](type_int const & base) const {
						return m_val[base];
					}
					
					//other function:
					inline type_int Get_num_base(void) const {
						return type_int(m_val.size());
					}
					inline type_int Get_num_power(void) const {
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
									m_val[i_base][i_pow]=
										m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
								}
							}
						}
						if(nn_base>n_base){
							m_val.resize(nn_base);
							for(i_base=n_base;i_base<nn_base;++i_base){
								m_val[i_base].resize(nn_pow);
								m_val[i_base][0]=type_real(1);
								for(i_pow=1;i_pow<nn_pow;++i_pow){
									m_val[i_base][i_pow]=
										m_val[i_base][i_pow-type_int(1)]*type_real(i_base);
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
					type_arr2d_real m_val;
			}; //class type_tab_powval;
			
			//member values:
			type_tab_powval m_tab_powval;
			type_map_addr_valtab m_addr_tab_polyval;
	}; //class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_polyval;
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_haarval{
		public:
			//types:
			typedef Type_ArrayTemp<type_int> type_arr_int;
			typedef Type_MatTemp<type_real> type_mat_real;
			typedef std::tuple<type_mat_real,type_arr_int> type_valtab;
			typedef std::array<type_int,(size_t)2> type_mapkey;
			class type_mapkey_less{
				public:
				typedef type_mapkey first_argument_type;
				typedef type_mapkey second_argument_type;
				typedef bool result_type;
				bool operator ()(type_mapkey const & k0,type_mapkey const & k1) const {
					for(size_t i=0;i<2;++i){
						if(k0[i]<k1[i]) return true;
						if(k0[i]>k1[i]) return false;
					}
					return false;
				}
			}; //class type_mapkey_less;
		protected:
			//types:
			typedef std::map<type_mapkey,type_valtab*,type_mapkey_less> type_map_addr_valtab;
			
		public:
			//Constructor:
			type_tab_haarval(type_tab_polyval & tab_polyval):
			m_addr_valtab(),ref_tab_polyval(tab_polyval){
			}
			type_tab_haarval(type_tab_haarval const &)=delete;
			type_tab_haarval(type_tab_haarval &&)=delete;
			
			//Destructor:
			~type_tab_haarval(void){
				for(typename type_map_addr_valtab::iterator iter=m_addr_valtab.begin();
				iter!=m_addr_valtab.end();++iter){
					if(iter->second) delete iter->second;
				}
			}
			
			//operator:
			type_tab_haarval & operator =(type_tab_haarval const &)=delete;
			type_tab_haarval & operator =(type_tab_haarval &&)=delete;
			
			//other function:
			static inline type_mapkey make_key(type_int const & num_level,type_int const & base){
				return type_mapkey({num_level,base});
			}
			inline type_valtab const & get(type_int const & num_level,type_int const & base) const {
				return *(m_addr_valtab.at(make_key(num_level,base)));
			}
			void Cal_val(type_int const & num_level,type_int const & base){
				if(num_level<type_int(1) || base<type_int(2)) return;
				type_mapkey key=make_key(num_level,base);
				if(m_addr_valtab.find(key)!=m_addr_valtab.end()) return;
				//checkings end.
				m_addr_valtab.insert
					(std::pair<type_mapkey const,type_valtab*>(key,
					new type_valtab(type_mat_real::Zero(num_level,num_level),type_arr_int(num_level))
					));
				type_mat_real & m_valtab_val=std::get<0>(*m_addr_valtab[key]);
				type_arr_int & m_valtab_ord=std::get<1>(*m_addr_valtab[key]);
				type_int i_level,i_order,i_fun,i_part,i_free;
				type_int i_fun_0;
				type_real norm;
				type_real norm0=type_real(sqrt(type_real(num_level)));
				typename type_tab_polyval::type_mat_real const * ptr_tab_polyval;
				type_real haarval;
				//Haar-like function of order 0{
				haarval=type_real(1)/norm0;
				for(i_level=0;i_level<num_level;++i_level) m_valtab_val(0,i_level)=haarval;
				m_valtab_ord[0]=type_int(0);
				//}
				//getting other Haar-like functions by a pre-order traversal of a tree:
				type_int order_max=type_int(ceil(log(type_real(num_level))/log(type_real(base))));
				Type_ArrayTemp<Type_ArrayTemp<type_int> > tree_part
					(order_max+type_int(1),Type_ArrayTemp<type_int>(base+type_int(1)));
				tree_part[0]={0,num_level};
				type_arr_int i_tree_part(order_max+type_int(1),type_int(0));
				type_int part_len,part_sublen,part_n_minor,part_base;
				i_order=i_fun=type_int(1);
				while(true){
					//the length of the discrete interval to be partitioned:
					part_len=
						tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]+type_int(1)]
						-tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]];
					//turning back if the interval cannot be partitioned:
					if(part_len==type_int(1)){
						--i_order;
						while(i_order>type_int(0)){
							++i_tree_part[i_order];
							if(i_tree_part[i_order]<type_int(tree_part[i_order].size())-type_int(1)){
								break;
							}
							--i_order;
						}
						if(i_order==type_int(0)) break;
						++i_order;
						continue;
					}
					part_sublen=part_len/base; //the length of the shorter sub-interval;
					part_n_minor=part_len%base; //the number of the longer sub-interval;
					//the actual base number for this interval:
					part_base= part_sublen==type_int(0) ? part_n_minor : base;
					part_n_minor=part_base-part_n_minor; //the number of the shorter sub-interval;
					//partitioning the interval (recording the begin and end indices){
					tree_part[i_order].resize(part_base+type_int(1));
					tree_part[i_order][0]=
						tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]];
					for(i_part=1;i_part<=part_n_minor;++i_part){
						tree_part[i_order][i_part]=
							tree_part[i_order][i_part-type_int(1)]+part_sublen;
					}
					++part_sublen; //the length of the longer sub-interval;
					for(;i_part<=part_base;++i_part){
						tree_part[i_order][i_part]=
							tree_part[i_order][i_part-type_int(1)]+part_sublen;
					}
					//}
					//calculating the Haar-like function values{
					//functions orthogonal to lower-order Haar-like functions{
					ref_tab_polyval.Cal_polyval(part_base);
					ptr_tab_polyval=std::addressof(ref_tab_polyval[part_base]);
					i_fun_0=i_fun;
					for(i_free=1;i_free<part_base;++i_free,++i_fun){
						for(i_part=0;i_part<part_base;++i_part){
							haarval=(*ptr_tab_polyval)(i_free,i_part);
							haarval/=type_real(tree_part[i_order][i_part+type_int(1)]
								-tree_part[i_order][i_part]);
							for(i_level=tree_part[i_order][i_part];
							i_level<tree_part[i_order][i_part+type_int(1)];++i_level){
								m_valtab_val(i_fun,i_level)=haarval;
							}
						}
						m_valtab_ord[i_fun]=i_order;
					}
					//}
					//orthonormalisation of the new functions{
					i_fun=i_fun_0;
					for(i_free=1;i_free<part_base;++i_free,++i_fun){
						for(i_level=i_fun_0;i_level<i_fun;++i_level){
							m_valtab_val.block(i_fun,tree_part[i_order][0],
							1,tree_part[i_order][part_base]-tree_part[i_order][0])
								-=
								(m_valtab_val.block(i_fun,tree_part[i_order][0],
								1,tree_part[i_order][part_base]-tree_part[i_order][0])
								*m_valtab_val.block(i_level,tree_part[i_order][0],
								1,tree_part[i_order][part_base]-tree_part[i_order][0]).transpose())
								(0,0)
								*m_valtab_val.block(i_level,tree_part[i_order][0],
								1,tree_part[i_order][part_base]-tree_part[i_order][0]);
						}
						norm=m_valtab_val.block(i_fun,tree_part[i_order][0],
							1,tree_part[i_order][part_base]-tree_part[i_order][0]).norm();
						m_valtab_val.block(i_fun,tree_part[i_order][0],
						1,tree_part[i_order][part_base]-tree_part[i_order][0])
							/=norm;
					}
					//}
					//}
					if(i_fun==num_level) break;
					//locating the node of the tree:
					i_tree_part[i_order]=type_int(0);
					++i_order;
				}
				m_valtab_val*=norm0;
			}
		protected:
			//member values:
			type_map_addr_valtab m_addr_valtab;
			type_tab_polyval & ref_tab_polyval;
	}; //class WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_haarval;
	
	//static member values:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_polyval
	WordlengthPattern<RealType,IntType,TraitsCheck>::m_tab_polyval;
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename WordlengthPattern<RealType,IntType,TraitsCheck>::type_tab_haarval
	WordlengthPattern<RealType,IntType,TraitsCheck>::m_tab_haarval
	(WordlengthPattern<RealType,IntType,TraitsCheck>::m_tab_polyval);
	
//}(class WordlengthPattern)
/****************************************************************************************************/
//class StratificationPattern{
//public:
	//other function:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,typename>
	typename StratificationPattern<RealType,IntType,TraitsCheck>::type_res
	StratificationPattern<RealType,IntType,TraitsCheck>::Stratification
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base) const {
		typedef Type_CompTemp<type_real> type_comp;
		typedef typename std::make_signed<type_int>::type type_sint;
		type_int dim=num_level.size();
		type_int i_dim,i_ord;
		type_int ord_max=0;
		//determining the maximum order:
		if(base.size()==0){
			for(i_dim=0;i_dim<dim;++i_dim){
				ord_max+=type_int(ceil(
					log(type_real(num_level[i_dim]))/log(type_real(2))));
			}
		} else if(type_int(base.size())<dim){
			i_ord=type_int(base.size());
			for(i_dim=0;i_dim<i_ord;++i_dim){
				ord_max+=type_int(ceil(
					log(type_real(num_level[i_dim]))/log(type_real(base[i_dim]))));
			}
			for(--i_ord;i_dim<dim;++i_dim){
				ord_max+=type_int(ceil(
					log(type_real(num_level[i_dim]))/log(type_real(base[i_ord]))));
			}
		} else {
			for(i_dim=0;i_dim<dim;++i_dim){
				ord_max+=type_int(ceil(
					log(type_real(num_level[i_dim]))/log(type_real(base[i_dim]))));
			}
		}
		if(ord_max==type_int(0)) return type_res(1,type_real(1));
		type_res res(ord_max+type_int(1)); //the result;
		res[0]=type_real(1);
		type_real grid=type_real(LZ_DEF_const_math_2pi)/type_real(ord_max);
		Type_ArrayTemp<type_comp> list_root(ord_max); //the unit roots;
		for(i_ord=0;i_ord<ord_max;++i_ord){
			list_root[i_ord]=type_comp
				(type_real(cos(type_real(i_ord)*grid)),type_real(sin(type_real(i_ord)*grid)));
		}
		Type_ArrayTemp<type_comp> enum_val(ord_max); //the values of enumerators;
		for(i_ord=0;i_ord<ord_max;++i_ord){
			enum_val[i_ord]=this->enumerator_Stratification
				(des_begin,des_end,weight_begin,weight_end,num_level,base,list_root[i_ord])
				-type_comp(1);
		}
		grid=type_real(ord_max);
		type_comp res0;
		for(i_ord=1;i_ord<=ord_max;++i_ord){
			res0=0;
			for(i_dim=0;i_dim<ord_max;++i_dim){
				res0+=enum_val[i_dim]*
					list_root[Liuze::math::mod(-type_sint(i_ord*i_dim),type_sint(ord_max))];
			}
			res[i_ord]=res0.real()/grid;
		}
		return res;
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,typename>
	typename StratificationPattern<RealType,IntType,TraitsCheck>::type_res
	StratificationPattern<RealType,IntType,TraitsCheck>::Stratification
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base,
	type_int order_max,type_real err) const {
		typedef Type_CompTemp<type_real> type_comp;
		typedef typename std::make_signed<type_int>::type type_sint;
		//the full stratification pattern:
		type_res res=
			this->Stratification(des_begin,des_end,weight_begin,weight_end,num_level,base);
		type_int order_full=type_int(res.size())-type_int(1); //the full maximum order;
		if(order_max<=type_int(0) || order_max>=order_full) return res;
		if(err<=type_real(0) || err>type_real(1)) err=type_real(LZ_DEF_const_default_precision);
		type_int i_ord;
		for(i_ord=0;i_ord<=order_full;++i_ord) res[i_ord]=Liuze::math::abs(res[i_ord]);
		//a bisearch to find the largest feasible decay{
		//the function determining whether a decay is enough (return a non-positive value) or not:
		std::function<type_real(type_real const &,type_int const &)> fun_errdiff=
			[&res,&order_full,&order_max,&err]
			(type_real const & decay,type_int order)->type_real{
				type_real errdiff= res[order]>err ? (-err*res[order]) : (-err);
				order+=order_max;
				type_real pow=decay;
				for(;order<=order_full;order+=order_max,pow*=decay){
					errdiff+=pow*res[order];
				}
				return errdiff;
			};
		type_real decay=type_real(1);
		type_real errdiff,decayL,decayR,decayM;
		type_real grid=std::numeric_limits<type_real>::epsilon();
		Type_Bool flag_found=1;
		for(i_ord=1;i_ord<=order_max;++i_ord){
			errdiff=fun_errdiff(type_real(1),i_ord);
			if(errdiff<=type_real(0)) continue;
			decayL=type_real(0);
			decayR=type_real(1);
			flag_found=0;
			while(true){
				decayM=decayR-decayL;
				if(decayR<=grid || (flag_found && decayM<=err)) break;
				decayM=decayL+decayM/type_real(2);
				errdiff=fun_errdiff(decayM,i_ord);
				if(errdiff>type_real(0)){
					decayR=decayM;
				} else {
					flag_found=1;
					decayL=decayM;
					if(errdiff==type_real(0)) break;
				}
			}
			if(flag_found){
				decay=std::min(decay,decayL);
			} else {
				break;
			}
		}
		if(!flag_found) return type_res(0);
		decay=std::pow(decay,type_real(1)/type_real(order_max));
		//}(a bisearch to find the largest feasible decay)
		res.resize(order_max+type_int(1)); //the result;
		res[0]=type_real(1);
		grid=type_real(LZ_DEF_const_math_2pi)/type_real(order_max);
		Type_ArrayTemp<type_comp> list_root(order_max); //the unit roots;
		for(i_ord=0;i_ord<order_max;++i_ord){
			list_root[i_ord]=type_comp
				(type_real(cos(type_real(i_ord)*grid)),type_real(sin(type_real(i_ord)*grid)));
		}
		Type_ArrayTemp<type_comp> enum_val(order_max); //the values of enumerators;
		for(i_ord=0;i_ord<order_max;++i_ord){
			enum_val[i_ord]=this->enumerator_Stratification
				(des_begin,des_end,weight_begin,weight_end,num_level,base,decay*list_root[i_ord])
				-type_comp(1);
		}
		grid=type_real(order_max);
		type_comp res0;
		type_int i_ord_1;
		for(i_ord=1;i_ord<=order_max;++i_ord){
			res0=0;
			for(i_ord_1=0;i_ord_1<order_max;++i_ord_1){
				res0+=enum_val[i_ord_1]*
					list_root[Liuze::math::mod(-type_sint(i_ord*i_ord_1),type_sint(order_max))];
			}
			grid*=decay;
			res[i_ord]=res0.real()/grid;
		}
		return res;
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	template<typename ForwardIterType_des,typename ForwardIterType_weight,
		typename ValType,typename>
	ValType
	StratificationPattern<RealType,IntType,TraitsCheck>::enumerator_Stratification
	(ForwardIterType_des const & des_begin,ForwardIterType_des const & des_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
	Type_ArrayTemp<type_int> const & num_level,Type_ArrayTemp<type_int> const & base,
	ValType const & val) const {
		typedef ValType type_val;
		typedef typename type_tab_haarval::type_mapkey type_key;
		type_int dim=num_level.size();
		type_val cv1=Liuze::math::one(val);
		if(dim==type_int(0)) return cv1;
		type_int i_dim,n_key;
		//making keys{
		Type_ArrayTemp<type_key> list_key(dim);
		if(base.size()==0){
			for(i_dim=0;i_dim<dim;++i_dim){
				list_key[i_dim]=type_tab_haarval::make_key(num_level[i_dim],type_int(2));
			}
		} else if(type_int(base.size())<dim){
			n_key=type_int(base.size());
			for(i_dim=0;i_dim<n_key;++i_dim){
				list_key[i_dim]=type_tab_haarval::make_key(num_level[i_dim],base[i_dim]);
			}
			for(--n_key;i_dim<dim;++i_dim){
				list_key[i_dim]=type_tab_haarval::make_key(num_level[i_dim],base[n_key]);
			}
		} else {
			for(i_dim=0;i_dim<dim;++i_dim){
				list_key[i_dim]=type_tab_haarval::make_key(num_level[i_dim],base[i_dim]);
			}
		}
		//}
		//re-ordering the factors according to their keys{
		Type_ArrayTemp<type_int> id_fac(dim);
		for(i_dim=0;i_dim<dim;++i_dim) id_fac[i_dim]=i_dim;
		std::sort(id_fac.begin(),id_fac.end(),
			[&list_key](type_int const & i0,type_int const & i1)->bool
			{return typename type_tab_haarval::type_mapkey_less()(list_key[i0],list_key[i1]);});
		type_key key=type_tab_haarval::make_key(type_int(0),type_int(1));
		n_key=type_int(0);
		for(i_dim=0;i_dim<dim;++i_dim){
			if(typename type_tab_haarval::type_mapkey_less()(key,list_key[id_fac[i_dim]])){
				++n_key;
				key=list_key[id_fac[i_dim]];
			}
		}
		//`n_key' is the number of different keys;
		if(n_key==type_int(0)) return std::numeric_limits<type_val>::quiet_NaN();
		//}
		//getting the powers of `val'{
		type_int i_fun;
		type_int deg_max=0;
		for(i_dim=0;i_dim<dim;++i_dim){
			i_fun=type_int(ceil(
				log(type_real(list_key[i_dim][0]))/log(type_real(list_key[i_dim][1]))));
			if(i_fun>deg_max) deg_max=i_fun;
		}
		Type_ArrayTemp<type_val> list_val(deg_max+type_int(1));
		list_val[0]=cv1;
		for(i_fun=1;i_fun<=deg_max;++i_fun){
			list_val[i_fun]=list_val[i_fun-type_int(1)]*val;
		}
		//}
		//contrast similarity{
		Type_ArrayTemp<Type_ArrayTemp<type_val> > ContrastSimil(n_key);
		type_int i0_level,i1_level,i_key,i_simil;
		typename type_tab_haarval::type_valtab const * ptr_valtab;
		key=type_tab_haarval::make_key(type_int(0),type_int(1));
		i_dim=i_key=0;
		while(i_key<n_key){
			while(list_key[id_fac[i_dim]]==key) ++i_dim;
			key=list_key[id_fac[i_dim]];
			ContrastSimil[i_key].resize(key[0]*key[0]);
			type_this::m_tab_haarval.Cal_val(key[0],key[1]);
			ptr_valtab=std::addressof(type_this::m_tab_haarval.get(key[0],key[1]));
			for(i0_level=0;i0_level<key[0];++i0_level){
				for(i1_level=0;i1_level<i0_level;++i1_level){
					ContrastSimil[i_key][i0_level*key[0]+i1_level]=
						ContrastSimil[i_key][i1_level*key[0]+i0_level];
				}
				for(;i1_level<key[0];++i1_level){
					i_simil=i0_level*key[0]+i1_level;
					ContrastSimil[i_key][i_simil]=cv1;
					for(i_fun=1;i_fun<key[0];++i_fun){
						ContrastSimil[i_key][i_simil]+=type_val(
							std::get<0>(*ptr_valtab)(i_fun,i0_level)
							*std::get<0>(*ptr_valtab)(i_fun,i1_level)
							*list_val[std::get<1>(*ptr_valtab)[i_fun]]);
					}
				}
			}
			++i_key;
		}
		//}
		type_val res_1;
		type_val res_ndiag=Liuze::math::zero(val),res_diag=res_ndiag;
		type_real wtsum=0;
		ForwardIterType_des iter0_des=des_begin,iter1_des;
		ForwardIterType_weight iter0_wt=weight_begin,iter1_wt;
		while(iter0_des!=des_end && iter0_wt!=weight_end){
			iter1_des=des_begin;
			iter1_wt=weight_begin;
			while(true){
				res_1=cv1;
				key=type_tab_haarval::make_key(type_int(0),type_int(1));
				i_key=0;
				for(i_dim=0;i_dim<dim;++i_dim){
					if(typename type_tab_haarval::type_mapkey_less()(key,list_key[id_fac[i_dim]])){
						key=list_key[id_fac[i_dim]];
						++i_key;
					}
					res_1*=ContrastSimil[i_key-type_int(1)]
						[(*iter0_des)(id_fac[i_dim],0)*key[0]+(*iter1_des)(id_fac[i_dim],0)];
				}
				if(iter1_des!=iter0_des){
					res_ndiag+=(*iter0_wt)*(*iter1_wt)*res_1;
				} else {
					res_diag+=(*iter0_wt)*(*iter1_wt)*res_1;
					break;
				}
				++iter1_des;
				++iter1_wt;
			}
			wtsum+=(*iter0_wt);
			++iter0_des;
			++iter0_wt;
		}
		return type_val((type_real(1)/(wtsum*wtsum))*(res_ndiag+res_ndiag+res_diag));
	}
	
	//static function:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename StratificationPattern<RealType,IntType,TraitsCheck>::type_tab_haarval::type_mat_real
	StratificationPattern<RealType,IntType,TraitsCheck>::get_contrast_Haar_value
	(type_int const & num_level,type_int const & base){
		if(num_level<type_int(1) || base<type_int(2)){
			return typename type_tab_haarval::type_mat_real(0,0);
		}
		type_this::m_tab_haarval.Cal_val(num_level,base);
		return std::get<0>(type_this::m_tab_haarval.get(num_level,base));
	}
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename StratificationPattern<RealType,IntType,TraitsCheck>::type_tab_haarval::type_arr_int
	StratificationPattern<RealType,IntType,TraitsCheck>::get_contrast_Haar_order
	(type_int const & num_level,type_int const & base){
		if(num_level<type_int(1) || base<type_int(2)){
			return typename type_tab_haarval::type_arr_int(0);
		}
		type_this::m_tab_haarval.Cal_val(num_level,base);
		return std::get<1>(type_this::m_tab_haarval.get(num_level,base));
	}
	
//protected:
	//types:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	class StratificationPattern<RealType,IntType,TraitsCheck>::type_tab_haarval{
		public:
			//types:
			typedef Type_ArrayTemp<type_int> type_arr_int;
			typedef Type_MatTemp<type_real> type_mat_real;
			typedef std::tuple<type_mat_real,type_arr_int> type_valtab;
			typedef std::array<type_int,(size_t)2> type_mapkey;
			class type_mapkey_less{
				public:
				typedef type_mapkey first_argument_type;
				typedef type_mapkey second_argument_type;
				typedef bool result_type;
				bool operator ()(type_mapkey const & k0,type_mapkey const & k1) const {
					for(size_t i=0;i<2;++i){
						if(k0[i]<k1[i]) return true;
						if(k0[i]>k1[i]) return false;
					}
					return false;
				}
			}; //class type_mapkey_less;
		protected:
			//types:
			typedef std::map<type_mapkey,type_valtab*,type_mapkey_less> type_map_addr_valtab;
			
		public:
			//Constructor:
			type_tab_haarval(void):m_addr_valtab(){}
			type_tab_haarval(type_tab_haarval const &)=delete;
			type_tab_haarval(type_tab_haarval &&)=delete;
			
			//Destructor:
			~type_tab_haarval(void){
				for(typename type_map_addr_valtab::iterator iter=m_addr_valtab.begin();
				iter!=m_addr_valtab.end();++iter){
					if(iter->second) delete iter->second;
				}
			}
			
			//operator:
			type_tab_haarval & operator =(type_tab_haarval const &)=delete;
			type_tab_haarval & operator =(type_tab_haarval &&)=delete;
			
			//other function:
			static inline type_mapkey make_key(type_int const & num_level,type_int const & base){
				return type_mapkey({num_level,base});
			}
			inline type_valtab const & get(type_int const & num_level,type_int const & base) const {
				return *(m_addr_valtab.at(make_key(num_level,base)));
			}
			void Cal_val(type_int const & num_level,type_int const & base){
				if(num_level<type_int(1) || base<type_int(2)) return;
				type_mapkey key=make_key(num_level,base);
				if(m_addr_valtab.find(key)!=m_addr_valtab.end()) return;
				//checkings end.
				m_addr_valtab.insert
					(std::pair<type_mapkey const,type_valtab*>(key,
					new type_valtab(type_mat_real::Zero(num_level,num_level),type_arr_int(num_level))
					));
				type_mat_real & m_valtab_val=std::get<0>(*m_addr_valtab[key]);
				type_arr_int & m_valtab_ord=std::get<1>(*m_addr_valtab[key]);
				type_int i_level,i_order,i_fun,i_part,i_free;
				type_int i_fun_0;
				type_real norm;
				type_real norm0=type_real(sqrt(type_real(num_level)));
				type_real haarval;
				//Haar-like function of order 0{
				haarval=type_real(1)/norm0;
				for(i_level=0;i_level<num_level;++i_level) m_valtab_val(0,i_level)=haarval;
				m_valtab_ord[0]=type_int(0);
				//}
				//getting other Haar-like functions by a pre-order traversal of a tree:
				type_int order_max=type_int(ceil(log(type_real(num_level))/log(type_real(base))));
				Type_ArrayTemp<Type_ArrayTemp<type_int> > tree_part
					(order_max+type_int(1),Type_ArrayTemp<type_int>(base+type_int(1)));
				tree_part[0]={0,num_level};
				type_arr_int i_tree_part(order_max+type_int(1),type_int(0));
				type_int part_len,part_sublen,part_n_minor,part_base;
				i_order=i_fun=type_int(1);
				while(true){
					//the length of the discrete interval to be partitioned:
					part_len=
						tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]+type_int(1)]
						-tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]];
					//turning back if the interval cannot be partitioned:
					if(part_len==type_int(1)){
						--i_order;
						while(i_order>type_int(0)){
							++i_tree_part[i_order];
							if(i_tree_part[i_order]<type_int(tree_part[i_order].size())-type_int(1)){
								break;
							}
							--i_order;
						}
						if(i_order==type_int(0)) break;
						++i_order;
						continue;
					}
					part_sublen=part_len/base; //the length of the shorter sub-interval;
					part_n_minor=part_len%base; //the number of the longer sub-interval;
					//the actual base number for this interval:
					part_base= part_sublen==type_int(0) ? part_n_minor : base;
					part_n_minor=part_base-part_n_minor; //the number of the shorter sub-interval;
					//partitioning the interval (recording the begin and end indices){
					tree_part[i_order].resize(part_base+type_int(1));
					tree_part[i_order][0]=
						tree_part[i_order-type_int(1)][i_tree_part[i_order-type_int(1)]];
					for(i_part=1;i_part<=part_n_minor;++i_part){
						tree_part[i_order][i_part]=
							tree_part[i_order][i_part-type_int(1)]+part_sublen;
					}
					++part_sublen; //the length of the longer sub-interval;
					for(;i_part<=part_base;++i_part){
						tree_part[i_order][i_part]=
							tree_part[i_order][i_part-type_int(1)]+part_sublen;
					}
					//}
					//calculating the Haar-like function values{
					i_fun_0=i_fun;
					for(i_free=1;i_free<part_base;++i_free,++i_fun){
						haarval=type_real(i_free)
							/type_real(tree_part[i_order][i_free]-tree_part[i_order][0]);
						for(i_level=tree_part[i_order][0];
						i_level<tree_part[i_order][i_free];++i_level){
							m_valtab_val(i_fun,i_level)=haarval;
						}
						haarval=-type_real(i_free)
							/type_real(tree_part[i_order][i_free+type_int(1)]
							-tree_part[i_order][i_free]);
						for(i_level=tree_part[i_order][i_free];
						i_level<tree_part[i_order][i_free+type_int(1)];++i_level){
							m_valtab_val(i_fun,i_level)=haarval;
						}
						norm=m_valtab_val.block(i_fun,tree_part[i_order][0],
							1,tree_part[i_order][i_free+type_int(1)]-tree_part[i_order][0]).norm();
						m_valtab_val.block(i_fun,tree_part[i_order][0],
						1,tree_part[i_order][i_free+type_int(1)]-tree_part[i_order][0])
							/=norm;
						m_valtab_ord[i_fun]=i_order;
					}
					//}
					if(i_fun==num_level) break;
					//locating the node of the tree:
					i_tree_part[i_order]=type_int(0);
					++i_order;
				}
				m_valtab_val*=norm0;
			}
		protected:
			//member values:
			type_map_addr_valtab m_addr_valtab;
	}; //class StratificationPattern<RealType,IntType,TraitsCheck>::type_tab_haarval;
	
	//static member values:
	
	template<typename RealType,typename IntType,typename TraitsCheck>
	typename StratificationPattern<RealType,IntType,TraitsCheck>::type_tab_haarval
	StratificationPattern<RealType,IntType,TraitsCheck>::m_tab_haarval;
	
//}(class StratificationPattern)
/****************************************************************************************************/

} //namespace crit;
} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_aberration_CPP
#endif //#if LZ_DEF_LZ_H_stats_DOE_crit_aberration!=202105L