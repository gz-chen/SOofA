#if LZ_DEF_LZ_H_stats_DOE_crit_uniformity!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_uniformity_CPP
#define LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_uniformity_CPP 202105L

namespace Liuze{
namespace stats{
namespace crit{

/****************************************************************************************************/
//class Discrepancy_Leb_infty::CDF_Uniform_std{
//public:
	//Constructor:
	template<typename RealType>
	Discrepancy_Leb_infty<RealType>::CDF_Uniform_std::CDF_Uniform_std
	(type_size dimension):
	dim(dimension>(type_size)0 ? dimension : (type_size)1){
	}
	//CDF:
	template<typename RealType>
	typename Discrepancy_Leb_infty<RealType>::type_real
	Discrepancy_Leb_infty<RealType>::CDF_Uniform_std::operator ()
	(type_matrix const & x) const {
		if(x.rows()!=dim || x.cols()!=(type_size)1) return (type_real)0;
		type_real res=1;
		for(type_size i_dim=0;i_dim<dim;++i_dim){
			if(x(i_dim,0)<=(type_real)0) return (type_real)0;
			if(x(i_dim,0)<(type_real)1) res*=x(i_dim,0);
		}
		return res;
	}
//}(class Discrepancy_Leb_infty::CDF_Uniform_std)
//class Discrepancy_Leb_infty{
//public:
	//Constructor:
	
	template<typename RealType>
	Discrepancy_Leb_infty<RealType>::Discrepancy_Leb_infty
	(type_size dimension,type_real precision,type_stat cdf_def):
	pre(precision>(type_real)0 ? precision : LZ_DEF_const_default_precision),
	CDF_def(cdf_def){
		switch(CDF_def){
			case type_this::val_CDF_def_less: break;
			case type_this::val_CDF_def_lesseq: break;
			default: CDF_def=type_this::val_CDF_def_less;
		}
		dimension= dimension>(type_size)0 ? dimension : (type_size)1;
		fun_cdf=type_this::CDF_Uniform_std(dimension);
	}
	
	template<typename RealType>
	template<typename FuncDistrType>
	Discrepancy_Leb_infty<RealType>::Discrepancy_Leb_infty
	(FuncDistrType & func_distr,type_real precision,type_stat cdf_def):
	pre(precision>(type_real)0 ? precision : LZ_DEF_const_default_precision),
	CDF_def(cdf_def){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<FuncDistrType>::type(type_matrix)>::type,
			type_real>::value));
		switch(CDF_def){
			case type_this::val_CDF_def_less: break;
			case type_this::val_CDF_def_lesseq: break;
			default: CDF_def=type_this::val_CDF_def_less;
		}
		fun_cdf=func_distr;
	}
	
	//Calculate:
	
	template<typename RealType>
	template<typename ForwardIterType_sam,typename ForwardIterType_weight>
	typename Discrepancy_Leb_infty<RealType>::type_real
	Discrepancy_Leb_infty<RealType>::operator ()
	(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
		//following is the type checking.
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::iterator_traits<ForwardIterType_sam>::value_type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_real,
			typename std::iterator_traits<ForwardIterType_weight>::value_type>::value));
		//following is the value checking.
		if(sam_begin==sam_end || sam_begin->rows()==0 || weight_begin==weight_end){
			return type_real(0); //invalid parameters;
		}
		type_size dim=sam_begin->rows(),size_sam=0;
		type_size i_dim,i_knot; //id vars;
		ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
		ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
		while(iter_sam!=sam_end && iter_wt!=weight_end){
			if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
				return type_real(0); //invalid parameters;
			}
			++size_sam;
			++iter_sam;
			++iter_wt;
		}
		type_real precision=pre; //the precision;
		Type_ArrayTemp<Type_ArrayTemp<type_real> > knots(dim); //the lists of knots on all dimensions;
		Type_ArrayTemp<type_real> ta0; //temp_array_0;
		type_real tr0; //temp var;
		//the following constructs the knots on all dimensions.
		ta0.resize(size_sam);
		for(i_dim=0;i_dim<dim;++i_dim){
			i_knot=0;
			for(iter_sam=sam_begin;iter_sam!=sam_end;++iter_sam){
				ta0[i_knot]=(*iter_sam)(i_dim,0);
				++i_knot;
			}
			std::sort(ta0.begin(),ta0.end());
			knots[i_dim].resize(1,ta0[0]);
			for(i_knot=1;i_knot<size_sam;++i_knot){
				//filtering the duplicate coordinates:
				if(ta0[i_knot]>ta0[i_knot-1]){
					knots[i_dim].push_back(ta0[i_knot]);
					tr0=(ta0[i_knot]-ta0[i_knot-1])/(type_real)2;
					if(tr0<precision) precision=tr0; //re-tuning the precision;
				}
			}
			//the following determines the approximate positive infinity of this dimension.
			tr0=(type_real)5*std::max((type_real)1,precision);
			if(knots[i_dim].size()>(type_size)1){
				tr0*=(knots[i_dim][knots[i_dim].size()-1]-knots[i_dim][0]);
			} else {
				tr0*=std::max(Liuze::math::abs(knots[i_dim][0]),(type_real)1);
			}
			tr0+=(knots[i_dim][knots[i_dim].size()-1]);
			knots[i_dim].push_back(tr0);
		}
		//this part deals with the upper and lower values of the cdf.
		type_matrix disturb=type_matrix::Ones(dim,1)*(CDF_def==(Type_Stat)0 ? precision : -precision);
		type_real ecdf0,ecdf1,cdf0,cdf1; //store the values of ecdf and cdf under 2 definitions;
		//the following calculates the result.
		type_real res=0; //the return result;
		type_matrix knot_pos(dim,1); //one knot position;
		//indicating the knot is a key position (sam: at a sample; crd: the coordinates):
		std::vector<bool> is_key_pos_crd(dim),is_key_pos_crd_sam(dim),is_positive_infty(dim);
		Type_Bool is_key_pos_sam;
		//Following: the loop for each knot point.
		Type_ArrayTemp<type_size> loop_knot(dim,(type_size)0); //the loop var;
		while(true){
			//renew the knot point:
			for(i_dim=0;i_dim<dim;++i_dim) knot_pos(i_dim,0)=knots[i_dim][loop_knot[i_dim]];
			ecdf0=ecdf1=(type_real)0; //initialisation;
			is_key_pos_crd.assign(dim,false); //initialisation;
			//Following: the loop for each sample.
			iter_sam=sam_begin;
			iter_wt=weight_begin;
			while(iter_sam!=sam_end && iter_wt!=weight_end){
				is_key_pos_crd_sam.assign(dim,false); //initialisation;
				is_key_pos_sam=false; //initialisation;
				//Following: the loop for each coordinate.
				for(i_dim=0;i_dim<dim;++i_dim){
					//this sample is not dominated by the knot point:
					if((*iter_sam)(i_dim,0)>knot_pos(i_dim,0)) goto GTS_loop_end_crd;
					if((*iter_sam)(i_dim,0)==knot_pos(i_dim,0)){
						is_key_pos_sam=true;
						is_key_pos_crd_sam[i_dim]=true;
					}
				}
				//come here only when (*iter_sam) is component-wise "leq" knot_pos.
				ecdf1+=(*iter_wt);
				//if (*iter_sam) is component-wise less than knot_pos:
				if(is_key_pos_sam==false){
					ecdf0+=(*iter_wt);
				} else {
					//mark the dimension at which there are samples dominated by knot_pos:
					std::transform(is_key_pos_crd.begin(),is_key_pos_crd.end(),is_key_pos_crd_sam.begin(),
						is_key_pos_crd.begin(),std::logical_or<bool>());
				}
				GTS_loop_end_crd:
				++iter_sam;
				++iter_wt;
			}
			is_positive_infty.assign(dim,false);
			for(i_dim=0;i_dim<dim;++i_dim){
				//if knot_pos reaches the approximate positive infinity of this dimension:
				if(loop_knot[i_dim]+1==knots[i_dim].size()){
					is_positive_infty[i_dim]=is_key_pos_crd[i_dim]=true;
				}
			}
			//if knot_pos is a key knot point:
			if(std::accumulate(is_key_pos_crd.begin(),is_key_pos_crd.end(),true,std::logical_and<bool>()) &&
				!(std::accumulate(is_positive_infty.begin(),is_positive_infty.end(),true,std::logical_and<bool>()))){
				//Calculate the corresponding cdf:
				switch(CDF_def){
					case (Type_Stat)1:
						cdf0=fun_cdf(knot_pos+disturb);
						cdf1=fun_cdf(knot_pos);
						break;
					case (Type_Stat)0:
					default:
						cdf0=fun_cdf(knot_pos);
						cdf1=fun_cdf(knot_pos+disturb);
				}
				tr0=std::max(math::abs(cdf0-ecdf0),math::abs(cdf1-ecdf1));
				if(tr0>res) res=tr0;
			}
			//the following renews loop_knot;
			i_dim=dim-1;
			while(loop_knot[i_dim]+1==knots[i_dim].size()){
				if(i_dim==0) goto GTS_loop_end_knot;
				--i_dim;
			}
			++(loop_knot[i_dim]);
			while((++i_dim)<dim) loop_knot[i_dim]=(type_size)0;
		}
		GTS_loop_end_knot:
		return res;
	}
	
	//Get:
	
	template<typename RealType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_stat
	Discrepancy_Leb_infty<RealType>::Get_CDF_def
	() const {
		return CDF_def;
	}
	
	template<typename RealType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_cdf
	Discrepancy_Leb_infty<RealType>::Get_CDF
	() const {
		return fun_cdf;
	}
	
	template<typename RealType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_real
	Discrepancy_Leb_infty<RealType>::Get_precision
	() const {
		return pre;
	}
	
	//Set:
	
	template<typename RealType>
	typename Discrepancy_Leb_infty<RealType>::type_stat
	Discrepancy_Leb_infty<RealType>::Set_CDF_def
	(type_stat cdf_def){
		switch(cdf_def){
			case type_this::val_CDF_def_less:
			case type_this::val_CDF_def_lesseq:
				CDF_def=cdf_def;
		}
		return CDF_def;
	}
	
	template<typename RealType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_cdf
	Discrepancy_Leb_infty<RealType>::Set_CDF
	(type_size dimension){
		if(dimension>(type_size)0) fun_cdf=type_this::CDF_Uniform_std(dimension);
		return fun_cdf;
	}
	
	template<typename RealType>
	template<typename FuncDistrType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_cdf
	Discrepancy_Leb_infty<RealType>::Set_CDF
	(FuncDistrType & func_distr){
		LZ_DEF_func_check_traits((std::is_same<
			typename std::result_of<typename std::decay<FuncDistrType>::type(type_matrix)>::type,
			type_real>::value));
		return fun_cdf=func_distr;
	}
	
	template<typename RealType>
	inline
	typename Discrepancy_Leb_infty<RealType>::type_real
	Discrepancy_Leb_infty<RealType>::Set_precision
	(type_real precision){
		if(precision>(type_real)0) pre=precision;
		return pre;
	}
	
//}(class Discrepancy_Leb_infty)
/****************************************************************************************************/
//class Discrepancy_canon_Leb_infty_approxLambda{
//public:
	//Calculate:
	
	template<typename RealType>
	template<typename ForwardIterType_sam,typename ForwardIterType_weight>
	typename Discrepancy_canon_Leb_infty_approxLambda<RealType>::type_real
	Discrepancy_canon_Leb_infty_approxLambda<RealType>::operator ()
	(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
		//following is the type checking.
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::iterator_traits<ForwardIterType_sam>::value_type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_real,
			typename std::iterator_traits<ForwardIterType_weight>::value_type>::value));
		//following is the value checking.
		if(sam_begin==sam_end || sam_begin->rows()==0 || weight_begin==weight_end){
			return type_real(0); //invalid parameters;
		}
		type_size dim=sam_begin->rows();
		type_size i_dim; //id vars;
		ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
		ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
		while(iter_sam!=sam_end && iter_wt!=weight_end){
			if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
				return type_real(0); //invalid parameters;
			}
			++iter_sam;
			++iter_wt;
		}
		//the following calculates the result.
		type_real constexpr cpi=LZ_DEF_const_math_pi;
		type_real constexpr c2=2;
		type_real constexpr c1=1;
		type_real constexpr c2_pi=c2/cpi;
		type_real res=0,tr0;
		iter_wt=weight_begin;
		for(iter_sam=sam_begin;iter_sam!=sam_end;++iter_sam){
			tr0=c1;
			for(i_dim=(type_size)0;i_dim<dim;++i_dim){
				tr0*=(c1-c2_pi*log(c2*sin(cpi*((*iter_sam)(i_dim,0)))));
			}
			res+=(tr0*(*iter_wt));
			++iter_wt;
		}
		return res;
	}
	
//}(class Discrepancy_canon_Leb_infty_approxLambda)
/****************************************************************************************************/
//class Discrepancy_canon_Leb_2_kernel::type_kernel{
//public:
	//Constructor:
	template<typename RealType>
	Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel::type_kernel
	(type_size const & dimension):
	dim(dimension){
	}
	//Get:
	template<typename RealType>
	inline
	typename Discrepancy_canon_Leb_2_kernel<RealType>::type_size
	Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel::Get_dim
	() const {
		return dim;
	}
//}(class Discrepancy_canon_Leb_2_kernel::type_kernel)
//class Discrepancy_canon_Leb_2_kernel::type_kernel_centered{
template<typename RealType>
class Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel_centered:
public Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel{
	public:
		typedef type_kernel_centered type_this;
		typedef type_kernel type_base;
		
		//Constructor:
		type_kernel_centered(type_size dimension=1,type_real weight_projective=1):
		type_base(dimension>(type_size)0 ? dimension : (type_size)1),
		weight_project(weight_projective){
			val_exp=math::powint((type_real)1/(type_real)12+weight_project,this->dim);
		}
		//operator ():
		virtual type_real operator ()(type_matrix const & x0,type_matrix const & x1) const {
			//without check:
			//\if(x0.rows()<this->dim || x1.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1,cr1f2=0.5;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				res*=(weight_project+cr1f2*
					(math::abs(x0(i_dim,0)-cr1f2)+math::abs(x1(i_dim,0)-cr1f2)
					-math::abs(x0(i_dim,0)-x1(i_dim,0))));
			}
			return res;
		}
		virtual type_real operator ()(type_matrix const & x) const {
			//without check:
			//\if(x.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1,tr0,cr1f2=0.5;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				tr0=math::abs(x(i_dim,0)-cr1f2);
				res*=(weight_project+cr1f2*(tr0-tr0*tr0));
			}
			return res;
		}
		virtual type_real operator ()() const {
			return val_exp;
		}
	protected:
		type_real weight_project;
		type_real val_exp; //the expectation;
};
//}(class Discrepancy_canon_Leb_2_kernel::type_kernel_centered)
//class Discrepancy_canon_Leb_2_kernel::type_kernel_wrap{
template<typename RealType>
class Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel_wrap:
public Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel{
	public:
		typedef type_kernel_wrap type_this;
		typedef type_kernel type_base;
		
		//Constructor:
		type_kernel_wrap(type_size dimension=1,type_real weight_projective=1):
		type_base(dimension>(type_size)0 ? dimension : (type_size)1),
		weight_project(weight_projective){
			val_exp=math::powint((type_real)1/(type_real)3+weight_project,this->dim);
			coef2=(type_real)1/(type_real)2+weight_project;
		}
		//operator ():
		virtual type_real operator ()(type_matrix const & x0,type_matrix const & x1) const {
			//without check:
			//\if(x0.rows()<this->dim || x1.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1,tr0;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				tr0=math::abs(x0(i_dim,0)-x1(i_dim,0));
				res*=(coef2-tr0+tr0*tr0);
			}
			return res;
		}
		virtual type_real operator ()(type_matrix const & x) const {
			return val_exp;
		}
		virtual type_real operator ()() const {
			return val_exp;
		}
	protected:
		type_real weight_project;
		type_real coef2;
		type_real val_exp; //the expectation;
};
//}(class Discrepancy_canon_Leb_2_kernel::type_kernel_wrap)
//class Discrepancy_canon_Leb_2_kernel::type_kernel_mixture{
template<typename RealType>
class Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel_mixture:
public Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel{
	public:
		typedef type_kernel_mixture type_this;
		typedef type_kernel type_base;
		
		//Constructor:
		type_kernel_mixture(type_size dimension=1,type_real weight_projective=1):
		type_base(dimension>(type_size)0 ? dimension : (type_size)1),
		weight_project(weight_projective){
			val_exp=math::powint((type_real)7/(type_real)12+weight_project,this->dim);
			coef1=(type_real)2/(type_real)3+weight_project;
			coef2=(type_real)7/(type_real)8+weight_project;
		}
		//operator ():
		virtual type_real operator ()(type_matrix const & x0,type_matrix const & x1) const {
			//without check:
			//\if(x0.rows()<this->dim || x1.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1,tr0,cr1f2=0.5,cr1f4=0.25,cr3f4=0.75;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				tr0=math::abs(x0(i_dim,0)-x1(i_dim,0));
				res*=(coef2-cr1f4*(math::abs(x0(i_dim,0)-cr1f2)+math::abs(x1(i_dim,0)-cr1f2))
					-cr3f4*tr0+cr1f2*tr0*tr0);
			}
			return res;
		}
		virtual type_real operator ()(type_matrix const & x) const {
			//without check:
			//\if(x.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1,tr0,cr1f2=0.5,cr1f4=0.25;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				tr0=math::abs(x(i_dim,0)-cr1f2);
				res*=(coef1-cr1f4*(tr0+tr0*tr0));
			}
			return res;
		}
		virtual type_real operator ()() const {
			return val_exp;
		}
	protected:
		type_real weight_project;
		type_real coef1,coef2;
		type_real val_exp; //the expectation;
};
//}(class Discrepancy_canon_Leb_2_kernel::type_kernel_mixture)
//class Discrepancy_canon_Leb_2_kernel::type_kernel_discrete{
template<typename RealType>
class Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel_discrete:
public Discrepancy_canon_Leb_2_kernel<RealType>::type_kernel{
	public:
		typedef type_kernel_discrete type_this;
		typedef type_kernel type_base;
		typedef Type_ArrayTemp<type_real> type_real_arr;
		
		//Constructor:
		type_kernel_discrete(type_size dimension=1):
		type_base(dimension>(type_size)0 ? dimension : (type_size)1){
			type_real rdim=this->dim;
			n_level.assign(this->dim,this->dim);
			par_a.assign(this->dim,rdim-type_real(1)+type_real(1)/rdim);
			par_b.assign(this->dim,type_real(1)/rdim);
			par_err.assign(this->dim,type_real(0.5)/rdim);
			val_exp=type_real(1);
		}
		//Set:
		void Set_level(type_size const & num_level){
			if(num_level<=type_size(0)) return;
			n_level.assign(this->dim,num_level);
		}
		template<typename SizeAterType,
			typename=typename std::enable_if<std::is_convertible<
			typename SizeAterType::value_type,type_size>::value>::type>
		void Set_level(SizeAterType const & num_level){
			type_size i_dim;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				if(type_size(num_level[i_dim])<=type_size(0)) continue;
				n_level[i_dim]=type_size(num_level[i_dim]);
			}
		}
		void Set_par_a(type_real const & a){
			if(a<=type_real(0)) return;
			par_a.assign(this->dim,a);
		}
		template<typename RealAterType,
			typename=typename std::enable_if<std::is_convertible<
			typename RealAterType::value_type,type_real>::value>::type>
		void Set_par_a(RealAterType const & a){
			type_size i_dim;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				if(type_real(a[i_dim])<=type_real(0)) continue;
				par_a[i_dim]=type_real(a[i_dim]);
			}
		}
		void Set_par_b(type_real const & b){
			par_b.assign(this->dim,b);
		}
		template<typename RealAterType,
			typename=typename std::enable_if<std::is_convertible<
			typename RealAterType::value_type,type_real>::value>::type>
		void Set_par_b(RealAterType const & b){
			type_size i_dim;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				par_b[i_dim]=type_real(b[i_dim]);
			}
		}
		void Set_par_ab(void){
			type_size i_dim;
			type_real rlv;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				rlv=type_real(n_level[i_dim]);
				par_b[i_dim]=type_real(1)/rlv;
				par_a[i_dim]=rlv-type_real(1)+par_b[i_dim];
			}
		}
		void Set_init(void){
			type_size i_dim;
			val_exp=type_real(1);
			for(i_dim=0;i_dim<this->dim;++i_dim){
				val_exp*=((par_a[i_dim]+par_b[i_dim]*type_real(n_level[i_dim]-type_size(1)))
					/type_real(n_level[i_dim]));
				par_err[i_dim]=type_real(0.5)/type_real(n_level[i_dim]);
			}
		}
		//operator ():
		virtual type_real operator ()(type_matrix const & x0,type_matrix const & x1) const {
			//without check:
			//\if(x0.rows()<this->dim || x1.rows()<this->dim) return (type_real)0;
			type_size i_dim;
			type_real res=1;
			for(i_dim=0;i_dim<this->dim;++i_dim){
				res*=(math::abs(x0(i_dim,0)-x1(i_dim,0))>=par_err[i_dim] ?
					par_b[i_dim] : par_a[i_dim]);
			}
			return res;
		}
		virtual type_real operator ()(type_matrix const & x) const {
			return val_exp;
		}
		virtual type_real operator ()() const {
			return val_exp;
		}
	protected:
		Type_ArrayTemp<type_size> n_level;
		type_real_arr par_a,par_b,par_err;
		type_real val_exp; //the expectation;
};
//}(class Discrepancy_canon_Leb_2_kernel::type_kernel_discrete)
//class Discrepancy_canon_Leb_2_kernel{
//public:
	//Constructor:
	
	template<typename RealType>
	Discrepancy_canon_Leb_2_kernel<RealType>::Discrepancy_canon_Leb_2_kernel
	(type_kernel & kernel):
	pt_ker(std::addressof(kernel)),stat_ker(type_this::val_stat_ker_nonew){
	}
	
	template<typename RealType>
	Discrepancy_canon_Leb_2_kernel<RealType>::Discrepancy_canon_Leb_2_kernel
	(type_size dimension,type_stat kernel,type_real weight_projective){
		dimension= dimension>(type_size)0 ? dimension : (type_size)0;
		switch(kernel){
			case type_this::val_stat_ker_centered:
				pt_ker=new (std::nothrow)type_kernel_centered(dimension,weight_projective);
				break;
			case type_this::val_stat_ker_wrap:
				pt_ker=new (std::nothrow)type_kernel_wrap(dimension,weight_projective);
				break;
			case type_this::val_stat_ker_mixture:
				pt_ker=new (std::nothrow)type_kernel_mixture(dimension,weight_projective);
				break;
			default:
				pt_ker=NULL;
		}
		if(pt_ker==NULL) stat_ker=type_this::val_stat_ker_NULL;
		else stat_ker=kernel;
	}
	
	//Destructor:
	
	template<typename RealType>
	Discrepancy_canon_Leb_2_kernel<RealType>::~Discrepancy_canon_Leb_2_kernel
	(){
		if(stat_ker!=val_stat_ker_nonew && stat_ker!=val_stat_ker_NULL) delete pt_ker;
	}
	
	//Calculate:
	
	template<typename RealType>
	template<typename ForwardIterType_sam,typename ForwardIterType_weight>
	typename Discrepancy_canon_Leb_2_kernel<RealType>::type_real
	Discrepancy_canon_Leb_2_kernel<RealType>::operator ()
	(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
	ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
		//following is the type checking.
		LZ_DEF_func_check_traits((std::is_same<type_matrix,
			typename std::iterator_traits<ForwardIterType_sam>::value_type>::value));
		LZ_DEF_func_check_traits((std::is_same<type_real,
			typename std::iterator_traits<ForwardIterType_weight>::value_type>::value));
		//following is the value checking.
		if(stat_ker==type_this::val_stat_ker_NULL) return type_real(0);
		if(sam_begin==sam_end || weight_begin==weight_end) return type_real(0); //invalid parameters;
		ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
		ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
		while(iter_sam!=sam_end && iter_wt!=weight_end){
			if(iter_sam->rows()!=pt_ker->Get_dim() || iter_sam->cols()!=1 || (*iter_wt)<=(type_real)0){
				return type_real(0); //invalid parameters;
			}
			++iter_sam;
			++iter_wt;
		}
		if(iter_sam!=sam_end) return type_real(0);
		//the following calculates the result.
		type_real inte_samsam=0,inte_popsam=0; //the integrations of the kernel about sample and population;
		ForwardIterType_sam iter_sam0=sam_begin; //iterator to samples;
		ForwardIterType_weight iter_wt0=weight_begin; //iterator to weights;
		while(iter_sam0!=sam_end){
			inte_popsam+=((*pt_ker)(*iter_sam0)*(*iter_wt0));
			inte_samsam+=((*pt_ker)(*iter_sam0,*iter_sam0)*(*iter_wt0)*(*iter_wt0));
			iter_sam=iter_sam0;
			iter_wt=iter_wt0;
			++iter_sam;
			++iter_wt;
			while(iter_sam!=sam_end){
				inte_samsam+=((*pt_ker)(*iter_sam0,*iter_sam)*(*iter_wt0)*(*iter_wt)*type_real(2));
				++iter_sam;
				++iter_wt;
			}
			++iter_sam0;
			++iter_wt0;
		}
		inte_samsam=(*pt_ker)()-(type_real)2*inte_popsam+inte_samsam;
		return inte_samsam>(type_real)0 ? type_real(sqrt(inte_samsam)) : (type_real)0;
	}
	
	//Set:
	
	template<typename RealType>
	inline
	typename Discrepancy_canon_Leb_2_kernel<RealType>::type_stat
	Discrepancy_canon_Leb_2_kernel<RealType>::Set_kernel
	(type_kernel & kernel){
		if(stat_ker!=val_stat_ker_nonew && stat_ker!=val_stat_ker_NULL) delete pt_ker;
		pt_ker=std::addressof(kernel);
		stat_ker=type_this::val_stat_ker_nonew;
		return stat_ker;
	}
	
	template<typename RealType>
	inline
	typename Discrepancy_canon_Leb_2_kernel<RealType>::type_stat
	Discrepancy_canon_Leb_2_kernel<RealType>::Set_kernel
	(type_size dimension,type_stat kernel,type_real weight_projective){
		if(stat_ker!=val_stat_ker_nonew && stat_ker!=val_stat_ker_NULL) delete pt_ker;
		dimension= dimension>(type_size)0 ? dimension : (type_size)0;
		switch(kernel){
			case type_this::val_stat_ker_centered:
				pt_ker=new (std::nothrow)type_kernel_centered(dimension,weight_projective);
				break;
			case type_this::val_stat_ker_wrap:
				pt_ker=new (std::nothrow)type_kernel_wrap(dimension,weight_projective);
				break;
			case type_this::val_stat_ker_mixture:
				pt_ker=new (std::nothrow)type_kernel_mixture(dimension,weight_projective);
				break;
			default:
				pt_ker=NULL;
		}
		if(pt_ker==NULL) stat_ker=type_this::val_stat_ker_NULL;
		else stat_ker=kernel;
		return stat_ker;
	}
	
//}(class Discrepancy_canon_Leb_2_kernel)
/****************************************************************************************************/
} //namespace crit;
} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_stats_src_DOE_DOE_crit_uniformity_CPP
#endif //#if LZ_DEF_LZ_H_stats_DOE_crit_uniformity!=202105L