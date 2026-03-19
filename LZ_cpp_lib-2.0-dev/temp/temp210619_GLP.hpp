#ifndef LZ_DEF_temp210619
#define LZ_DEF_temp210619

#include<Eigen/QR>
#include<LZ_H/math/distribution.hpp>

namespace Liuze{
namespace stats{

/*
//TODO: iterator for Type_MatTemp: for entry, for row and for column;
template<typename Tp>
class EigenIterator: public std::iterator<std::random_access_iterator_tag,Tp>{
	public:
		
};
*/

//The lexicographical order for Type_MatTemp:
//(in the future it may be modified to be compatible with the std::lexicographical_compare in <algorithm>)
template<typename RealType=Type_Real>
Type_Bool Order_prec_lexicographical
(Type_MatTemp<RealType> const & x0,Type_MatTemp<RealType> const & x1){
	LZ_DEF_func_check_traits((std::is_arithmetic<RealType>::value));
	Type_Size nrow=std::min(x0.rows(),x1.rows()),ncol=std::min(x0.cols(),x1.cols());
	if(nrow==0 || ncol==0) return (Type_Bool)false;
	Type_Size i,j;
	for(i=0;i<nrow;++i){
		for(j=0;j<ncol;++j){
			if(x0(i,j)<x1(i,j)) return (Type_Bool)true;
			if(x0(i,j)>x1(i,j)) return (Type_Bool)false;
		}
	}
	return (Type_Bool)false;
}

} //namespace stats;
} //namespace Liuze;

namespace Liuze{
namespace stats{

//Cal: the L^(+infinity) distance between the empirical CDF of some samples and a given CDF:
//this is the version for the multi-dimensional case:
//parameters:{
//sam_begin, sam_end: the iterators for providing the samples;
//weight_begin, weight_end: the iterators for providing the weights of the samples;
//func_distr: the target distribution function;
//CDF_def: 0: F(x)=P(X<x); 1: F(x)=P(X<=x);
//pre: the precision;}
template<typename ForwardIterType_sam,typename ForwardIterType_weight,typename FuncType>
typename ForwardIterType_weight::value_type
Cal_discrepancy_Leb_infin_multi
(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end,
FuncType & func_distr,Type_Stat const CDF_def=0,
typename ForwardIterType_weight::value_type pre=-1){
	typedef typename ForwardIterType_weight::value_type type_real;
	typedef typename ForwardIterType_sam::value_type type_matrix;
	//following is the type checking.
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::result_of<typename std::decay<FuncType>::type(type_matrix)>::type,
		type_real>::value));
	//following is the value checking.
	if(sam_begin==sam_end || sam_begin->rows()==0 || weight_begin==weight_end){
		return (type_real)(-1); //invalid parameters;
	}
	if(pre<=(type_real)0) pre=(type_real)LZ_DEF_const_default_precision; //the precision;
	Type_Size dim=sam_begin->rows(),size_sam=0;
	Type_Size i_dim,i_knot; //id vars;
	ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
	ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
	while(iter_sam!=sam_end && iter_wt!=weight_end){
		if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
			return (type_real)(-1); //invalid parameters;
		}
		++size_sam;
		++iter_sam;
		++iter_wt;
	}
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
				if(tr0<pre) pre=tr0; //re-tuning the precision;
			}
		}
		//the following determines the approximate positive infinity of this dimension.
		tr0=(type_real)5*std::max((type_real)1,pre);
		if(knots[i_dim].size()>(Type_Size)1){
			tr0*=(knots[i_dim][knots[i_dim].size()-1]-knots[i_dim][0]);
		} else {
			tr0*=std::max(Liuze::math::abs(knots[i_dim][0]),(type_real)1);
		}
		tr0+=(knots[i_dim][knots[i_dim].size()-1]);
		knots[i_dim].push_back(tr0);
	}
	//this part deals with the upper and lower values of the cdf.
	type_matrix disturb=type_matrix::Ones(dim,1)*(CDF_def==(Type_Stat)0 ? pre : -pre);
	type_real ecdf0,ecdf1,cdf0,cdf1; //store the values of ecdf and cdf under 2 definitions;
	//the following calculates the result.
	type_real res=0; //the return result;
	type_matrix knot_pos(dim,1); //one knot position;
	//indicating the knot is a key position (sam: at a sample; crd: the coordinates):
	std::vector<bool> is_key_pos_crd(dim),is_key_pos_crd_sam(dim),is_positive_infty(dim);
	Type_Bool is_key_pos_sam;
	//Following: the loop for each knot point.
	Type_ArrayTemp<Type_Size> loop_knot(dim,(Type_Size)0); //the loop var;
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
					cdf0=func_distr(knot_pos+disturb);
					cdf1=func_distr(knot_pos);
					break;
				case (Type_Stat)0:
				default:
					cdf0=func_distr(knot_pos);
					cdf1=func_distr(knot_pos+disturb);
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
		while((++i_dim)<dim) loop_knot[i_dim]=(Type_Size)0;
	}
	GTS_loop_end_knot:
	return res;
} //Cal_Discrepancy_Leb_infin_multi;

//Cal: the wrap-around L^2 discrepancy with respect to the uniform distribution for canonical designs:
template<typename RealType=Type_Real>
class Cal_discrepancy_canon_Wrap_2{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//Type: the real number:
		typedef RealType type_real;
		
		Cal_discrepancy_canon_Wrap_2():
		weight_projective(1){
		}
		Cal_discrepancy_canon_Wrap_2(type_real const & weight_projective_uniform):
		weight_projective(std::max(weight_projective_uniform,(type_real)0)){
		}
		
		void Set_weight_projective(){
			weight_projective=(type_real)1;
		}
		void Set_weight_projective(type_real const & weight_projective_uniform){
			weight_projective=std::max(weight_projective_uniform,(type_real)0);
		}
		
		//fun: Cal: L^2 wrap-around discrepancy for canonical designs:
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
			typedef typename ForwardIterType_sam::value_type type_matrix;
			//following is the type checking.
			LZ_DEF_func_check_traits((std::is_same<typename ForwardIterType_weight::value_type,
				type_real>::value));
			LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
			//following is the value checking.
			if(sam_begin==sam_end || sam_begin->rows()==0 || weight_begin==weight_end){
				return (type_real)(-1); //invalid parameters;
			}
			Type_Size dim=sam_begin->rows();
			Type_Size i_dim; //id vars;
			ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
			ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
			while(iter_sam!=sam_end && iter_wt!=weight_end){
				if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
					return (type_real)(-1); //invalid parameters;
				}
				++iter_sam;
				++iter_wt;
			}
			//the following calculates the result.
			//cr1=1/2+weight_projective, cr1 is 3/2 by default:
			type_real res=0,cr0=1,cr1=(type_real)1/(type_real)2+weight_projective,cr2=1;
			type_real tr0=(type_real)1/(type_real)3+weight_projective,tr1; //temp var;
			//cr0=(1/3+weight_projective)^(dim), cr2=(1/2+weight_projective)^(dim):
			for(i_dim=0;i_dim<dim;++i_dim){
				cr0*=tr0;
				cr2*=cr1;
			}
			ForwardIterType_sam iter_sam0=sam_begin; //iterator to samples;
			ForwardIterType_weight iter_wt0=weight_begin; //iterator to weights;
			iter_sam=iter_sam0;
			iter_wt=iter_wt0;
			while(iter_sam0!=sam_end){
				res+=(cr2*(*iter_wt0)*(*iter_wt0));
				++iter_sam;
				++iter_wt;
				while(iter_sam!=sam_end){
					tr0=(type_real)1;
					for(i_dim=0;i_dim<dim;++i_dim){
						tr1=math::abs((*iter_sam0)(i_dim,0)-(*iter_sam)(i_dim,0));
						tr0*=(cr1-tr1+tr1*tr1);
					}
					tr0*=((*iter_wt0)*(*iter_wt));
					res+=(tr0+tr0);
					++iter_sam;
					++iter_wt;
				}
				iter_sam=++iter_sam0;
				iter_wt=++iter_wt0;
			}
			res=sqrt(res-cr0);
			return res;
		}
	protected:
		type_real weight_projective; //the weight of the projective uniformity, default is 1;
}; //class Cal_discrepancy_canon_Wrap_2;

//Cal: the approximated L^infinity discrepancy, based on p97 of Hua L.-K. and Wang Y. (1981):
template<typename RealType=Type_Real>
class Cal_discrepancy_canon_approxLambda{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//Type: the real number:
		typedef RealType type_real;
		
		Cal_discrepancy_canon_approxLambda()=default;
		
		//fun: Cal: the approximated Lambda function for canonical designs:
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
			typedef typename ForwardIterType_sam::value_type type_matrix;
			//following is the type checking.
			LZ_DEF_func_check_traits((std::is_same<typename ForwardIterType_weight::value_type,
				type_real>::value));
			LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
			//following is the value checking.
			if(sam_begin==sam_end || sam_begin->rows()==0 || weight_begin==weight_end){
				return (type_real)(-1); //invalid parameters;
			}
			Type_Size dim=sam_begin->rows();
			Type_Size i_dim; //id vars;
			ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
			ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
			while(iter_sam!=sam_end && iter_wt!=weight_end){
				if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
					return (type_real)(-1); //invalid parameters;
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
				for(i_dim=(Type_Size)0;i_dim<dim;++i_dim){
					tr0*=(c1-c2_pi*log(c2*sin(cpi*((*iter_sam)(i_dim,0)))));
				}
				res+=(tr0*(*iter_wt));
				++iter_wt;
			}
			return res;
		}
}; //class Cal_discrepancy_canon_approxLambda;

template<typename RealType=Type_Real>
class Cal_orthogonality_canon{
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//Type: the real number:
		typedef RealType type_real;
		//Type: the state:
		typedef Type_Stat type_stat;
		
		//static const values:
		static type_stat constexpr val_ortho_type_max=0; //the max correlation;
		static type_stat constexpr val_ortho_type_mean=1; //the average correlation;
		static type_stat constexpr val_ortho_type_median=2; //the median of the correlations;
		
		Cal_orthogonality_canon():
		stat_ortho_type(val_ortho_type_max),list_corr(){
		}
		Cal_orthogonality_canon(type_stat const & ortho_type):
		stat_ortho_type(ortho_type),list_corr(){
		}
		
		void Set_ortho_type(){
			stat_ortho_type=val_ortho_type_max;
		}
		void Set_ortho_type(type_stat const & ortho_type){
			stat_ortho_type=ortho_type;
		}
		
		//fun: Cal: the correlations between each pair of columns:
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_stat Cal_corr
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
			typedef typename ForwardIterType_sam::value_type type_matrix;
			//following is the type checking.
			LZ_DEF_func_check_traits((std::is_same<typename ForwardIterType_weight::value_type,
				type_real>::value));
			LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
			//following is the value checking.
			Type_Size dim=sam_begin->rows();
			if(dim==(Type_Size)0) return (type_stat)1;
			if(dim==(Type_Size)1){
				list_corr.resize(0);
				return (type_stat)0;
			}
			ForwardIterType_sam iter_sam=sam_begin; //iterator to samples;
			ForwardIterType_weight iter_wt=weight_begin; //iterator to weights;
			while(iter_sam!=sam_end && iter_wt!=weight_end){
				if(iter_sam->rows()<dim || iter_sam->cols()==0 || (*iter_wt)<=(type_real)0){
					return (type_stat)1; //invalid parameters;
				}
				++iter_sam;
				++iter_wt;
			}
			//the following calculates this->list_corr.
			list_corr.assign(dim*(dim-(Type_Size)1)/(Type_Size)2,(type_real)0);
			Type_ArrayTemp<type_real> list_len(dim,(type_real)0); //the L^2 lengths of the columns;
			Type_Size i_dim,i0_dim,i1_dim; //id vars;
			type_real tr0,tr1; //temp real;
			iter_sam=sam_begin;
			iter_wt=weight_begin;
			while(iter_sam!=sam_end && iter_wt!=weight_end){
				for(i0_dim=i_dim=0;i0_dim<dim;++i0_dim){
					tr1=(*iter_sam)(i0_dim,0)-(type_real)0.5;
					tr0=(*iter_wt)*(*iter_wt)*tr1;
					list_len[i0_dim]+=(tr0*tr1);
					for(i1_dim=i0_dim+1;i1_dim<dim;++i1_dim){
						list_corr[i_dim]+=(tr0*((*iter_sam)(i1_dim,0)-(type_real)0.5));
						++i_dim;
					}
				}
				++iter_sam;
				++iter_wt;
			}
			for(i0_dim=0;i0_dim<dim;++i0_dim) list_len[i0_dim]=sqrt(list_len[i0_dim]);
			for(i0_dim=i_dim=0;i0_dim<dim;++i0_dim){
				for(i1_dim=i0_dim+1;i1_dim<dim;++i1_dim){
					list_corr[i_dim]=math::abs(list_corr[i_dim])/(list_len[i0_dim]*list_len[i1_dim]);
					++i_dim;
				}
			}
			return (type_stat)0;
		}
		
		//fun: Cal: orthogonality for canonical designs:
		template<typename ForwardIterType_sam,typename ForwardIterType_weight>
		type_real operator ()
		(ForwardIterType_sam const & sam_begin,ForwardIterType_sam const & sam_end,
		ForwardIterType_weight const & weight_begin,ForwardIterType_weight const & weight_end) const {
			if(Cal_corr(sam_begin,sam_end,weight_begin,weight_end)!=(type_stat)0) return (type_real)(-1);
			if(list_corr.size()==(Type_Size)0) return (type_real)0;
			if(stat_ortho_type==val_ortho_type_max){
				return *(std::max_element(list_corr.begin(),list_corr.end()));
			}
			if(stat_ortho_type==val_ortho_type_mean){
				return std::accumulate(list_corr.begin(),list_corr.end(),(type_real)0)/
					(type_real)(list_corr.size());
			}
			if(stat_ortho_type==val_ortho_type_median){
				std::sort(list_corr.begin(),list_corr.end());
				return list_corr[list_corr.size()/(Type_Size)2];
			}
			return (type_real)(-1);
		}
	protected:
		type_stat stat_ortho_type; //the type of the representor for the distribution of list_corr;
		Type_ArrayTemp<type_real> mutable list_corr; //the correlations between each pair of columns;
}; //class Cal_orthogonality_canon;

template<typename IntType=Type_UInt,typename RealType=Type_Real>
class Iterator_design_canon_sequential{
	LZ_DEF_func_check_traits((std::is_integral<IntType>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//Type: this class:
		typedef Iterator_design_canon_sequential<IntType,RealType> type_this;
		//Type: the base class:
		typedef void type_base;
		//Type: the integer:
		typedef IntType type_int;
		//Type: the real number:
		typedef RealType type_real;
		//Type: the matrix:
		typedef Type_MatTemp<type_real> type_matrix;
		//Type: the state:
		typedef Type_Stat type_stat;
		//Type: its iterator:
		typedef type_this iterator;
		//Types: for iterator:
		typedef std::iterator<std::forward_iterator_tag,type_matrix const> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//static values:
		static type_stat constexpr val_stat_deref=0;
		static type_stat constexpr val_stat_end=1;
		static type_stat constexpr val_stat_null=10;
		
		//Constructor:
		Iterator_design_canon_sequential():
		dim(0),stat(this->val_stat_null),point(0,0){
		}
		Iterator_design_canon_sequential
		(type_int const & dimension){
			if(dimension<=(type_int)0){
				stat=val_stat_null;
				dim=(type_int)0;
				point.resize(0,0);
			} else {
				stat=val_stat_deref;
				dim=dimension;
				point.resize(dim,1);
			}
		}
		Iterator_design_canon_sequential(type_this const & iter):
		dim(iter.dim),point(iter.point),stat(iter.stat){
		}
		Iterator_design_canon_sequential(type_this && iter):
		dim(iter.dim),point(iter.point),stat(iter.stat){
		}
		
		//Destructor:
		virtual ~Iterator_design_canon_sequential()=default;
		
		//operator:
		type_this& operator =(type_this const & iter){
			if(&iter==this) return *this;
			dim=iter.dim;
			point=iter.point;
			stat=iter.stat;
			return *this;
		}
		type_this& operator =(type_this && iter){
			dim=iter.dim;
			point=iter.point;
			stat=iter.stat;
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return stat!=iter.stat ? false :
				(stat==val_stat_null ? true :
				(dim!=iter.dim ? false : (stat==val_stat_end ? true : point==iter.point)));
		}
		virtual value_type operator *() const final {
			return point;
		}
		virtual value_type * operator ->() const final {
			return std::addressof(point);
		}
		virtual type_this& operator ++()=0;
		
		//other functions:
		virtual type_stat Set_begin()=0;
		
		type_int Get_dim() const {
			return dim;
		}
		type_stat Get_stat() const {
			return stat;
		}
	protected:
		type_int dim;
		type_matrix point; //the current design point;
		type_stat stat; //the state;
}; //class Iterator_design_canon_sequential;

template<typename IntType=Type_UInt,typename RealType=Type_Real,
	typename RandUIntType=typename Liuze::math::Type_RandEngine::result_type>
class Iterator_design_canon_RandomUniform:
public Iterator_design_canon_sequential<IntType,RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_RandomUniform<IntType,RealType,RandUIntType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<IntType,RealType> type_base;
		//Type: the integer:
		typedef typename type_base::type_int type_int;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the random engine:
		typedef Liuze::math::RandomEngine<RandUIntType> type_rand_eng;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef typename type_base::iterator_category iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_RandomUniform():
		type_base(),engine(NULL),distr_U((type_real)0,(type_real)1){
		}
		Iterator_design_canon_RandomUniform(type_rand_eng & rand_engine):
		type_base(),engine(&rand_engine),distr_U((type_real)0,(type_real)1){
		}
		Iterator_design_canon_RandomUniform
		(type_int const & dimension,type_rand_eng & rand_engine):
		type_base(dimension),engine(&rand_engine),distr_U((type_real)0,(type_real)1){
			this->Cal_renew();
		}
		Iterator_design_canon_RandomUniform(type_this const & iter):
		engine(iter.engine),distr_U((type_real)0,(type_real)1),type_base((type_base const &)iter){
		}
		Iterator_design_canon_RandomUniform(type_this && iter):
		engine(iter.engine),distr_U((type_real)0,(type_real)1),type_base((type_base &&)iter){
		}
		
		//Destructor:
		virtual ~Iterator_design_canon_RandomUniform()=default;
		
		//operator:
		type_this& operator =(type_this const & iter){
			if(&iter==this) return *this;
			((type_base*)this)->operator =((type_base const &)iter);
			engine=iter.engine;
			return *this;
		}
		type_this& operator =(type_this && iter){
			engine=iter.engine;
			((type_base*)this)->operator =((type_base &&)iter);
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return (*((type_base*)this))==(type_base const &)iter && *engine==*(iter.engine);
		}
		virtual type_this& operator ++(){
			this->Cal_renew();
			return *this;
		}
		
		//other functions:
		virtual type_stat Set_begin(){
			if(this->stat!=this->val_stat_null) this->stat=this->val_stat_deref;
			this->Cal_renew();
			return this->stat;
		}
		type_stat Set_begin(type_int const & dimension){
			if(engine==NULL || dimension<=(type_int)0){
				this->stat=this->val_stat_null;
				this->point.resize(0,0);
				return this->stat;
			}
			this->stat=this->val_stat_deref;
			if(dimension!=this->dim){
				this->dim=dimension;
				this->point.resize(this->dim,1);
			}
			this->Cal_renew();
			return this->stat;
		}
		
		void Set_rand_engine(type_rand_eng & rand_engine){
			engine=&rand_engine;
		}
	protected:
		std::uniform_real_distribution<type_real> distr_U; //the uniform random variable generator;
		type_rand_eng * engine; //the random engine;
		
		void Cal_renew(){
			if(this->stat==this->val_stat_deref){
				for(type_int i_dim=0;i_dim<this->dim;++i_dim) this->point(i_dim,0)=distr_U(*engine);
			}
		}
}; //class Iterator_design_canon_RandomUniform;

template<typename IntType=Type_UInt,typename RealType=Type_Real>
class Iterator_design_canon_GoodPoint_cyclotomic:
public Iterator_design_canon_sequential<IntType,RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_GoodPoint_cyclotomic<IntType,RealType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<IntType,RealType> type_base;
		//Type: the integer:
		typedef typename type_base::type_int type_int;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef typename type_base::iterator_category iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_GoodPoint_cyclotomic():
		type_base(){
		}
		Iterator_design_canon_GoodPoint_cyclotomic
		(type_int const & dimension,type_int const & prime_number):
		type_base(){
			if(dimension>(type_int)0 && prime_number>=(type_int)2*dimension+(type_int)3){
				this->stat=this->val_stat_deref;
				this->dim=dimension;
				prime=prime_number;
				gp.resize(this->dim,1);
				type_real angle0=(type_real)LZ_DEF_const_math_2pi/(type_real)prime;
				for(type_int i_dim=0;i_dim<this->dim;++i_dim){
					gp(i_dim,0)=(type_real)2*cos(angle0*(type_real)(i_dim+1));
					gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
				}
				this->point=gp;
			}
		}
		Iterator_design_canon_GoodPoint_cyclotomic(type_this const & iter):
		prime(iter.prime),gp(iter.gp),type_base((type_base const &)iter){
		}
		Iterator_design_canon_GoodPoint_cyclotomic(type_this && iter):
		prime(iter.prime),gp(iter.gp),type_base((type_base &&)iter){
		}
		
		//Destructor:
		virtual ~Iterator_design_canon_GoodPoint_cyclotomic()=default;
		
		//operator:
		type_this& operator =(type_this const & iter){
			if(&iter==this) return *this;
			((type_base*)this)->operator =((type_base const &)iter);
			prime=iter.prime;
			gp=iter.gp;
			return *this;
		}
		type_this& operator =(type_this && iter){
			prime=iter.prime;
			gp=iter.gp;
			((type_base*)this)->operator =((type_base &&)iter);
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return (*((type_base*)this))!=(type_base const &)iter ? false :
				(this->stat==this->val_stat_null || prime==iter.prime);
		}
		virtual type_this& operator ++(){
			if(this->stat!=this->val_stat_deref) return *this;
			for(type_int i_dim=0;i_dim<this->dim;++i_dim){
				this->point(i_dim,0)+=gp(i_dim,0);
				if(this->point(i_dim,0)>=(type_real)1) this->point(i_dim,0)-=(type_real)1;
			}
			return *this;
		}
		
		//other functions:
		virtual type_stat Set_begin(){
			if(this->stat==this->val_stat_null) return this->stat;
			this->stat=this->val_stat_deref;
			this->point=gp;
			return this->stat;
		}
		type_stat Set_begin(type_int const & dimension,type_int const & prime_number){
			if(dimension<=(type_int)0 || prime_number<(type_int)2*dimension+(type_int)3){
				this->stat=this->val_stat_null;
				this->point.resize(0,0);
				return this->stat;
			}
			this->stat=this->val_stat_deref;
			if(dimension!=this->dim || prime_number!=prime){
				this->dim=dimension;
				prime=prime_number;
				gp.resize(this->dim,1);
				type_real angle0=(type_real)LZ_DEF_const_math_2pi/(type_real)prime;
				for(type_int i_dim=0;i_dim<this->dim;++i_dim){
					gp(i_dim,0)=(type_real)2*cos(angle0*(type_real)(i_dim+1));
					gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
				}
			}
			this->point=gp;
			return this->stat;
		}
	protected:
		type_int prime;
		type_matrix gp; //the good point;
}; //class Iterator_design_canon_GoodPoint_cyclotomic;

template<typename IntType=Type_UInt,typename RealType=Type_Real>
class Iterator_design_canon_GoodPoint_power:
public Iterator_design_canon_sequential<IntType,RealType>{
	public:
		//Type: this class:
		typedef Iterator_design_canon_GoodPoint_power<IntType,RealType> type_this;
		//Type: the base class:
		typedef Iterator_design_canon_sequential<IntType,RealType> type_base;
		//Type: the integer:
		typedef typename type_base::type_int type_int;
		//Type: the real number:
		typedef typename type_base::type_real type_real;
		//Type: the matrix:
		typedef typename type_base::type_matrix type_matrix;
		//Type: the state:
		typedef typename type_base::type_stat type_stat;
		//Type: its iterator:
		typedef typename type_base::iterator iterator;
		//Types: for iterator:
		typedef typename type_base::iterator_category iterator_category;
		typedef typename type_base::value_type value_type;
		typedef typename type_base::difference_type difference_type;
		typedef typename type_base::pointer pointer;
		typedef typename type_base::reference reference;
		
		//Constructor:
		Iterator_design_canon_GoodPoint_power():
		type_base(){
		}
		Iterator_design_canon_GoodPoint_power
		(type_int const & dimension,type_int const prime_number=2,type_stat const modify=0):
		//When dimension=d, prime_number=p, let q=p^(1/(d+1)), then gp=(q^1,q^2,...,q^d);
		//but when modify=1, gp=(a1,b1,a2,b2,...) with (a1,a2,...,b1,b2,...)=(q^1,q^2,...,q^d),
		//see J. G. Liao (1998) on Journal of Computational and Graphical Statistics.
		type_base(){
			if(dimension>(type_int)0 && prime_number>=(type_int)2){
				this->stat=this->val_stat_deref;
				this->dim=dimension;
				prime=prime_number;
				type_int i_dim;
				gp.resize(this->dim,1);
				gp(0,0)=type_real(pow((type_real)prime,(type_real)1/(type_real(this->dim)+(type_real)1)));
				if(modify==(type_stat)1){
					type_real tr=gp(0,0);
					type_int i_dim_0;
					for(i_dim_0=2;i_dim_0>=1;--i_dim_0){
						i_dim=i_dim_0;
						while(i_dim<this->dim){
							gp(i_dim,0)=(tr*=gp(0,0));
							i_dim+=2;
						}
					}
				} else {
					for(i_dim=1;i_dim<this->dim;++i_dim) gp(i_dim,0)=gp(i_dim-1,0)*gp(0,0);
				}
				for(i_dim=0;i_dim<this->dim;++i_dim) gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
				this->point=gp;
			}
		}
		Iterator_design_canon_GoodPoint_power(type_this const & iter):
		prime(iter.prime),gp(iter.gp),type_base((type_base const &)iter){
		}
		Iterator_design_canon_GoodPoint_power(type_this && iter):
		prime(iter.prime),gp(iter.gp),type_base((type_base &&)iter){
		}
		
		//Destructor:
		virtual ~Iterator_design_canon_GoodPoint_power()=default;
		
		//operator:
		type_this& operator =(type_this const & iter){
			if(&iter==this) return *this;
			((type_base*)this)->operator =((type_base const &)iter);
			prime=iter.prime;
			gp=iter.gp;
			return *this;
		}
		type_this& operator =(type_this && iter){
			prime=iter.prime;
			gp=iter.gp;
			((type_base*)this)->operator =((type_base &&)iter);
			return *this;
		}
		Type_Bool operator ==(type_this const & iter) const {
			return (*((type_base*)this))!=(type_base const &)iter ? false :
				(this->stat==this->val_stat_null || prime==iter.prime);
		}
		virtual type_this& operator ++(){
			if(this->stat!=this->val_stat_deref) return *this;
			for(type_int i_dim=0;i_dim<this->dim;++i_dim){
				this->point(i_dim,0)+=gp(i_dim,0);
				if(this->point(i_dim,0)>=(type_real)1) this->point(i_dim,0)-=(type_real)1;
			}
			return *this;
		}
		
		//other functions:
		virtual type_stat Set_begin(){
			if(this->stat==this->val_stat_null) return this->stat;
			this->stat=this->val_stat_deref;
			this->point=gp;
			return this->stat;
		}
		type_stat Set_begin
		(type_int const & dimension,type_int const prime_number=2,type_stat const modify=0){
			if(dimension<=(type_int)0 || prime_number<(type_int)2){
				this->stat=this->val_stat_null;
				this->point.resize(0,0);
				return this->stat;
			}
			this->stat=this->val_stat_deref;
			if(dimension!=this->dim || prime_number!=prime){
				this->dim=dimension;
				prime=prime_number;
				type_int i_dim;
				gp.resize(this->dim,1);
				gp(0,0)=type_real(pow((type_real)prime,(type_real)1/(type_real(this->dim)+(type_real)1)));
				if(modify==(type_stat)1){
					type_real tr=gp(0,0);
					type_int i_dim_0;
					for(i_dim_0=2;i_dim_0>=1;--i_dim_0){
						i_dim=i_dim_0;
						while(i_dim<this->dim){
							gp(i_dim,0)=(tr*=gp(0,0));
							i_dim+=2;
						}
					}
				} else {
					for(i_dim=1;i_dim<this->dim;++i_dim) gp(i_dim,0)=gp(i_dim-1,0)*gp(0,0);
				}
				for(i_dim=0;i_dim<this->dim;++i_dim) gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
			}
			this->point=gp;
			return this->stat;
		}
	protected:
		type_int prime;
		type_matrix gp; //the good point;
}; //class Iterator_design_canon_GoodPoint_power;

//fun: Get_design_canon_GoodPoint_cyclotomic:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_GoodPoint_cyclotomic
(IntType const & size,IntType const & dim,IntType prime,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	if(size<=(type_int)0 || dim<=(type_int)0) return result_begin;
	type_int i_size,i_dim;
	i_dim=(type_int)2*dim+(type_int)3;
	if(prime<i_dim) prime=i_dim; //prime should be a prime no less than 2*dim+3;
	//checking ends here.
	type_matrix gp(dim,1); //the basic good point vector;
	type_real tr=(type_real)LZ_DEF_const_math_2pi/(type_real)prime;
	for(i_dim=0;i_dim<dim;++i_dim){
		gp(i_dim,0)=(type_real)2*cos(tr*(type_real)(i_dim+1));
		gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
	}
	type_matrix point_pre=type_matrix::Zero(dim,1);
	for(i_size=0;i_size<size;++i_size){
		(*result_begin).resize(dim,1);
		for(i_dim=0;i_dim<dim;++i_dim){
			(*result_begin)(i_dim,0)=point_pre(i_dim,0)+gp(i_dim,0);
			if((*result_begin)(i_dim,0)>=(type_real)1) (*result_begin)(i_dim,0)-=(type_real)1;
			point_pre(i_dim,0)=(*result_begin)(i_dim,0);
		}
		++result_begin;
	}
	return result_begin;
} //(fun: Get_design_canon_GoodPoint_cyclotomic)

//fun: Get_design_canon_GoodPoint_power:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_GoodPoint_power
(IntType const & size,IntType const & dim,IntType prime,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	if(size<=(type_int)0 || dim<=(type_int)0) return result_begin;
	type_int i_size,i_dim;
	if(prime<=(type_int)1) prime=(type_int)2;
	//checking ends here.
	type_matrix gp(dim,1); //the basic good point vector;
	gp(0,0)=type_real(pow((type_real)prime,(type_real)1/((type_real)dim+(type_real)1)));
	for(i_dim=1;i_dim<dim;++i_dim) gp(i_dim,0)=gp(i_dim-1,0)*gp(0,0);
	for(i_dim=0;i_dim<dim;++i_dim) gp(i_dim,0)-=type_real(floor(gp(i_dim,0)));
	type_matrix point_pre=type_matrix::Zero(dim,1);
	for(i_size=0;i_size<size;++i_size){
		(*result_begin).resize(dim,1);
		for(i_dim=0;i_dim<dim;++i_dim){
			(*result_begin)(i_dim,0)=point_pre(i_dim,0)+gp(i_dim,0);
			if((*result_begin)(i_dim,0)>=(type_real)1) (*result_begin)(i_dim,0)-=(type_real)1;
			point_pre(i_dim,0)=(*result_begin)(i_dim,0);
		}
		++result_begin;
	}
	return result_begin;
} //(fun: Get_design_canon_GoodPoint_power)

//fun: Get_design_factor_GoodLatticePoint:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<IntType> >::iterator>
OutIterType Get_design_factor_GoodLatticePoint
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_int>,type_matrix>::value));
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return result_begin;
	//checking ends.
	type_matrix point=type_matrix::Zero(dim,1);
	type_int i_dim,i_level=0;
	while(true){
		*result_begin=point;
		++result_begin;
		++i_level;
		if(i_level==num_level) return result_begin;
		for(i_dim=0;i_dim<dim;++i_dim){
			point(i_dim,0)=(point(i_dim,0)+vector_gen(i_dim,0))%num_level;
		}
	}
} //(fun: Get_design_factor_GoodLatticePoint)
//fun: Get_design_canon_GoodLatticePoint:
//applying the transformation: x |-> (x+(0.5,...,0.5))/num_level
//for all points in the corresponding factor good lattice point design:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_GoodLatticePoint
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return result_begin;
	//checking ends here.
	type_matrix point(dim,1);
	Type_ArrayTemp<type_real> crd(num_level);
	Type_ArrayTemp<type_int> point_id(dim,(type_int)0); //the point in the factor design;
	type_int i_dim,i_level=0;
	//calculating the 1-dim coordinates:
	for(i_level=0;i_level<num_level;++i_level){
		crd[i_level]=((type_real)i_level+(type_real)(0.5))/((type_real)num_level);
	}
	//transformation from factor design to canonical design:
	auto fun_Get_point=[&i_dim,&dim,&point,&crd,&point_id]()->type_matrix&{
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
		if(i_level==num_level) return result_begin;
		for(i_dim=0;i_dim<dim;++i_dim){
			point_id[i_dim]=(point_id[i_dim]+vector_gen(i_dim,0))%num_level;
		}
	}
} //(fun: Get_design_canon_GoodLatticePoint)

//fun: Get_design_factor_GoodLatticePoint_LeaveOneOut:
//excluding the point (0,...,0) from the GLP(num_level+1,vector_gen)
//and then subtracting (1,...,1) from the rested points:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<IntType> >::iterator>
OutIterType Get_design_factor_GoodLatticePoint_LeaveOneOut
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_int>,type_matrix>::value));
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return result_begin;
	//checking ends.
	type_int num_level_full=num_level+(type_int)1;
	type_matrix point=type_matrix::Zero(dim,1),point_shift=type_matrix::Ones(dim,1);
	type_int i_dim,i_level;
	for(i_level=(type_int)0;i_level<num_level;++i_level){
		for(i_dim=0;i_dim<dim;++i_dim){
			point(i_dim,0)=(point(i_dim,0)+vector_gen(i_dim,0))%num_level_full;
		}
		*result_begin=point-point_shift;
		++result_begin;
	}
	return result_begin;
} //(fun: Get_design_factor_GoodLatticePoint_LeaveOneOut)
//fun: Get_design_canon_GoodLatticePoint_LeaveOneOut:
//applying the transformation: x |-> (x+(0.5,...,0.5))/num_level
//for all points in the corresponding leave-one-out factor good lattice point design:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_GoodLatticePoint_LeaveOneOut
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return result_begin;
	//checking ends here.
	type_int num_level_full=num_level+(type_int)1;
	type_matrix point(dim,1);
	Type_ArrayTemp<type_real> crd(num_level);
	Type_ArrayTemp<type_int> point_id(dim,(type_int)0); //the point in the factor design;
	type_int i_dim,i_level=0;
	//calculating the 1-dim coordinates:
	for(i_level=0;i_level<num_level;++i_level){
		crd[i_level]=((type_real)i_level+(type_real)(0.5))/((type_real)num_level);
	}
	//transformation from the full factor design to leave-one-out canonical design:
	auto fun_Get_point=[&i_dim,&dim,&point,&crd,&point_id]()->type_matrix&{
			for(i_dim=(type_int)0;i_dim<dim;++i_dim){
				point(i_dim,0)=crd[point_id[i_dim]-(type_int)1];
			}
			return point;
		};
	for(i_level=(type_int)0;i_level<num_level;++i_level){
		for(i_dim=0;i_dim<dim;++i_dim){
			point_id[i_dim]=(point_id[i_dim]+vector_gen(i_dim,0))%num_level_full;
		}
		*result_begin=fun_Get_point();
		++result_begin;
	}
	return result_begin;
} //(fun: Get_design_canon_GoodLatticePoint_LeaveOneOut)

//class: cls_Get_design_canon_GoodLatticePoint_optimal:
template<typename IntType=Type_UInt,typename RealType=Type_Real>
class cls_Get_design_GoodLatticePoint{
	LZ_DEF_func_check_traits((std::is_integral<IntType>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<RealType>::value));
	public:
		//type: this class:
		typedef cls_Get_design_GoodLatticePoint<IntType,RealType> type_this;
		//type: integer:
		typedef IntType type_int;
		//type: real number:
		typedef RealType type_real;
		//type: point for factorial design:
		typedef Type_MatTemp<type_int> type_pt_fac;
		//type: point for canonical design:
		typedef Type_MatTemp<type_real> type_pt_can;
		
		//Constructor:
		cls_Get_design_GoodLatticePoint():
		size(0){
			dim_max=(type_int)0;
			stat=(Type_Stat)1;
		}
		cls_Get_design_GoodLatticePoint(type_int const & num_size):
		size(num_size){
			if(size<=(type_int)0){
				size=(type_int)0;
				dim_max=(type_int)0;
				stat=(Type_Stat)1;
			} else {
				//the following calculates dim_max and unit.
				dim_max=(type_int)1;
				Type_ArrayTemp<type_int> unit_init(size-1);
				unit_init[0]=(type_int)1;
				type_int i;
				for(i=2;i<size;++i){
					if(math::Is_coprime(size,i)) unit_init[dim_max++]=i;
				}
				unit.assign(unit_init.begin(),unit_init.begin()+dim_max);
				//the following calculates crd_can.
				crd_can.resize(size);
				type_real tr0=size,tr1=(type_real)1/((type_real)(2*size));
				for(i=0;i<size;++i) crd_can[i]=(type_real)i/tr0+tr1;
			}
		}
		
		//Cal: find the best GLP over all vector_gen:
		template<typename CritType>
		Type_Stat Cal_best(type_int const & num_dim,CritType & fun_crit_can){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			dim=num_dim;
			if(dim<=(type_int)0 || dim>dim_max){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			typedef typename std::make_unsigned<type_int>::type type_uint;
			vec_gen.resize(dim,1);
			Type_Bool tag_first_des=true;
			type_real v_crit; //inner value of criterion;
			Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
			math::Comb_enum_comb<type_uint> iter_comb(dim_max,dim);
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
			type_int i_dim,i_size;
			for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
			while((*iter_comb)[0]==(type_uint)0){
				for(i_size=0;i_size<size;++i_size){
					for(i_dim=0;i_dim<dim;++i_dim){
						des[i_size](i_dim,0)=crd_can[(unit[(*iter_comb)[i_dim]]*i_size)%size];
					}
				}
				v_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
				if(v_crit<val_crit || tag_first_des){
					tag_first_des=false;
					val_crit=v_crit;
					for(i_dim=0;i_dim<dim;++i_dim) vec_gen(i_dim,0)=unit[(*iter_comb)[i_dim]];
				}
				++iter_comb;
			}
			stat=(Type_Stat)0;
			return stat;
		}
		//Cal: find the best GLP over all power vector_gen:
		template<typename CritType>
		Type_Stat Cal_best_power(type_int const & num_dim,CritType & fun_crit_can){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			dim=num_dim;
			if(dim<=(type_int)0 || dim>dim_max){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			Type_Bool tag_first_des=true;
			type_real v_crit; //inner value of criterion;
			Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
			type_pt_fac v_gen(dim,1); //inner vec_gen;
			Type_ArrayTemp<type_int> val_v_gen(dim); //the sorted entries of v_gen;
			val_v_gen[0]=v_gen(0,0)=(type_int)1;
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
			type_int i_dim,i_size,i_unit;
			for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
			for(i_unit=1;i_unit<dim_max;++i_unit){
				//first, check if unit[i_unit] can generate a feasible generating vector.
				for(i_dim=1;i_dim<dim;++i_dim) val_v_gen[i_dim]=(val_v_gen[i_dim-1]*unit[i_unit])%size;
				std::sort(val_v_gen.begin(),val_v_gen.end());
				for(i_dim=1;i_dim<dim;++i_dim){
					if(val_v_gen[i_dim]==val_v_gen[i_dim-1]) i_dim=dim;
				}
				if(i_dim>dim) continue;
				for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
				for(i_size=0;i_size<size;++i_size){
					for(i_dim=0;i_dim<dim;++i_dim){
						des[i_size](i_dim,0)=crd_can[(v_gen(i_dim,0)*i_size)%size];
					}
				}
				v_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
				if(v_crit<val_crit || tag_first_des){
					tag_first_des=false;
					val_crit=v_crit;
					vec_gen=v_gen;
				}
			}
			stat= tag_first_des ? (Type_Stat)2 : (Type_Stat)0;
			return stat;
		}
		//Cal: find the best GLP over all power vector_gen:
		//(without checking whether the components of the generating vector are different from each other):
		template<typename CritType>
		Type_Stat Cal_best_power_nocheck(type_int const & num_dim,CritType & fun_crit_can){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			dim=num_dim;
			if(dim<=(type_int)0){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			Type_Bool tag_first_des=true;
			type_real v_crit; //inner value of criterion;
			Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
			type_pt_fac v_gen(dim,1); //inner vec_gen;
			v_gen(0,0)=(type_int)1;
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
			type_int i_dim,i_size,i_unit;
			for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
			for(i_unit=(dim_max>(type_int)1 ? (type_int)1 : (type_int)0);i_unit<dim_max;++i_unit){
				for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
				for(i_size=0;i_size<size;++i_size){
					for(i_dim=0;i_dim<dim;++i_dim){
						des[i_size](i_dim,0)=crd_can[(v_gen(i_dim,0)*i_size)%size];
					}
				}
				v_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
				if(v_crit<val_crit || tag_first_des){
					tag_first_des=false;
					val_crit=v_crit;
					vec_gen=v_gen;
				}
			}
			stat= tag_first_des ? (Type_Stat)2 : (Type_Stat)0;
			return stat;
		}
		//Cal: find the best leave-one-out GLP over all power vector_gen:
		//(without checking whether the components of the generating vector are different from each other):
		template<typename CritType>
		Type_Stat Cal_best_power_LeaveOneOut_nocheck(type_int const & num_dim,CritType & fun_crit_can){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			dim=num_dim;
			if(dim<=(type_int)0){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			type_int size_loo=size-(type_int)1; //size of the leave-one-out (``loo'') GLP;
			Type_ArrayTemp<type_real> crd_can_loo(size_loo);
			Type_Bool tag_first_des=true;
			type_real v_crit; //inner value of criterion;
			Type_ArrayTemp<type_pt_can> des(size_loo); //inner canonical design;
			type_pt_fac v_gen(dim,1); //inner vec_gen;
			v_gen(0,0)=(type_int)1;
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size_loo,size_loo);
			type_int i_dim,i_size,i_unit;
			v_crit=(type_real)size_loo;
			for(i_size=0;i_size<size_loo;++i_size){
				des[i_size].resize(dim,1);
				crd_can_loo[i_size]=((type_real)i_size+(type_real)0.5)/v_crit;
			}
			for(i_unit=1;i_unit<dim_max;++i_unit){
				for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=(v_gen(i_dim-1,0)*unit[i_unit])%size;
				for(i_size=0;i_size<size_loo;++i_size){
					for(i_dim=0;i_dim<dim;++i_dim){
						des[i_size](i_dim,0)=
							crd_can_loo[(v_gen(i_dim,0)*(i_size+(type_int)1))%size-(type_int)1];
					}
				}
				v_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
				if(v_crit<val_crit || tag_first_des){
					tag_first_des=false;
					val_crit=v_crit;
					vec_gen=v_gen;
				}
			}
			stat= tag_first_des ? (Type_Stat)2 : (Type_Stat)0;
			return stat;
		}
		//Cal: find a rough GLP fast:
		template<typename CritType>
		Type_Stat Cal_roughfast(type_int const & num_dim,CritType & fun_crit_can){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			dim=num_dim;
			if(dim<=(type_int)0 || dim>dim_max){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			type_int i_dim,i_size;
			vec_gen.resize(dim,1);
			Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
			for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
			type_real tr0=(type_real)dim_max/(type_real)dim;
			for(i_dim=0;i_dim<dim;++i_dim){
				vec_gen(i_dim,0)=
					unit[std::min((type_int)(tr0*(type_real)i_dim+(type_real)0.5),dim_max-(type_int)1)];
			}
			for(i_size=0;i_size<size;++i_size){
				for(i_dim=0;i_dim<dim;++i_dim){
					des[i_size](i_dim,0)=crd_can[(vec_gen(i_dim,0)*i_size)%size];
				}
			}
			val_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
			stat=(Type_Stat)0;
			return stat;
		}
		//Cal: find the best GLP over all vector_gen which is compatible with the cell of another GLP:
		//cellbase_fac: level_another times of the basis matrix of the cell of
		//another GLP with level_another levels:
		template<typename CritType,typename MatRealType>
		Type_Stat Cal_best_cellbasecompat
		(type_int const & num_dim,CritType & fun_crit_can,Type_MatTemp<MatRealType> const & cellbase_fac){
			LZ_DEF_func_check_traits((std::is_same<
				typename std::result_of<typename std::decay<CritType>::type(
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					typename Type_ArrayTemp<type_pt_can>::const_iterator,
					Iterator_constval<type_real>,
					Iterator_constval<type_real>
				)>::type,
				type_real>::value));
			LZ_DEF_func_check_traits((std::is_arithmetic<MatRealType>::value));
			LZ_DEF_func_check_traits((std::is_convertible<type_int,MatRealType>::value));
			LZ_DEF_func_check_traits((std::is_convertible<MatRealType,type_int>::value));
			dim=num_dim;
			if(dim<=(type_int)0 || dim>dim_max ||
				cellbase_fac.rows()!=dim || cellbase_fac.cols()!=dim){
				stat=(Type_Stat)1;
				return stat;
			}
			//the following begins to search for the best design.
			typedef MatRealType type_real;
			typedef typename std::make_unsigned<type_int>::type type_uint;
			typedef typename std::make_signed<type_int>::type type_sint;
			vec_gen.resize(dim,1);
			Type_Bool tag_first_des=true;
			type_real v_crit; //inner value of criterion;
			Type_ArrayTemp<type_pt_can> des(size); //inner canonical design;
			type_pt_fac v_gen(dim,1); //inner vec_gen;
			v_gen(0,0)=(type_int)1;
			math::Comb_enum_comb<type_uint> iter_comb(dim_max,dim);
			Iterator_constval<type_real> iter_wt((type_real)1/(type_real)size,size);
			type_int i_dim,j_dim,i_size;
			type_real v_gen_new0; //component of the transformed v_gen by cellbase_fac;
			for(i_size=0;i_size<size;++i_size) des[i_size].resize(dim,1);
			while((*iter_comb)[0]==(type_uint)0){
				for(i_dim=1;i_dim<dim;++i_dim) v_gen(i_dim,0)=unit[(*iter_comb)[i_dim]];
				++iter_comb;
				//check if v_gen is compatible with cellbase_fac:
				for(i_dim=0;i_dim<dim;++i_dim){
					v_gen_new0=(type_real)0;
					for(j_dim=0;j_dim<dim;++j_dim){
						v_gen_new0+=(cellbase_fac(i_dim,j_dim)*(type_real)(v_gen(j_dim,0)));
					}
					i_size=math::mod((type_sint)(std::round(v_gen_new0)),size);
					if(!(std::binary_search(unit.begin(),unit.end(),i_size))) break; //not compatible;
				}
				if(i_dim<dim) continue; //not compatible;
				for(i_size=0;i_size<size;++i_size){
					for(i_dim=0;i_dim<dim;++i_dim){
						des[i_size](i_dim,0)=crd_can[(v_gen(i_dim,0)*i_size)%size];
					}
				}
				v_crit=fun_crit_can(des.begin(),des.end(),iter_wt,iter_wt.end());
				if(v_crit<val_crit || tag_first_des){
					tag_first_des=false;
					val_crit=v_crit;
					vec_gen=v_gen;
				}
			}
			stat= tag_first_des ? (Type_Stat)2 : (Type_Stat)0;
			return stat;
		}
		
		//Get:
		Type_Stat Get_stat() const {
			return stat;
		}
		type_pt_fac Get_vec_gen() const {
			return vec_gen;
		}
		type_real Get_val_crit() const {
			return stat==(Type_Stat)0 ? val_crit : (type_real)(-1);
		}
		type_int Get_size() const {
			return size;
		}
		type_int Get_dim() const {
			return dim;
		}
		Type_ArrayTemp<type_int> Get_unit() const {
			return unit;
		}
		Type_ArrayTemp<type_real> Get_design_unary() const {
			return crd_can;
		}
		
		//Set:
		Type_Stat Set_reset(type_int const & num_size){
			if(num_size<=(type_int)0) return (Type_Stat)1;
			if(num_size==size) return (Type_Stat)0;
			size=num_size;
			//the following calculates dim_max and unit.
			dim_max=(type_int)1;
			Type_ArrayTemp<type_int> unit_init(size-1);
			unit_init[0]=(type_int)1;
			type_int i;
			for(i=2;i<size;++i){
				if(math::Is_coprime(size,i)) unit_init[dim_max++]=i;
			}
			unit.assign(unit_init.begin(),unit_init.begin()+dim_max);
			//the following calculates crd_can.
			crd_can.resize(size);
			type_real tr0=size,tr1=(type_real)1/((type_real)(2*size));
			for(i=0;i<size;++i) crd_can[i]=(type_real)i/tr0+tr1;
			return (Type_Stat)0;
		}
	protected:
		//val:
		type_int size,dim; //size and dim of this design;
		type_int dim_max; //the max dim under this size;
		Type_Stat stat; //the state: 0: valid parameter; 1: invalid parameter;
		Type_ArrayTemp<type_int> unit; //{x : 1<x<size and gcd(x,size)=1};
		Type_ArrayTemp<type_real> crd_can; //coordinates for the canonical design;
		
		//val: results:
		type_pt_fac vec_gen; //the generation vector;
		type_real val_crit; //the value of the criterion;
}; //class cls_Get_design_GoodLatticePoint;

//fun: Get_design_canon_Augment_slice:
template<typename UIntType=Type_UInt,
	typename InIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::const_iterator,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_Augment_slice
(Type_ArrayTemp<UIntType> const & design_size,Type_ArrayTemp<InIterType> const & design_begin,
OutIterType result_begin){
	typedef UIntType type_uint;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	LZ_DEF_func_check_traits((std::is_unsigned_integral<type_uint>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::remove_cv<typename InIterType::value_type>::type,type_matrix>::value));
	//value check{
	type_uint n_slice=design_size.size(),i_slice;
	if(n_slice==(type_uint)0 || n_slice!=design_begin.size()) return result_begin;
	type_uint dim=(*(design_begin[0])).rows(),i_dim;
	if(dim<=(type_uint)0) return result_begin;
	for(i_slice=0;i_slice<n_slice;++i_slice){
		if(design_size[i_slice]<=(type_uint)0 ||
			(*(design_begin[i_slice])).rows()!=dim || (*(design_begin[i_slice])).cols()<(type_uint)1){
			return result_begin;
		}
	}
	//}
	Type_ArrayTemp<InIterType> iter_des=design_begin;
	Type_ArrayTemp<type_uint> i_size(n_slice,(type_uint)0);
	while(true){
		//augmentation{
		(*result_begin).resize(dim,1);
		for(i_dim=0;i_dim<dim;++i_dim){
			(*result_begin)(i_dim,0)=(type_real)0;
			for(i_slice=0;i_slice<n_slice;++i_slice){
				(*result_begin)(i_dim,0)+=(*(iter_des[i_slice]))(i_dim,0);
			}
			//mapping from [0,+infinity) to [0,1):
			while((*result_begin)(i_dim,0)>=(type_real)1) (*result_begin)(i_dim,0)-=(type_real)1;
		}
		++result_begin;
		//}
		//getting the next point (similar to Liuze::math::Comb_enum_cart::operator ++){
		i_slice=n_slice-(type_uint)1;
		while(i_size[i_slice]+(type_uint)1==design_size[i_slice]){
			if(i_slice==(type_uint)0) return result_begin;
			--i_slice;
		}
		++(i_size[i_slice]);
		++(iter_des[i_slice]);
		while((++i_slice)<n_slice){
			i_size[i_slice]=(type_uint)0;
			iter_des[i_slice]=design_begin[i_slice];
		}
		//}
	}
} //(fun: Get_design_canon_Augment_slice)

//fun: infun_Get_design_canon_Augment_GLPcell_matrixbase:
//used in function Get_design_canon_Augment_GLPcell,
//for getting the basis vectors of the cell of the GLP set:
template<typename IntType=Type_UInt,typename RealType=Type_Real>
Type_MatTemp<RealType> infun_Get_design_canon_Augment_GLPcell_matrixbase
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen){
	typedef IntType type_int;
	typedef RealType type_real;
	typedef Type_MatTemp<type_real> type_matrix;
	typedef Type_ArrayTemp<type_matrix> type_des;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	//value check{
	type_int const dim=vector_gen.rows();
	type_matrix mat_base(0,0);
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return mat_base;
	//}
	mat_base=type_matrix::Zero(dim,dim);
	type_matrix point(dim,1);
	type_matrix & mat_base_part=point;
	Type_ArrayTemp<type_real> crd(num_level); //all the 1-dim coordinates;
	Type_ArrayTemp<type_int> point_id(dim,(type_int)0);
	type_des GLP(num_level),GLP_shift(num_level); //GLP_shift is 0-centred and shifted to close to 0;
	type_real tr0;
	Type_ArrayTemp<type_real> dist(num_level,(type_real)0);
	Type_ArrayTemp<Type_Size> dist_id(num_level);
	Eigen::ColPivHouseholderQR<type_matrix> QRsolver;
	type_int i_dim,i_level,i_size;
	for(i_level=0;i_level<num_level;++i_level){
		crd[i_level]=((type_real)i_level+(type_real)(0.5))/((type_real)num_level);
	}
	auto fun_Get_point=[&i_dim,&dim,&point,&crd,&point_id]()->type_matrix&{
			for(i_dim=(type_int)0;i_dim<dim;++i_dim){
				point(i_dim,0)=crd[point_id[i_dim]];
			}
			return point;
		};
	for(i_level=(type_int)0;i_level<num_level;++i_level){
		GLP[i_level]=fun_Get_point();
		for(i_dim=0;i_dim<dim;++i_dim){
			point_id[i_dim]=(point_id[i_dim]+vector_gen(i_dim,0))%num_level;
		}
		dist_id[i_level]=(Type_Size)i_level;
		GLP_shift[i_level]=GLP[i_level]-GLP[0];
		for(i_dim=0;i_dim<dim;++i_dim){
			tr0=GLP_shift[i_level](i_dim,0);
			if(tr0>(type_real)0.5){
				GLP_shift[i_level](i_dim,0)-=(type_real)1;
				tr0=(type_real)1-tr0;
			}
			dist[i_level]+=(tr0*tr0);
		}
	}
	auto fun_compare_dist_id=[&dist](Type_Size const & id_0,Type_Size const & id_1)->bool{
			return dist[id_0]<dist[id_1];
		};
	std::sort(dist_id.begin(),dist_id.end(),fun_compare_dist_id);
	//the precision:
	tr0=std::min((type_real)LZ_DEF_const_default_precision,(type_real)0.1/(type_real)num_level);
	mat_base.col(0)=GLP_shift[dist_id[1]];
	i_dim=1;
	i_size=1;
	i_level=0;
	GTS_findbase:
	while(i_dim<dim){
		mat_base_part=mat_base.block(0,0,dim,i_dim);
		QRsolver.compute(mat_base_part);
		while((++i_size)<num_level &&
			(mat_base_part*type_matrix(QRsolver.solve(GLP_shift[dist_id[i_size]]))).isApprox
			(GLP_shift[dist_id[i_size]],tr0)){}
		if(i_size<num_level){
			mat_base.col(i_dim)=GLP_shift[dist_id[i_size]];
			++i_dim;
		} else break;
	}
	if(i_size==num_level && i_level==(type_int)0){
		for(i_level=(type_int)0;i_level<num_level;++i_level){
			dist_id[i_level]=(Type_Size)i_level;
			GLP_shift[i_level]=GLP[i_level]-GLP[0];
			for(i_size=0;i_size<dim;++i_size){
				tr0=GLP_shift[i_level](i_size,0);
				dist[i_level]+=(tr0*tr0);
			}
		}
		std::sort(dist_id.begin(),dist_id.end(),fun_compare_dist_id);
		i_size=0;
		i_level=1;
		goto GTS_findbase;
	}
	return mat_base;
} //(fun: infun_Get_design_canon_Augment_GLPcell_matrixbase)
//fun: Get_design_canon_Augment_GLPcell:
template<typename IntType=Type_UInt,
	typename InIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::const_iterator,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_Augment_GLPcell
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,
InIterType const & design_begin,InIterType const & design_end,OutIterType result_begin){
	typedef IntType type_int;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	typedef Type_ArrayTemp<type_matrix> type_des;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	LZ_DEF_func_check_traits((std::is_same<
		typename std::remove_cv<typename InIterType::value_type>::type,type_matrix>::value));
	//value check{
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0) return result_begin;
	if(design_begin==design_end) return result_begin;
	//}
	type_int i_dim,i_level,i_size;
	type_matrix point(dim,1);
	type_des GLP(num_level),GLP_shift(0);
	Get_design_canon_GoodLatticePoint(num_level,vector_gen,GLP.begin());
	type_matrix mat_base=
		infun_Get_design_canon_Augment_GLPcell_matrixbase<type_int,type_real>(num_level,vector_gen);
	InIterType iter_des=design_begin;
	for(iter_des=design_begin;iter_des!=design_end;++iter_des){
		GLP_shift.push_back(mat_base*(*iter_des));
	}
	for(i_level=0;i_level<num_level;++i_level){
		for(i_size=0;i_size<GLP_shift.size();++i_size){
			for(i_dim=0;i_dim<dim;++i_dim){
				point(i_dim,0)=GLP[i_level](i_dim,0)+GLP_shift[i_size](i_dim,0);
				if(point(i_dim,0)<(type_real)0) point(i_dim,0)+=(type_real)1;
				if(point(i_dim,0)>=(type_real)1) point(i_dim,0)-=(type_real)1;
			}
			*result_begin=point;
			++result_begin;
		}
	}
	return result_begin;
} //(fun: Get_design_canon_Augment_GLPcell)
//fun: Get_design_canon_Augment_GLPcell_test:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_Augment_GLPcell_test
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,
IntType const & num_level_fill,Type_MatTemp<IntType> const & vector_gen_fill,
OutIterType result_begin){
	typedef IntType type_int;
	typedef typename std::make_signed<type_int>::type type_sint;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	typedef Type_ArrayTemp<type_matrix> type_des;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	//value check{
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0 ||
		num_level_fill<=0 || dim!=vector_gen_fill.rows() || vector_gen_fill.cols()<=0){
		return result_begin;
	}
	//}
	type_int i_dim,j_dim;
	type_matrix mat_base=(type_real)num_level*
		infun_Get_design_canon_Augment_GLPcell_matrixbase<type_int,type_real>(num_level,vector_gen);
	type_matrix mat_base_fill=(type_real)num_level_fill*
		infun_Get_design_canon_Augment_GLPcell_matrixbase<type_int,type_real>(num_level_fill,vector_gen_fill);
	type_int num_level_inv,num_level_fill_inv;
	for(num_level_inv=1;num_level_inv<num_level_fill;++num_level_inv){
		if((num_level_inv*num_level)%num_level_fill==(type_int)1) break;
	}
	for(num_level_fill_inv=1;num_level_fill_inv<num_level;++num_level_fill_inv){
		if((num_level_fill_inv*num_level_fill)%num_level==(type_int)1) break;
	}
	if(num_level_inv==num_level_fill) return result_begin;
	Type_MatTemp<type_int> vec_gen(dim,1);
	type_sint ti0,ti1;
	type_int num_level_total=num_level*num_level_fill;
	for(i_dim=0;i_dim<dim;++i_dim){
		ti1=ti0=(type_sint)0;
		for(j_dim=0;j_dim<dim;++j_dim){
			ti0+=((type_sint)(std::round(mat_base(i_dim,j_dim)))*(type_sint)(vector_gen_fill(j_dim,0)));
			ti1+=((type_sint)(std::round(mat_base_fill(i_dim,j_dim)))*(type_sint)(vector_gen(j_dim,0)));
		}
		ti0*=(type_sint)(num_level_inv*num_level);
		ti1*=(type_sint)(num_level_fill_inv*num_level_fill);
		vec_gen(i_dim,0)=math::mod(ti0+ti1,num_level_total);
	}
	return Get_design_canon_GoodLatticePoint(num_level_total,vec_gen,result_begin);
} //(fun: Get_design_canon_Augment_GLPcell_test)
//fun: Get_design_canon_Augment_GLPcell_vecgen:
template<typename IntType=Type_UInt,
	typename OutIterType=typename Type_ArrayTemp<Type_MatTemp<Type_Real> >::iterator>
OutIterType Get_design_canon_Augment_GLPcell_vecgen
(IntType const & num_level,Type_MatTemp<IntType> const & vector_gen,
IntType const & num_level_fill,Type_MatTemp<IntType> const & vector_gen_fill,
OutIterType result_begin){
	typedef IntType type_int;
	typedef typename std::make_signed<type_int>::type type_sint;
	typedef typename OutIterType::value_type type_matrix;
	typedef typename type_matrix::value_type type_real;
	typedef Type_ArrayTemp<type_matrix> type_des;
	LZ_DEF_func_check_traits((std::is_integral<type_int>::value));
	LZ_DEF_func_check_traits((std::is_floating_point<type_real>::value));
	LZ_DEF_func_check_traits((std::is_same<Type_MatTemp<type_real>,type_matrix>::value));
	//value check{
	type_int const dim=vector_gen.rows();
	if(num_level<=0 || dim<=0 || vector_gen.cols()<=0 ||
		num_level_fill<=0 || dim!=vector_gen_fill.rows() || vector_gen_fill.cols()<=0){
		return result_begin;
	}
	//}
	type_int i_dim,j_dim;
	type_matrix mat_base=(type_real)num_level*
		infun_Get_design_canon_Augment_GLPcell_matrixbase<type_int,type_real>(num_level,vector_gen);
	type_int num_level_inv;
	for(num_level_inv=1;num_level_inv<num_level_fill;++num_level_inv){
		if((num_level_inv*num_level)%num_level_fill==(type_int)1) break;
	}
	if(num_level_inv==num_level_fill) return result_begin;
	Type_MatTemp<type_int> vec_gen(dim,1);
	type_sint ti0;
	type_int num_level_total=num_level*num_level_fill;
	for(i_dim=0;i_dim<dim;++i_dim){
		ti0=(type_sint)0;
		for(j_dim=0;j_dim<dim;++j_dim){
			ti0+=((type_sint)(std::round(mat_base(i_dim,j_dim)))*(type_sint)(vector_gen_fill(j_dim,0)));
		}
		ti0*=(type_sint)(num_level_inv*num_level);
		ti0+=(type_sint)(vector_gen(i_dim,0)*num_level_fill);
		vec_gen(i_dim,0)=math::mod(ti0,num_level_total);
	}
	return Get_design_canon_GoodLatticePoint(num_level_total,vec_gen,result_begin);
} //(fun: Get_design_canon_Augment_GLPcell_vecgen)

} //namespace stats;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_temp210619