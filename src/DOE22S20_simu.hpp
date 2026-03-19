//DOE22S20=bib_pp_394

#ifndef LZ_DEF_hpp_DOE22S20_simu
#define LZ_DEF_hpp_DOE22S20_simu

#include<random>
#include<functional>
#include<LZ_H/stats/DOE_crit_uniformity.hpp>
#include"OofA_des_strat.hpp"

template<typename RealType=Liuze::Type_Real,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
class DOE22S20_data_generator{
	public:
		//type:
		typedef DOE22S20_data_generator<RealType> type_this;
		typedef RealType type_real;
		typedef Liuze::Type_ArrayTemp<type_real> type_arr_real;
		typedef Liuze::Type_MatTemp<type_real> type_mat_real;
		typedef Liuze::Type_Size type_size;
		
		//static value:
		static std::default_random_engine rand;
		
		//constructor:
		template<typename RandEngType=std::default_random_engine>
		DOE22S20_data_generator
		(type_size n_job=0,type_real pow_cost=1,RandEngType * ptr_rand=NULL):
			m_n_job(n_job>=type_size(0) ? n_job : type_size(0)),
			m_powcost(pow_cost>type_real(0) ? pow_cost : type_real(1)){
			m_weight.resize(m_n_job);
			m_time.resize(m_n_job);
			std::chi_squared_distribution<type_real> distr_chisq(type_real(1));
			type_size i_job;
			for(i_job=0;i_job<m_n_job;++i_job){
				m_weight[i_job]= ptr_rand==NULL ? distr_chisq(rand) : distr_chisq(*ptr_rand);
				m_time[i_job]= ptr_rand==NULL ? distr_chisq(rand) : distr_chisq(*ptr_rand);
			}
		}
		
		//function:
		template<typename DesMatType>
		type_mat_real response(DesMatType const & mat_OofA) const {
			typedef DesMatType type_mat_des;
			type_size n_run=mat_OofA.rows();
			type_size n_step=mat_OofA.cols();
			type_size i_run,i_step;
			type_real resp0;
			type_mat_real mat_resp(n_run,1);
			for(i_run=0;i_run<n_run;++i_run){
				mat_resp(i_run,0)=resp0=type_real(0);
				for(i_step=0;i_step<n_step;++i_step){
					resp0+=m_time[mat_OofA(i_run,i_step)];
					mat_resp(i_run,0)+=
						m_weight[mat_OofA(i_run,i_step)]*type_real(pow(resp0,m_powcost));
				}
				if(m_powcost>type_real(1)){
					mat_resp(i_run,0)=pow(mat_resp(i_run,0),type_real(1)/m_powcost);
				}
			}
			return mat_resp;
		}
		
		type_arr_real const & weight(void) const {
			return m_weight;
		}
		type_arr_real const & time(void) const {
			return m_time;
		}
	protected:
		//value:
		type_size m_n_job;
		type_real m_powcost;
		type_arr_real m_weight;
		type_arr_real m_time;
}; //class DOE22S20_data_generator;
template<typename RealType,typename TraitsCheck>
std::default_random_engine
DOE22S20_data_generator<RealType,TraitsCheck>::rand;

template<typename IntType=Liuze::Type_Int,typename RealType=Liuze::Type_Real,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value &&
	std::is_integral<IntType>::value>::type>
class DesCritVal_data_generator{
	public:
		//type:
		typedef DesCritVal_data_generator<IntType,RealType> type_this;
		typedef RealType type_real;
		typedef IntType type_int;
		typedef Liuze::Type_MatTemp<type_real> type_mat_real;
		typedef Liuze::Type_MatTemp<type_int> type_mat_int;
		typedef std::function<type_real(type_mat_int const &)> type_crit;
		typedef Liuze::Type_Stat type_stat;
		
		//static value:
		static type_stat const stat_leftact=0;
		static type_stat const stat_rightact=1;
		
		//constructor:
		DesCritVal_data_generator
		(type_mat_int const & matdes,type_crit const & crit,
		type_stat act_type=type_this::stat_leftact):
		m_matdes(matdes),m_crit(crit),
		m_stat_act
		(act_type==type_this::stat_rightact ? type_this::stat_rightact : type_this::stat_leftact)
		{}
		
		//function:
		type_mat_real response(type_mat_int const & mat_PofC) const {
			if(m_matdes.cols()==0 || m_matdes.rows()<2 ||
				m_matdes.cols()!=mat_PofC.cols() || mat_PofC.rows()==0){
				return type_mat_real(0,0);
			}
			type_int n_run=mat_PofC.rows();
			type_int i_run;
			type_mat_real mat_resp(n_run,1);
			for(i_run=0;i_run<n_run;++i_run){
				if(m_stat_act==type_this::stat_leftact){
					mat_resp(i_run,0)=m_crit(
						act_perm_OofADesign_left(type_mat_int(mat_PofC.row(i_run)),m_matdes));
				} else if(m_stat_act==type_this::stat_rightact){
					mat_resp(i_run,0)=m_crit(
						act_perm_OofADesign_left(m_matdes,type_mat_int(mat_PofC.row(i_run))));
				} else {
					mat_resp(i_run,0)=type_real(0);
				}
			}
			return mat_resp;
		}
	protected:
		//value:
		type_mat_int m_matdes;
		type_crit m_crit;
		type_stat m_stat_act;
}; //class DesCritVal_data_generator;

template<typename RealType=Liuze::Type_Real,typename DesMatType,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
Liuze::Type_MatTemp<RealType> Get_mat_model_OofA_poly1(DesMatType const & mat_PofC){
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<RealType> type_mat_real;
	typedef Liuze::Type_Size type_size;
	type_size n_run=mat_PofC.rows();
	type_size n_comp=mat_PofC.cols();
	switch(n_comp){
		case type_size(0):
			return type_mat_real(n_run,0);
		case type_size(1):
			return type_mat_real::Ones(n_run,1);
	}
	type_real coef1=sqrt(type_real(12)/(type_real(n_comp)*type_real(n_comp)-type_real(1)));
	type_real pos_mid=(type_real(n_comp)-type_real(1))/type_real(2);
	type_size i_run,i_comp;
	type_mat_real mat_model(n_run,n_comp+1);
	for(i_run=0;i_run<n_run;++i_run){
		mat_model(i_run,0)=type_real(1);
		for(i_comp=0;i_comp<n_comp;++i_comp){
			mat_model(i_run,i_comp+type_size(1))=coef1*(type_real(mat_PofC(i_run,i_comp))-pos_mid);
		}
	}
	return mat_model;
} //fun: Get_mat_model_OofA_poly1;

template<typename IntType=Liuze::Type_Size,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Get_rank_model_OofA_poly1(IntType const & n_comp){
	typedef IntType type_int;
	return n_comp>=type_int(0) ? n_comp : type_int(0);
} //fun: Get_rank_model_OofA_poly1;

template<typename RealType=Liuze::Type_Real,typename DesMatType,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
Liuze::Type_MatTemp<RealType> Get_mat_model_OofA_poly2(DesMatType const & mat_PofC){
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<RealType> type_mat_real;
	typedef Liuze::Type_Size type_size;
	type_size n_run=mat_PofC.rows();
	type_size n_comp=mat_PofC.cols();
	switch(n_comp){
		case type_size(0):
			return type_mat_real(n_run,0);
		case type_size(1):
			return type_mat_real::Ones(n_run,1);
		case type_size(2):
			return Get_mat_model_OofA_poly1(mat_PofC);
	}
	type_real coef2=type_real(n_comp);
	coef2*=coef2;
	type_real coef1=sqrt(type_real(12)/(coef2-type_real(1)));
	coef2=sqrt((type_real(5)*(coef2-type_real(1)))/(type_real(4)*(coef2-type_real(4))));
	type_real pos_mid=(type_real(n_comp)-type_real(1))/type_real(2);
	type_size i_run,i_col;
	type_mat_real mat_model(n_run,type_size(1)+n_comp+n_comp);
	for(i_run=0;i_run<n_run;++i_run){
		mat_model(i_run,0)=type_real(1);
		for(i_col=0;i_col<n_comp;++i_col){
			mat_model(i_run,i_col+type_size(1))=coef1*(type_real(mat_PofC(i_run,i_col))-pos_mid);
			mat_model(i_run,n_comp+type_size(1)+i_col)=coef2
				*(mat_model(i_run,i_col+type_size(1))*mat_model(i_run,i_col+type_size(1))
				-type_real(1));
		}
	}
	return mat_model;
} //fun: Get_mat_model_OofA_poly2;

template<typename IntType=Liuze::Type_Size,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Get_rank_model_OofA_poly2(IntType const & n_comp){
	typedef IntType type_int;
	return n_comp>=type_int(3) ? n_comp+n_comp-type_int(1) : Get_rank_model_OofA_poly1(n_comp);
} //fun: Get_rank_model_OofA_poly2;

template<typename RealType=Liuze::Type_Real,typename DesMatType,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
Liuze::Type_MatTemp<RealType> Get_mat_model_OofA_poly2var2(DesMatType const & mat_PofC){
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<RealType> type_mat_real;
	typedef Liuze::Type_Size type_size;
	type_size n_run=mat_PofC.rows();
	type_size n_comp=mat_PofC.cols();
	switch(n_comp){
		case type_size(0):
			return type_mat_real(n_run,0);
		case type_size(1):
			return type_mat_real::Ones(n_run,1);
		case type_size(2):
			return Get_mat_model_OofA_poly1(mat_PofC);
	}
	type_real coef2=type_real(n_comp);
	coef2*=coef2;
	type_real coef1=sqrt(type_real(12)/(coef2-type_real(1)));
	coef2=sqrt((type_real(5)*(coef2-type_real(1)))/(type_real(4)*(coef2-type_real(4))));
	type_real pos_mid=(type_real(n_comp)-type_real(1))/type_real(2);
	type_size i_run,i_col;
	type_mat_real mat_model(n_run,
		type_size(1)+n_comp+n_comp+Liuze::math::Comb_num::comb<type_size>(n_comp,type_size(2)));
	for(i_run=0;i_run<n_run;++i_run){
		mat_model(i_run,0)=type_real(1);
		for(i_col=0;i_col<n_comp;++i_col){
			mat_model(i_run,i_col+type_size(1))=coef1*(type_real(mat_PofC(i_run,i_col))-pos_mid);
			mat_model(i_run,n_comp+type_size(1)+i_col)=coef2
				*(mat_model(i_run,i_col+type_size(1))*mat_model(i_run,i_col+type_size(1))
				-type_real(1));
		}
	}
	i_col=n_comp+n_comp+type_size(1);
	Liuze::math::Comb_enum_comb<type_size> iter_fi2(n_comp,type_size(2));
	while(iter_fi2){
		for(i_run=0;i_run<n_run;++i_run){
			mat_model(i_run,i_col)=mat_model(i_run,(*iter_fi2)[0]+type_size(1))
				*mat_model(i_run,(*iter_fi2)[1]+type_size(1));
		}
		++iter_fi2;
		++i_col;
	}
	return mat_model;
} //fun: Get_mat_model_OofA_poly2var2;

template<typename IntType=Liuze::Type_Size,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Get_rank_model_OofA_poly2var2(IntType const & n_comp){
	typedef IntType type_int;
	return n_comp>=type_int(3) ? (n_comp-type_int(1))*(n_comp+type_int(2))/type_int(2) :
		Get_rank_model_OofA_poly1(n_comp);
} //fun: Get_rank_model_OofA_poly2var2;

template<typename RealType=Liuze::Type_Real,typename DesMatType,typename StratType,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
Liuze::Type_MatTemp<RealType> Get_mat_model_OofA_SCP1
(DesMatType const & mat_PofC,StratType const & vec_strat){
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<RealType> type_mat_real;
	typedef Liuze::Type_Size type_size;
	type_size n_run=mat_PofC.rows();
	type_size n_comp=mat_PofC.cols();
	type_size n_strat=0;
	type_size i_comp;
	for(i_comp=0;i_comp<n_comp;++i_comp){
		if(vec_strat[i_comp]>n_strat) n_strat=vec_strat[i_comp];
		if(vec_strat[i_comp]<0) return type_mat_real(n_run,0);
	}
	if(vec_strat.size()!=n_comp || n_comp==type_size(0) || n_strat>=n_comp){
		return type_mat_real(n_run,0);
	}
	++n_strat;
	if(n_comp==type_size(1) || n_strat==type_size(1)) return type_mat_real::Ones(n_run,1);
	type_size i_run,i_strat,i_col;
	type_mat_real mat_model(n_run,type_size(1)+n_comp*n_strat);
	for(i_run=0;i_run<n_run;++i_run){
		mat_model(i_run,0)=type_real(1);
		i_col=1;
		for(i_comp=0;i_comp<n_comp;++i_comp){
			for(i_strat=0;i_strat<n_strat;++i_strat){
				mat_model(i_run,i_col)= type_size(vec_strat[mat_PofC(i_run,i_comp)])==i_strat ?
					type_real(1) : type_real(0);
				++i_col;
			}
		}
	}
	return mat_model;
} //fun: Get_mat_model_OofA_SCP1;

template<typename IntType=Liuze::Type_Size,typename StratType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Get_rank_model_OofA_SCP1(StratType const & vec_strat){
	typedef IntType type_int;
	type_int n_comp=vec_strat.size();
	if(n_comp==type_int(0)) return type_int(0);
	type_int strat_max=0;
	type_int i_comp;
	for(i_comp=0;i_comp<n_comp;++i_comp){
		if(vec_strat[i_comp]>strat_max) strat_max=vec_strat[i_comp];
		if(vec_strat[i_comp]<0) return type_int(0);
	}
	if(strat_max>=n_comp) return type_int(0);
	return type_int(1)+(n_comp-type_int(1))*strat_max;
} //fun: Get_rank_model_OofA_SCP1;

template<typename RealType=Liuze::Type_Real,typename DesMatType,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
Liuze::Type_MatTemp<RealType> Get_mat_model_OofA_PWO(DesMatType const & mat_PofC){
	typedef RealType type_real;
	typedef Liuze::Type_MatTemp<RealType> type_mat_real;
	typedef Liuze::Type_Size type_size;
	type_size n_run=mat_PofC.rows();
	type_size n_comp=mat_PofC.cols();
	if(n_comp==type_size(0)) return type_mat_real(n_run,0);
	if(n_comp==type_size(1)) return type_mat_real::Ones(n_run,1);
	type_mat_real mat_model
		(n_run,Liuze::math::Comb_num::comb<type_size>(n_comp,type_size(2))+type_size(1));
	type_size i_run;
	for(i_run=0;i_run<n_run;++i_run) mat_model(i_run,0)=type_real(1);
	type_size i_col=1;
	Liuze::math::Comb_enum_comb<type_size> iter_PWO(n_comp,type_size(2));
	while(iter_PWO){
		for(i_run=0;i_run<n_run;++i_run){
			mat_model(i_run,i_col)= mat_PofC(i_run,(*iter_PWO)[0])<mat_PofC(i_run,(*iter_PWO)[1]) ?
				type_real(1) : type_real(-1);
		}
		++iter_PWO;
		++i_col;
	}
	return mat_model;
} //fun: Get_mat_model_OofA_PWO;

template<typename IntType=Liuze::Type_Size,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Get_rank_model_OofA_PWO(IntType const & n_comp){
	typedef IntType type_int;
	return n_comp<=type_int(0) ? type_int(0) : type_int(1)+n_comp*(n_comp-type_int(1))/type_int(2);
} //fun: Get_rank_model_OofA_PWO;

template<typename RealType=Liuze::Type_Real,
	typename=typename std::enable_if<std::is_floating_point<RealType>::value>::type>
class DOE22S20_simu_LinReg{
	public:
		//type:
		typedef DOE22S20_simu_LinReg<RealType> type_this;
		typedef RealType type_real;
		typedef Matrix_Inverse_MoorePenrose<type_real> type_MatInv;
		typedef typename type_MatInv::type_size type_size;
		typedef Liuze::Type_MatTemp<type_real> type_mat_real;
		
		//constructor:
		DOE22S20_simu_LinReg(void)=default;
		DOE22S20_simu_LinReg(type_mat_real const & mat_model,type_mat_real const & mat_resp){
			this->compute(mat_model,mat_resp);
		}
		
		//function:
		type_mat_real const & mat_model(void) const {
			return m_mat_model;
		}
		type_mat_real const & mat_resp(void) const {
			return m_mat_resp;
		}
		type_mat_real const & mat_coef(void) const {
			return m_mat_coef;
		}
		type_mat_real const & mat_infor(void) const {
			return m_mat_infor;
		}
		type_mat_real mat_infor_inv(void) const {
			return m_matinv.Get_inverse();
		}
		type_size rank(void) const {
			return m_matinv.Get_dim();
		}
		type_MatInv const & mat_infor_inverter(void) const {
			return m_matinv;
		}
		type_this & compute(type_mat_real const & mat_model,type_mat_real const & mat_resp){
			if(mat_model.rows()!=mat_resp.rows()) return *this;
			m_mat_model=mat_model;
			m_mat_resp=mat_resp;
			m_mat_infor=m_mat_model.transpose()*m_mat_model;
			m_matinv.compute(m_mat_infor);
			m_mat_coef=m_matinv.Get_inverse()*m_mat_model.transpose()*m_mat_resp;
			return *this;
		}
		type_mat_real predict(type_mat_real const & mat_model) const {
			return mat_model.cols()==m_mat_coef.rows() ? mat_model*m_mat_coef : type_mat_real(0,0);
		}
	protected:
		//value:
		type_MatInv m_matinv;
		type_mat_real m_mat_model,m_mat_resp,m_mat_coef,m_mat_infor;
}; //class DOE22S20_simu_LinReg;

#endif //#ifndef LZ_DEF_hpp_DOE22S20_simu