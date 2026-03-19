#ifndef LZ_DEF_hpp_OofA_des_strat
#define LZ_DEF_hpp_OofA_des_strat

#include"OofA.hpp"

/****************************************************************************************************/
//matrix operation{
template<typename EigenDerived,typename FuncType,
	typename=typename std::enable_if<std::is_convertible<
		typename std::remove_reference<typename std::remove_cv<
		typename std::result_of<FuncType(typename Eigen::DenseBase<EigenDerived>::Scalar)>::type
		>::type>::type,
		typename Eigen::DenseBase<EigenDerived>::Scalar
	>::value>::type>
void applyOnMat_entrywise
(FuncType && func,Eigen::DenseBase<EigenDerived> & X){
	typename Eigen::DenseBase<EigenDerived>::Index ir,ic;
	for(ir=0;ir<X.rows();++ir){
		for(ic=0;ic<X.cols();++ic) X(ir,ic)=func(X(ir,ic));
	}
} //fun: applyOnMat_entrywise;
template<typename EigenDerived,typename FuncType,
	typename=typename std::enable_if<std::is_arithmetic<
		typename std::remove_reference<typename std::remove_cv<
		typename std::result_of<FuncType(typename Eigen::DenseBase<EigenDerived>::Scalar)>::type
		>::type>::type
	>::value>::type>
Liuze::Type_MatTemp<
	typename std::remove_reference<typename std::remove_cv<
	typename std::result_of<FuncType(typename Eigen::DenseBase<EigenDerived>::Scalar)>::type
	>::type>::type
>
evalOnMat_entrywise
(FuncType && func,Eigen::DenseBase<EigenDerived> const & X){
	Liuze::Type_MatTemp<
		typename std::remove_reference<typename std::remove_cv<
		typename std::result_of<FuncType(typename Eigen::DenseBase<EigenDerived>::Scalar)>::type
		>::type>::type
	> res(X.rows(),X.cols());
	typename Eigen::DenseBase<EigenDerived>::Index ir,ic;
	for(ir=0;ir<X.rows();++ir){
		for(ic=0;ic<X.cols();++ic) res(ir,ic)=func(X(ir,ic));
	}
	return res;
} //fun: evalOnMat_entrywise;
template<typename Scalar=void,typename EigenDerived_0,typename EigenDerived_1>
Liuze::Type_MatTemp<
	typename std::conditional<std::is_same<Scalar,void>::value,
	typename std::remove_reference<typename std::remove_cv<
	typename Eigen::DenseBase<EigenDerived_0>::Scalar
	>::type>::type,
	typename std::remove_reference<typename std::remove_cv<Scalar>::type>::type
	>::type
>
Matrix_product_Kronecker
(Eigen::DenseBase<EigenDerived_0> const & mat_0,Eigen::MatrixBase<EigenDerived_1> const & mat_1){
	typedef
		typename std::conditional<std::is_same<Scalar,void>::value,
		typename std::remove_reference<typename std::remove_cv<
		typename Eigen::DenseBase<EigenDerived_0>::Scalar
		>::type>::type,
		typename std::remove_reference<typename std::remove_cv<Scalar>::type>::type
		>::type
		type_res_scalar;
	typedef Liuze::Type_MatTemp<type_res_scalar> type_res;
	typedef typename Eigen::DenseBase<EigenDerived_0>::Index type_size;
	type_size ir,ic,ir0,ic0;
	type_res res(mat_0.rows()*mat_1.rows(),mat_0.cols()*mat_1.cols());
	for(ir=0,ir0=0;ir0<mat_0.rows();ir+=mat_1.rows(),++ir0){
		for(ic=0,ic0=0;ic0<mat_0.cols();ic+=mat_1.cols(),++ic0){
			res.block(ir,ic,mat_1.rows(),mat_1.cols())=mat_0(ir0,ic0)*mat_1.array();
		}
	}
	return res;
} //fun: Matrix_product_Kronecker;
//}(matrix operation)
/****************************************************************************************************/
template<typename ValType=Liuze::Type_Int,
	typename=typename std::enable_if<std::is_signed<ValType>::value>::type>
Liuze::Type_MatTemp<ValType> Mat_Hadamard_Sylvester(Liuze::Type_Size const & num_exp){
	typedef Liuze::Type_Size type_size;
	typedef ValType type_val;
	typedef Liuze::Type_MatTemp<type_val> type_mat;
	if(num_exp==type_size(0)) return type_mat::Zero(1,1);
	type_mat doubler=type_mat::Ones(2,2);
	doubler(1,1)=-Liuze::math::one<type_val>();
	if(num_exp==type_size(1)) return doubler;
	type_mat res=doubler;
	type_size i_exp;
	for(i_exp=1;i_exp<num_exp;++i_exp){
		res=Matrix_product_Kronecker(doubler,res);
	}
	return res;
}
/****************************************************************************************************/

//fun: Get_rowBIBD_complement:
//des_rowBIBD: a balanced incomplete block design, in which each row represents a block.
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_rowBIBD_complement
(Liuze::Type_MatTemp<IntType> const & des_rowBIBD,IntType const & num_symbol){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	type_int num_row=des_rowBIBD.rows(),num_col=des_rowBIBD.cols();
	if(num_symbol<=type_int(0) || num_col>=num_symbol) return type_mat_int(0,0);
	if(num_row==type_int(0) || num_col==type_int(0)){
		type_mat_int des_rowBIBD_1(1,num_symbol);
		type_int i_col;
		for(i_col=0;i_col<num_symbol;++i_col) des_rowBIBD_1(0,i_col)=i_col;
		return des_rowBIBD_1;
	}
	//checkings end.
	type_int num_col_1=num_symbol-num_col;
	type_mat_int des_rowBIBD_1(num_row,num_col_1);
	type_int i_row,i_col,i_sym;
	std::vector<bool> flag_sym(num_symbol);
	for(i_row=0;i_row<num_row;++i_row){
		for(i_sym=0;i_sym<num_symbol;++i_sym) flag_sym[i_sym]=true;
		for(i_col=0;i_col<num_col;++i_col) flag_sym[des_rowBIBD(i_row,i_col)]=false;
		for(i_col=0,i_sym=0;i_col<num_col_1;++i_sym){
			if(flag_sym[i_sym]){
				des_rowBIBD_1(i_row,i_col)=i_sym;
				++i_col;
			}
		}
	}
	return des_rowBIBD_1;
} //fun: Get_rowBIBD_complement;

//fun: Get_rowBIBD_Hadamard:
//mat_Hadamard: a Hadamard matrix of order n with entries in {0,1};
//assign: the block size of the returned BIBD
//is n/2-1 if `assign' is 0, and is n/2 if `assign' is 1;
//return: a balanced incomplete block design, in which each row represents a block.
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_rowBIBD_Hadamard
(Liuze::Type_MatTemp<IntType> mat_Hadamard,IntType assign=0){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	if(mat_Hadamard.rows()<=type_int(3) || mat_Hadamard.rows()!=mat_Hadamard.cols()){
		return type_mat_int(0,0);
	}
	if(assign<type_int(0)) assign=type_int(0);
	else if(assign>type_int(1)) assign=type_int(1);
	//checkings end.
	type_int n_row=mat_Hadamard.rows(),n_col=n_row;
	type_int i_row,i_col,i_sym;
	for(i_col=0;i_col<n_col;++i_col){
		if(mat_Hadamard(0,i_col)!=assign){
			for(i_row=0;i_row<n_row;++i_row){
				mat_Hadamard(i_row,i_col)=
					mat_Hadamard(i_row,i_col)==type_int(0) ? type_int(1) : type_int(0);
			}
		}
	}
	for(i_row=0;i_row<n_row;++i_row){
		if(mat_Hadamard(i_row,0)!=assign){
			for(i_col=0;i_col<n_col;++i_col){
				mat_Hadamard(i_row,i_col)=
					mat_Hadamard(i_row,i_col)==type_int(0) ? type_int(1) : type_int(0);
			}
		}
	}
	type_mat_int mat_BIBD
		(n_row-type_int(1),n_col/type_int(2)-(assign==type_int(0) ? type_int(1) : type_int(0)));
	for(i_row=0;i_row<mat_BIBD.rows();++i_row){
		for(i_sym=1,i_col=0;i_col<mat_BIBD.cols();++i_sym){
			if(mat_Hadamard(i_row+type_int(1),i_sym)==type_int(0)){
				mat_BIBD(i_row,i_col)=i_sym-type_int(1);
				++i_col;
			}
		}
	}
	return mat_BIBD;
} //fun: Get_rowBIBD_Hadamard;

//fun: Get_rowBIBD_derived:
//des_rowBIBD: a symmetric BIBD, in which each row represents a block.
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_rowBIBD_derived
(Liuze::Type_MatTemp<IntType> const & des_rowBIBD){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	if(des_rowBIBD.rows()<=type_int(1) || des_rowBIBD.cols()<=type_int(1)){
		return type_mat_int(0,0);
	}
	type_int n_sym=des_rowBIBD.rows(),n_col=des_rowBIBD.cols(),n_col_1=0;
	type_int i_row,i_col,i_sym;
	std::vector<bool> flag_row0(n_sym,false),flag_row(n_sym,false);
	Liuze::Type_ArrayTemp<type_int> map_sym_1(n_sym,n_sym);
	for(i_col=0;i_col<n_col;++i_col) flag_row0[des_rowBIBD(0,i_col)]=true;
	for(i_col=0;i_col<n_col;++i_col){
		if(flag_row0[des_rowBIBD(1,i_col)]) ++n_col_1;
	}
	type_mat_int des_rowBIBD_1(n_sym-type_int(1),n_col_1);
	for(i_row=1;i_row<n_sym;++i_row){
		for(i_col=0;i_col<n_col;++i_col){
			if(flag_row0[des_rowBIBD(i_row,i_col)]) flag_row[des_rowBIBD(i_row,i_col)]=true;
		}
		for(i_col=0,i_sym=0;i_col<n_col_1;++i_sym){
			if(flag_row[i_sym]){
				des_rowBIBD_1(i_row-type_int(1),i_col)=i_sym;
				flag_row[i_sym]=false;
				map_sym_1[i_sym]=type_int(0);
				++i_col;
			}
		}
	}
	for(i_col=0,i_sym=0;i_sym<n_sym;++i_sym){
		if(map_sym_1[i_sym]<n_sym) map_sym_1[i_sym]=i_col++;
	}
	for(i_row=0;i_row<des_rowBIBD_1.rows();++i_row){
		for(i_col=0;i_col<n_col_1;++i_col){
			des_rowBIBD_1(i_row,i_col)=map_sym_1[des_rowBIBD_1(i_row,i_col)];
		}
	}
	return des_rowBIBD_1;
} //fun: Get_rowBIBD_derived;

//fun: Get_GYD_rowregular:
//des_BBD: a balanced block design, in which each row represents a block;
//see Colbourn and Dinitz (2007). Handbook of Combinatorial Designs. Construction VI.65.31.
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_GYD_rowregular
(Liuze::Type_MatTemp<IntType> const & des_BBD,IntType const & num_symbol){
	typedef IntType type_int;
	typedef typename std::make_signed<type_int>::type type_sint;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	typedef Liuze::Type_MatTemp<type_sint> type_mat_sint;
	type_int num_row=des_BBD.rows(),num_col=des_BBD.cols();
	if(num_symbol<=type_int(0) || num_row<=type_int(0) || num_col<=type_int(0)
	|| num_row%num_symbol!=type_int(0)){
		return type_mat_int(0,0);
	}
	//checkings end.
	//initialising the bipartite multigraph{
	type_int i_row,i_col,i_sym,i_col_0,i_col_1;
	type_mat_sint graph=type_mat_sint::Constant(num_symbol,num_row,-type_sint(1));
	type_int count0_sym_col=num_row/num_symbol;
	type_mat_int count_sym_col=type_mat_int::Zero(num_symbol,num_col);
	for(i_row=0;i_row<num_row;++i_row){
		for(i_col=0;i_col<num_col;++i_col){
			graph(des_BBD(i_row,i_col),i_row)=i_col;
			++count_sym_col(des_BBD(i_row,i_col),i_col);
		}
	}
	//}(initialising the bipartite multigraph)
	//adjusting colours of edges{
	type_mat_sint graph_temp;
	Liuze::Type_ArrayTemp<Liuze::Type_ArrayTemp<type_int> > path(0);
	while(true){
		//finding a colour whose number of appearance on some symbol is less than supposed:
		for(i_sym=0;i_sym<num_symbol;++i_sym){
			for(i_col=0;i_col<num_col;++i_col){
				if(count_sym_col(i_sym,i_col)<count0_sym_col) goto GTL_col_minor;
			}
		}
		break;
		GTL_col_minor:
		i_col_0=i_col;
		//finding a colour whose number of appearance on that symbol is greater than supposed:
		for(i_col=0;i_col<num_col;++i_col){
			if(count_sym_col(i_sym,i_col)>count0_sym_col) break;
		}
		i_col_1=i_col;
		//constructing an alternating path{
		path.resize(0);
		graph_temp=graph;
		while(count_sym_col(i_sym,i_col_0)<=count_sym_col(i_sym,i_col_1)){
			for(i_row=0;i_row<num_row;++i_row){
				if(graph_temp(i_sym,i_row)==i_col_1) break;
			}
			path.push_back({i_sym,i_row});
			graph_temp(i_sym,i_row)=-type_sint(1);
			for(i_sym=0;i_sym<num_symbol;++i_sym){
				if(graph_temp(i_sym,i_row)==i_col_0) break;
			}
			path.push_back({i_sym,i_row});
			graph_temp(i_sym,i_row)=-type_sint(1);
		}
		//}(constructing an alternating path)
		//interchanging the colours `i_col_0' and `i_col_1' on `path':
		for(i_col=0;i_col<path.size();){
			graph(path[i_col][0],path[i_col][1])=i_col_0;
			--count_sym_col(path[i_col][0],i_col_1);
			++count_sym_col(path[i_col][0],i_col_0);
			++i_col;
			graph(path[i_col][0],path[i_col][1])=i_col_1;
			--count_sym_col(path[i_col][0],i_col_0);
			++count_sym_col(path[i_col][0],i_col_1);
			++i_col;
		}
	}
	//}(adjusting colours of edges)
	type_mat_int des_GYD(num_row,num_col);
	for(i_sym=0;i_sym<num_symbol;++i_sym){
		for(i_row=0;i_row<num_row;++i_row){
			if(graph(i_sym,i_row)>=type_sint(0)){
				des_GYD(i_row,graph(i_sym,i_row))=i_sym;
			}
		}
	}
	return des_GYD;
} //fun: Get_GYD_rowregular;

/****************************************************************************************************/

template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
class LevelCollapsing{
	public:
		//type:
		typedef LevelCollapsing<IntType> type_this;
		typedef IntType type_int;
		typedef Liuze::Type_ArrayTemp<type_int> type_map;
		
		//Constructor:
		LevelCollapsing(void)=default;
		template<typename ArrayType>
		explicit LevelCollapsing(ArrayType const & map):m_map(map.size()){
			type_int i;
			for(i=0;i<m_map.size();++i) m_map[i]=map[i];
		}
		LevelCollapsing(type_this const & coll):m_map(coll.m_map){}
		LevelCollapsing(type_this && coll):m_map(std::move(coll.m_map)){}
		
		//Destructor:
		~LevelCollapsing(void)=default;
		
		//operator:
		type_this & operator =(type_this const & coll){
			if(&coll==this) return *this;
			m_map=coll.m_map;
			return *this;
		}
		type_this & operator =(type_this && coll){
			m_map=std::move(coll.m_map);
			return *this;
		}
		template<typename ArrayType>
		type_this & operator =(ArrayType const & map){
			m_map.resize(map.size());
			type_int i;
			for(i=0;i<m_map.size();++i) m_map[i]=map[i];
			return *this;
		}
		type_int operator ()(type_int const & x) const {
			return x>=type_int(m_map.size()) ? x : m_map[x];
		}
		type_int const & operator [](type_int const & x) const {
			return m_map[x];
		}
		type_int & operator [](type_int const & x){
			return m_map[x];
		}
		type_this operator *(type_this const & coll) const {
			type_this res;
			res.m_map.resize(coll.m_map.size());
			type_int i;
			for(i=0;i<res.m_map.size();++i) res.m_map[i]=(*this)(coll(i));
			return res;
		}
		type_this & operator *=(type_this const & coll){
			std::swap(((*this)*coll).m_map,this->m_map);
			return *this;
		}
		
		//other function:
		type_int size() const {
			return m_map.size();
		}
		type_this & resize(type_int const & size_new){
			m_map.resize(size_new);
			return *this;
		}
		Liuze::Type_ArrayTemp<type_int> sizes() const {
			typedef Liuze::Type_ArrayTemp<type_int> type_res;
			if(this->size()==type_int(0)) return type_res(0);
			type_int n_group;
			n_group=*(std::min_element(m_map.begin(),m_map.end()));
			if(n_group<type_int(0)) return type_res(0);
			n_group=*(std::max_element(m_map.begin(),m_map.end()));
			if(n_group>=this->size()) return type_res(0);
			++n_group;
			type_res res(n_group,type_int(0));
			type_int i;
			for(i=0;i<m_map.size();++i) ++res[m_map[i]];
			return res;
		}
		type_map & map(){
			return m_map;
		}
		type_map const & map() const {
			return m_map;
		}
		Liuze::Type_Bool Is_proper() const {
			if(m_map.size()==0) return true;
			Liuze::Type_ArrayTemp<type_int> size_group=this->sizes();
			if(size_group.size()==0) return false;
			type_int i;
			for(i=0;i<size_group.size();++i){
				if(size_group[i]==type_int(0)) return false;
			}
			return true;
		}
		type_this & Set_identity(type_int size_total=type_int(0)){
			if(size_total<=type_int(0)){
				if(m_map.size()==0) return *this;
			} else {
				m_map.resize(size_total);
			}
			type_int i;
			for(i=0;i<m_map.size();++i) m_map[i]=i;
			return *this;
		}
		type_this & Set_collapsing_std(type_int const & size_group,type_int size_total=type_int(0)){
			if(size_group<=type_int(0)) return *this;
			if(size_total<=type_int(0)){
				if(m_map.size()==0) return *this;
			} else {
				m_map.resize(size_total);
			}
			type_int i=0,ig=0,iig=0;
			while(i<m_map.size()){
				m_map[i]=ig;
				++i;
				++iig;
				if(iig==size_group){
					iig=type_int(0);
					++ig;
				}
			}
			return *this;
		}
		type_this & Set_collapsing_costd(type_int const & n_group,type_int size_total=type_int(0)){
			if(n_group<=type_int(0)) return *this;
			if(size_total<=type_int(0)){
				if(m_map.size()==0) return *this;
			} else {
				m_map.resize(size_total);
			}
			type_int i=0,ig=0;
			while(i<m_map.size() && ig<n_group){
				while(i<m_map.size()){
					m_map[i]=ig;
					i+=n_group;
				}
				++ig;
				i=ig;
			}
			return *this;
		}
		type_this & Set_reordering_costd(type_int const & n_group,type_int size_total=type_int(0)){
			if(n_group<=type_int(0)) return *this;
			if(size_total<=type_int(0)){
				if(m_map.size()==0) return *this;
			} else {
				m_map.resize(size_total);
			}
			type_int i=0,ig=0,j=0;
			while(i<m_map.size() && ig<n_group){
				while(i<m_map.size()){
					m_map[i]=j++;
					i+=n_group;
				}
				++ig;
				i=ig;
			}
			return *this;
		}
		type_this & Set_reverse(){
			type_map map(m_map.size());
			type_int i,ir;
			for(i=0,ir=m_map.size()-1;i<m_map.size();++i,--ir) map[i]=m_map[ir];
			std::swap(map,m_map);
			return *this;
		}
		type_this & Set_shiftleft(type_int const & offset){
			if(offset==type_int(0)) return *this;
			if(offset<type_int(0)) return this->Set_shiftright(-offset);
			type_map map(m_map.size());
			type_int i,is=Liuze::math::mod(offset,type_int(m_map.size()));
			for(i=0;i<m_map.size();++i){
				map[i]=m_map[is];
				++is;
				if(is==type_int(m_map.size())) is=type_int(0);
			}
			std::swap(map,m_map);
			return *this;
		}
		type_this & Set_shiftright(type_int const & offset){
			if(offset==type_int(0)) return *this;
			if(offset<type_int(0)) return this->Set_shiftleft(-offset);
			type_map map(m_map.size());
			type_int i,is=Liuze::math::mod(offset,type_int(m_map.size()));
			for(i=0;i<m_map.size();++i){
				map[is]=m_map[i];
				++is;
				if(is==type_int(m_map.size())) is=type_int(0);
			}
			std::swap(map,m_map);
			return *this;
		}
		template<typename FieldType,
			typename=typename std::enable_if<
				std::is_integral<typename FieldType::value_type>::value &&
				std::is_same<typename Liuze::math::GaloisField::tag_type_space<FieldType>::type,
				Liuze::math::GaloisField::tag_type_space_tabulated>::value
			>::type>
		type_this & Set_collapsing_GFlog(FieldType const & field){
			if(field.size()<=type_int(0)){
				m_map.resize(0);
				return *this;
			}
			if(field.size()==type_int(1)){
				m_map.assign(1,type_int(0));
				return *this;
			}
			m_map.resize(field.size());
			m_map[0]=type_int(0);
			type_int i;
			for(i=1;i<field.size();++i){
				m_map[i]=field.space_act().table_log(i-type_int(1));
			}
			return *this;
		}
	protected:
		//member value:
		type_map m_map;
}; //class LevelCollapsing;

//fun: Get_COAfrac_collapser:
//tran_level:
//0: origin; 1: ordered_post0; 2: collapsed_post0; 3: ordered_pre0; 4: collapsed_pre0:
template<typename FieldType,typename IntType=typename FieldType::value_type,
	typename=typename std::enable_if<
		std::is_integral<typename FieldType::value_type>::value &&
		std::is_integral<IntType>::value &&
		std::is_same<typename Liuze::math::GaloisField::tag_type_space<FieldType>::type,
		Liuze::math::GaloisField::tag_type_space_tabulated>::value
	>::type>
LevelCollapsing<IntType> Get_COAfrac_collapser
(FieldType const & field,IntType n_run_frac=0,Liuze::Type_Stat tran_level=0){
	typedef typename FieldType::value_type type_val;
	typedef typename FieldType::type_element type_eleGF;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	typedef LevelCollapsing<type_int> type_coll;
	type_int n_comp=field.size();
	type_coll coll;
	if(n_comp<=type_int(0)) return coll;
	if(n_comp==type_int(1)) return coll.Set_identity(n_comp);
	if(n_run_frac<=type_int(0)) n_run_frac=1;
	if(n_comp%n_run_frac==type_int(0)){
		if(tran_level==Type_Stat(2) || tran_level==Type_Stat(4)){
			coll.Set_collapsing_std(n_run_frac,n_comp);
		} else {
			coll.Set_identity(n_comp);
		}
		return coll;
	} else if((n_comp-type_int(1))%n_run_frac==type_int(0)){
		type_int n_mul=(n_comp-type_int(1))/n_run_frac,i_comp;
		switch(tran_level){
			case Type_Stat(1):
				coll=type_coll().Set_reordering_costd(n_mul,n_comp-type_int(1))
					*type_coll().Set_collapsing_GFlog(field);
				coll[0]=n_comp-type_int(1);
				return coll;
			case Type_Stat(2):
				coll=type_coll().Set_collapsing_costd(n_mul,n_comp-type_int(1))
					*type_coll().Set_collapsing_GFlog(field);
				coll[0]=n_mul;
				return coll;
			case Type_Stat(3):
				coll=type_coll().Set_reordering_costd(n_mul,n_comp-type_int(1))
					*type_coll().Set_collapsing_GFlog(field);
				for(i_comp=1;i_comp<n_comp;++i_comp) ++coll[i_comp];
				return coll;
			case Type_Stat(4):
				coll=type_coll().Set_collapsing_costd(n_mul,n_comp-type_int(1))
					*type_coll().Set_collapsing_GFlog(field);
				for(i_comp=1;i_comp<n_comp;++i_comp) ++coll[i_comp];
				return coll;
			default:
				return coll.Set_identity(n_comp);
		}
	} else {
		return coll;
	}
} //fun: Get_COAfrac_collapser;

//fun: Get_COA2frac:
template<typename FieldType,typename IntType=typename FieldType::value_type,
	typename=typename std::enable_if<
		std::is_integral<typename FieldType::value_type>::value &&
		std::is_integral<IntType>::value &&
		std::is_same<typename Liuze::math::GaloisField::tag_type_space<FieldType>::type,
		Liuze::math::GaloisField::tag_type_space_tabulated>::value
	>::type>
Liuze::Type_MatTemp<IntType> Get_COA2frac(FieldType const & field,IntType n_run_frac=0){
	typedef typename FieldType::value_type type_val;
	typedef typename FieldType::type_element type_eleGF;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_comp=field.size();
	if(n_comp<=type_int(0)) return type_res(0,0);
	if(n_comp==type_int(1)) return type_res::Zero(1,1);
	if(n_run_frac<=type_int(0)) n_run_frac=1;
	if(n_comp%n_run_frac==type_int(0)){
		type_int n_add=n_comp/n_run_frac;
		type_int n_run=n_add*(n_comp-type_int(1));
		type_res res(n_run,n_comp);
		type_int ic,ir0,ir1,im;
		type_eleGF t_prod;
		for(ir1=0,im=1;ir1<n_run;ir1+=n_add,++im){
			for(ic=0;ic<n_comp;++ic){
				t_prod=field(im)*field(ic);
				res(ir1,ic)=type_int(t_prod);
				for(ir0=1;ir0<n_add;++ir0){
					res(ir1+ir0,ic)=type_int(t_prod+field(ir0*n_run_frac));
				}
			}
		}
		return res;
	} else if((n_comp-type_int(1))%n_run_frac==type_int(0)){
		type_int n_mul=(n_comp-type_int(1))/n_run_frac;
		type_int n_run=n_comp*n_mul;
		type_res res(n_run,n_comp);
		type_int ic,ir0,ir1;
		type_eleGF eleGFprim(type_int(field.space_act().primitive()),field);
		type_eleGF t_sum;
		for(ir0=0;ir0<n_comp;++ir0){
			for(ic=0;ic<n_comp;++ic){
				t_sum=field(ir0)+field(ic);
				res(ir0,ic)=type_int(t_sum);
				for(ir1=1;ir1<n_mul;++ir1){
					res(ir1*n_comp+ir0,ic)=type_int(eleGFprim.pow(ir1)*t_sum);
				}
			}
		}
		return res;
	} else {
		return type_res(0,0);
	}
} //fun: Get_COA2frac;

//fun: Get_COA3frac:
template<typename FieldType,typename IntType=typename FieldType::value_type,
	typename=typename std::enable_if<std::is_integral<typename FieldType::value_type>::value &&
	std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_COA3frac
(FieldType const & field_init,IntType n_run_frac=0,Liuze::Type_MatTemp<IntType> * ptr_des_init=NULL){
	typedef typename FieldType::value_type type_val;
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_res;
	type_int n_comp_init=field_init.size();
	if(n_comp_init<=type_int(1)) return type_res(0,0);
	type_int i_row,i_col;
	type_int n_comp=n_comp_init+type_int(1);
	type_res augment_scheme(n_comp,n_comp);
	for(i_row=0;i_row<n_comp;++i_row) augment_scheme(i_row,i_row)=n_comp_init;
	for(i_col=1;i_col<n_comp;++i_col){
		augment_scheme(0,i_col)=type_int(field_init(i_col-type_int(1)));
	}
	augment_scheme(1,0)=type_int(field_init(0));
	for(i_col=2;i_col<n_comp;++i_col){
		augment_scheme(1,i_col)=type_int(field_init(i_col-type_int(1)).inv());
	}
	for(i_row=2;i_row<n_comp;++i_row){
		for(i_col=0;i_col<type_int(2);++i_col){
			augment_scheme(i_row,i_col)=type_int(field_init(i_col));
		}
		for(;i_col<n_comp;++i_col){
			if(i_col==i_row) continue;
			augment_scheme(i_row,i_col)=type_int(
				field_init(i_row-type_int(1))*
				Liuze::math::inverse(field_init(i_row-type_int(1))-field_init(i_col-type_int(1))));
		}
	}
	return Get_COA_addComponent_custom
		((ptr_des_init ? *ptr_des_init : Get_COA2frac(field_init,n_run_frac)),augment_scheme);
} //fun: Get_COA3frac;

//fun: Get_mat_PofC_isomorph_symbperm_rand:
template<typename IntType=Liuze::Type_UInt,typename RandEngType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_MatTemp<IntType> Get_mat_PofC_isomorph_symbperm_rand
(Liuze::Type_MatTemp<IntType> const & des_PofC,LevelCollapsing<IntType> const & strat,
RandEngType & rand_eng){
	typedef IntType type_int;
	typedef Eigen::RowVectorX<type_int> type_row_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	typedef Liuze::Type_ArrayTemp<type_int> type_arr_int;
	type_int n_run=des_PofC.rows(),n_comp=des_PofC.cols();
	type_arr_int vec_stratsize=strat.sizes();
	type_int n_strat=vec_stratsize.size();
	type_int i_run,i_strat,i_comp;
	type_arr_int vec_offset(n_strat);
	for(i_comp=i_strat=0;i_strat<n_strat;++i_strat){
		vec_offset[i_strat]=i_comp;
		i_comp+=vec_stratsize[i_strat];
	}
	type_mat_int mat_res(n_run,n_comp);
	for(i_strat=0;i_strat<n_strat;++i_strat){
		mat_res.block(0,vec_offset[i_strat],n_run,vec_stratsize[i_strat])=
			(Get_OofARandDesign(n_run,vec_stratsize[i_strat],rand_eng).array()
			+vec_offset[i_strat]).matrix();
	}
	for(i_run=0;i_run<n_run;++i_run){
		mat_res.row(i_run)=
			Liuze::math::FiniteSymmetricGroup<type_int>::template compose<type_row_int>
			(type_row_int(mat_res.row(i_run)),type_row_int(des_PofC.row(i_run)),n_comp);
	}
	return mat_res;
} //fun: Get_mat_PofC_isomorph_symbperm_rand;

//fun: Cal_SO_numrun_min:
template<typename IntType,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
IntType Cal_SO_numrun_min
(IntType const & n_comp,IntType const & strth,
Liuze::Type_ArrayTemp<IntType> const & vec_stratsize){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	typedef Liuze::Type_ArrayTemp<type_int> type_arr_int;
	typedef unsigned long long type_lint;
	type_int n_strat=vec_stratsize.size();
	if(n_comp<=type_int(0) || strth<type_int(0) || strth>n_comp
	|| n_strat==type_int(0) || n_strat>n_comp){
		return 0;
	}
	if(strth==type_int(0)) return 1;
	type_int i,sum=0;
	for(i=0;i<n_strat;++i) sum+=vec_stratsize[i];
	if(sum!=n_comp) return 0;
	type_arr_int vec_stratsize_1(n_strat);
	for(i=0;i<n_strat;++i) vec_stratsize_1[i]=vec_stratsize[i]+type_int(1);
	Liuze::math::Comb_enum_cart<type_int> iter_possize(vec_stratsize_1);
	type_lint g(0),prod;
	while(iter_possize){
		sum=0;
		for(i=0;i<n_strat && sum<=strth;++i) sum+=(*iter_possize)[i];
		if(sum!=strth) goto GTL_nextpossize;
		prod=1;
		for(i=0;i<n_strat;++i){
			prod*=Liuze::math::Comb_num::perm<type_lint>(vec_stratsize[i],(*iter_possize)[i]);
		}
		g=EuclideanDivision<type_lint>(g,prod);
		GTL_nextpossize:
		++iter_possize;
	}
	return type_int(Liuze::math::Comb_num::perm<type_lint>(n_comp,strth)/g);
} //fun: Cal_SO_numrun_min;

//fun: Check_PofC_SO:
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
Liuze::Type_Bool Check_PofC_SO
(Liuze::Type_MatTemp<IntType> const & des,IntType const & strth,
LevelCollapsing<IntType> strat=LevelCollapsing<IntType>().Set_identity()){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	typedef Liuze::Type_ArrayTemp<type_int> type_arr_int;
	type_int n_run=des.rows(),n_comp=des.cols();
	if(n_run<=type_int(0) || n_comp<=type_int(0) || strth<type_int(0) || strth>n_comp
	|| strat.size()>n_comp){
		return false;
	}
	if(n_comp==type_int(1) || strth==type_int(0)) return true;
	if(!strat.Is_proper()) return false;
	type_arr_int strat_sizes=strat.sizes();
	type_int strat_n_group=strat_sizes.size();
	type_int i_run,i_pos;
	if(strat.size()<n_comp){
		i_run=strat.size();
		i_pos=strat_n_group;
		strat.resize(n_comp);
		for(;i_run<n_comp;++i_run,++i_pos) strat[i_run]=i_pos;
		strat_sizes=strat.sizes();
		strat_n_group=strat_sizes.size();
	}
	//argument checkings end.
	//calculating the repetition numbers of stratum tuples{
	type_int denom_n_rep_tup_pos=Liuze::math::Comb_num::perm(n_comp,strth);
	//the repetition numbers of stratum tuples:
	type_arr_int n_rep_tup_pos(Liuze::math::powint(strat_n_group,strth));
	Liuze::math::Comb_enum_pow<type_int> iter_pos(strat_n_group,strth);
	type_arr_int iter_pos_sizes(strat_n_group);
	i_pos=type_int(0);
	while(iter_pos.Is_stat_deref()){
		for(i_run=0;i_run<strat_n_group;++i_run) iter_pos_sizes[i_run]=type_int(0);
		for(i_run=0;i_run<strth;++i_run) ++iter_pos_sizes[(*iter_pos)[i_run]];
		n_rep_tup_pos[i_pos]=type_int(1);
		for(i_run=0;i_run<strat_n_group;++i_run){
			n_rep_tup_pos[i_pos]*=
				Liuze::math::Comb_num::perm(strat_sizes[i_run],iter_pos_sizes[i_run]);
		}
		n_rep_tup_pos[i_pos]*=n_run;
		if(n_rep_tup_pos[i_pos]%denom_n_rep_tup_pos!=type_int(0)) return false;
		n_rep_tup_pos[i_pos]/=denom_n_rep_tup_pos;
		++i_pos;
		++iter_pos;
	}
	//}(calculating the repetition numbers of stratum tuples)
	iter_pos_sizes.resize(strth);
	type_arr_int n_rep_tup_pos_act(n_rep_tup_pos.size());
	Liuze::math::Comb_enum_comb<type_int> iter_comp(n_comp,strth);
	while(iter_comp.Is_stat_deref()){
		for(i_pos=0;i_pos<n_rep_tup_pos_act.size();++i_pos) n_rep_tup_pos_act[i_pos]=type_int(0);
		for(i_run=0;i_run<n_run;++i_run){
			for(i_pos=0;i_pos<strth;++i_pos){
				iter_pos_sizes[i_pos]=strat(des(i_run,(*iter_comp)[i_pos]));
			}
			Liuze::math::Comb_enum_pow<type_int>::value_to_index
				(strat_n_group,iter_pos_sizes,i_pos);
			++n_rep_tup_pos_act[i_pos];
			if(n_rep_tup_pos_act[i_pos]>n_rep_tup_pos[i_pos]) return false;
		}
		++iter_comp;
	}
	return true;
} //fun: Check_PofC_SO;

//fun: Find_COA2_normSO3strat3:
template<typename FieldType,typename IntType=typename FieldType::value_type,
	typename=typename std::enable_if<
		std::is_integral<typename FieldType::value_type>::value &&
		std::is_integral<IntType>::value &&
		std::is_same<typename Liuze::math::GaloisField::tag_type_space<FieldType>::type,
		Liuze::math::GaloisField::tag_type_space_tabulated>::value
	>::type>
IntType Find_COA2_normSO3strat3(FieldType const & field){
	typedef IntType type_int;
	typedef Liuze::Type_MatTemp<type_int> type_mat_int;
	typedef Liuze::Type_ArrayTemp<type_int> type_arr_int;
	type_int n_comp=field.size();
	if(n_comp<=type_int(7) || n_comp%type_int(9)!=type_int(8)) return type_int(0);
	type_int res=0;
	type_int n_strat=3,size_strat0=(n_comp+type_int(1))/type_int(3)-type_int(1);
	type_int i_pos,i_pos_1;
	type_mat_int des=Get_COA2frac(field);
	LevelCollapsing<type_int> stratmap;
	stratmap.resize(n_comp);
	type_arr_int arr_pos_rest(n_comp-size_strat0-type_int(1));
	Liuze::math::Comb_enum_comb<type_int>
		iter_strat0(n_comp,size_strat0),iter_strat1(type_int(arr_pos_rest.size()),size_strat0);
	while(iter_strat0.Is_stat_deref()){
		for(i_pos=0;i_pos<n_comp;++i_pos) stratmap[i_pos]=type_int(3);
		for(i_pos=0;i_pos<size_strat0;++i_pos) stratmap[(*iter_strat0)[i_pos]]=type_int(0);
		for(i_pos=0;i_pos<n_comp;++i_pos){
			if(stratmap[i_pos]==type_int(3)) break;
		}
		stratmap[i_pos]=type_int(1);
		for(++i_pos,i_pos_1=0;i_pos_1<arr_pos_rest.size();++i_pos){
			if(stratmap[i_pos]==type_int(3)){
				arr_pos_rest[i_pos_1]=i_pos;
				++i_pos_1;
			}
		}
		iter_strat1=iter_strat1.begin();
		while(iter_strat1.Is_stat_deref()){
			for(i_pos=0;i_pos<arr_pos_rest.size();++i_pos){
				stratmap[arr_pos_rest[i_pos]]=type_int(2);
			}
			for(i_pos=0;i_pos<size_strat0;++i_pos){
				stratmap[arr_pos_rest[(*iter_strat1)[i_pos]]]=type_int(1);
			}
			if(Check_PofC_SO(des,type_int(3),stratmap)){
				cout<<stratmap.map()<<endl;
				++res;
			}
			++iter_strat1;
		}
		++iter_strat0;
	}
	return res;
} //fun: Find_COA2_normSO3strat3;

/****************************************************************************************************/
template<typename IntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<IntType>::value>::type>
class Searcher_COA1_SO2{
	public:
		//type:
		typedef Searcher_COA1_SO2<IntType> type_this;
		typedef IntType type_int;
		typedef Liuze::Type_MatTemp<type_int> type_mat_des;
		
		static Liuze::Type_Bool push_result_null_false(type_mat_des const &){
			return false;
		}
		
		static type_int search
		(type_int const & n_component,type_int size_strat_0=0,
		std::function<Liuze::Type_Bool(type_mat_des const &)> push_result=
		type_this::push_result_null_false){
			type_int res=0;
			if(n_component<=type_int(1)) return res;
			if(size_strat_0<=type_int(0) || size_strat_0>=n_component){
				size_strat_0=n_component/type_int(2);
			}
			type_int size_strat_1=n_component-size_strat_0; //the size of the other stratum;
			type_int constexpr n_type_tup_pos=4; //num_type_tuple_position;
			type_int i_col_1,i_tup_pos;
			Liuze::Type_ArrayTemp<type_int> n_rep_type_tup_pos(n_type_tup_pos);
			n_rep_type_tup_pos[0]=size_strat_0*(size_strat_0-type_int(1));
			n_rep_type_tup_pos[1]=size_strat_0*size_strat_1;
			n_rep_type_tup_pos[2]=n_rep_type_tup_pos[1];
			n_rep_type_tup_pos[3]=size_strat_1*(size_strat_1-type_int(1));
			i_col_1=n_component-type_int(1);
			for(i_tup_pos=0;i_tup_pos<n_type_tup_pos;++i_tup_pos){
				if(n_rep_type_tup_pos[i_tup_pos]%i_col_1==type_int(0)){
					n_rep_type_tup_pos[i_tup_pos]/=i_col_1;
				} else {
					return res;
				}
			}
			//checkings end.
			type_int i_col,i_tup_comp;
			std::function<type_int(type_int const &,type_int const &)> infun_id_tup_comp=
				[](type_int const & c0,type_int const & c1)->type_int{
					return (c1-type_int(1))*c1/type_int(2)+c0;
				};
			std::function<type_int(type_int const &,type_int const &)> infun_type_tup_pos=
				[&size_strat_0](type_int const & x0,type_int const & x1)->type_int{
					return (x0<size_strat_0 ? type_int(0) : type_int(2))
						+ (x1<size_strat_0 ? type_int(0) : type_int(1));
				};
			Liuze::Type_ArrayTemp<Liuze::Type_ArrayTemp<type_int> > n_rep_type_tup_comppos
				(n_component*(n_component-type_int(1))/type_int(2),
				Liuze::Type_ArrayTemp<type_int>(n_type_tup_pos,type_int(0)));
			for(i_col=1;i_col<n_component;++i_col){
				for(i_col_1=0;i_col_1<i_col;++i_col_1){
					i_tup_comp=infun_id_tup_comp(i_col_1,i_col);
					i_tup_pos=infun_type_tup_pos(i_col_1,i_col);
					++n_rep_type_tup_comppos[i_tup_comp][i_tup_pos];
					if(n_rep_type_tup_comppos[i_tup_comp][i_tup_pos]>n_rep_type_tup_pos[i_tup_pos]){
						return res;
					}
				}
			}
			type_mat_des LS(n_component,n_component); //the design;
			type_int i_row,i_row_1,i_pos;
			for(i_col=0;i_col<n_component;++i_col) LS(0,i_col)=i_col;
			for(i_row=0;i_row<n_component;++i_row) LS(i_row,0)=i_row;
			std::vector<bool> flag_pos_cand_row(n_component),flag_pos_cand_ent(n_component);
			i_row=1;
			i_col=1;
			flag_pos_cand_row.assign(n_component,true);
			flag_pos_cand_row[1]=false;
			while(true){
				//initialise the candidate values for the current entry{
				flag_pos_cand_ent=flag_pos_cand_row;
				for(i_row_1=0;i_row_1<i_row;++i_row_1){
					flag_pos_cand_ent[LS(i_row_1,i_col)]=false;
				}
				//}
				//find the smallest candidate value for the current entry:
				for(i_pos=0;i_pos<n_component;++i_pos){
					if(flag_pos_cand_ent[i_pos]) break;
				}
				//if there is no candidate value:
				if(i_pos==n_component) goto GTL_rollback_gloabal;
				//test the feasibility of the value of the current entry:
				GTL_feasiblity_test:
				LS(i_row,i_col)=i_pos;
				for(i_col_1=0;i_col_1<i_col;++i_col_1){
					i_tup_comp=infun_id_tup_comp(i_col_1,i_col);
					i_tup_pos=infun_type_tup_pos(LS(i_row,i_col_1),LS(i_row,i_col));
					++n_rep_type_tup_comppos[i_tup_comp][i_tup_pos];
					if(n_rep_type_tup_comppos[i_tup_comp][i_tup_pos]>n_rep_type_tup_pos[i_tup_pos]){
						goto GTL_rollback_local;
					}
				}
				//the value of the current entry is feasible.
				flag_pos_cand_row[i_pos]=false;
				++i_col;
				if(i_col==n_component){
					++i_row;
					if(i_row==n_component){
						++res;
						if(!push_result(LS)) return res;
						--i_row;
						--i_col;
						flag_pos_cand_row[i_pos]=true;
						i_col_1=i_col-type_int(1);
						goto GTL_rollback_local;
					}
					i_col=type_int(1);
					flag_pos_cand_row.assign(n_component,true);
					flag_pos_cand_row[i_row]=false;
				}
				continue;
				//rollback for the current entry:
				GTL_rollback_local:
				--n_rep_type_tup_comppos[i_tup_comp][i_tup_pos];
				while(i_col_1>type_int(0)){
					--i_col_1;
					--n_rep_type_tup_comppos[infun_id_tup_comp(i_col_1,i_col)]
						[infun_type_tup_pos(LS(i_row,i_col_1),LS(i_row,i_col))];
				}
				//rollback for earlier entries:
				GTL_rollback_gloabal:
				while(true){
					//find the next smallest candidate value for the current entry:
					for(++i_pos;i_pos<n_component;++i_pos){
						if(flag_pos_cand_ent[i_pos]) break;
					}
					if(i_pos<n_component) goto GTL_feasiblity_test;
					//there is no candidate value for the current entry.
					--i_col;
					if(i_col==type_int(0)){
						--i_row;
						if(i_row==type_int(0)) return res;
						i_col=n_component-type_int(1);
						flag_pos_cand_row.assign(n_component,false);
					}
					//remove the value for the current entry{
					for(i_col_1=0;i_col_1<i_col;++i_col_1){
						--n_rep_type_tup_comppos[infun_id_tup_comp(i_col_1,i_col)]
							[infun_type_tup_pos(LS(i_row,i_col_1),LS(i_row,i_col))];
					}
					i_pos=LS(i_row,i_col);
					flag_pos_cand_row[i_pos]=true;
					flag_pos_cand_ent=flag_pos_cand_row;
					for(i_row_1=0;i_row_1<i_row;++i_row_1){
						flag_pos_cand_ent[LS(i_row_1,i_col)]=false;
					}
					//}
				}
			}
			return res;
		} //fun: search;
		
		static type_int search_onecycle
		(type_int const & n_component,type_int size_strat_0=0,
		std::function<Liuze::Type_Bool(type_mat_des const &)> push_result=
		type_this::push_result_null_false){
			type_int res=0;
			if(n_component<=type_int(1)) return res;
			if(n_component==type_int(2)){
				type_mat_des des=type_mat_des::Ones(2,2);
				des(0,0)=des(1,1)=type_int(0);
				push_result(des);
				return type_int(1);
			}
			if(size_strat_0<=type_int(0) || size_strat_0>=n_component){
				size_strat_0=n_component/type_int(2);
			}
			type_int size_strat_1=n_component-size_strat_0; //the size of the other stratum;
			type_int constexpr n_type_tup_pos=4; //num_type_tuple_position;
			type_int i_col_1,i_tup_pos;
			Liuze::Type_ArrayTemp<type_int> n_rep_type_tup_pos(n_type_tup_pos);
			n_rep_type_tup_pos[0]=size_strat_0*(size_strat_0-type_int(1));
			n_rep_type_tup_pos[1]=size_strat_0*size_strat_1;
			n_rep_type_tup_pos[2]=n_rep_type_tup_pos[1];
			n_rep_type_tup_pos[3]=size_strat_1*(size_strat_1-type_int(1));
			i_col_1=n_component-type_int(1);
			for(i_tup_pos=0;i_tup_pos<n_type_tup_pos;++i_tup_pos){
				if(n_rep_type_tup_pos[i_tup_pos]%i_col_1==type_int(0)){
					n_rep_type_tup_pos[i_tup_pos]/=i_col_1;
				} else {
					return res;
				}
			}
			//checkings end.
			std::function<type_int(type_int const &,type_int const &)> infun_type_tup_pos=
				[&size_strat_0](type_int const & x0,type_int const & x1)->type_int{
					return (x0<size_strat_0 ? type_int(0) : type_int(2))
						+ (x1<size_strat_0 ? type_int(0) : type_int(1));
				};
			Liuze::Type_ArrayTemp<type_int> n_rep_type_tup_pos_act(n_type_tup_pos);
			type_int i_row,i_col;
			type_mat_des des(n_component,n_component);
			for(i_col=0;i_col<n_component;++i_col) des(0,i_col)=i_col;
			Eigen::RowVectorX<type_int> des_row(n_component);
			Liuze::math::Comb_enum_perm<type_int> iter_cycle_head(n_component-type_int(2));
			while(iter_cycle_head.Is_stat_deref()){
				i_row=n_component-type_int(3);
				for(i_col=0;i_col<i_row;++i_col){
					des_row[(*iter_cycle_head)[i_col]]=(*iter_cycle_head)[i_col+type_int(1)];
				}
				des_row[(*iter_cycle_head)[i_col]]=n_component-type_int(2);
				des_row[n_component-type_int(2)]=n_component-type_int(1);
				des_row[n_component-type_int(1)]=(*iter_cycle_head)[0];
				i_row=des_row[0];
				des.row(i_row)=des_row;
				for(i_col=2;i_col<n_component;++i_col){
					des_row=Liuze::math::FiniteSymmetricGroup<type_int>::
						template compose<Eigen::RowVectorX<type_int> >(des_row,des.row(i_row));
					des.row(des_row[0])=des_row;
				}
				//checking feasibility of the current design `des':
				for(i_col=1;i_col<n_component;++i_col){
					for(i_col_1=0;i_col_1<i_col;++i_col_1){
						n_rep_type_tup_pos_act.assign(n_type_tup_pos,type_int(0));
						for(i_row=0;i_row<n_component;++i_row){
							i_tup_pos=infun_type_tup_pos(des(i_row,i_col_1),des(i_row,i_col));
							++n_rep_type_tup_pos_act[i_tup_pos];
							if(n_rep_type_tup_pos_act[i_tup_pos]>n_rep_type_tup_pos[i_tup_pos]){
								goto GTL_infeasible;
							}
						}
					}
				}
				//`des' is feasible.
				++res;
				if(!push_result(des)) return res;
				//`des' is infeasible:
				GTL_infeasible:
				++iter_cycle_head;
			}
			return res;
		} //fun: search_onecycle;
		
		static Liuze::Type_Bool check
		(type_mat_des const & design,type_int size_strat_0=0){
			type_int n_component=design.cols();
			if(n_component<=type_int(1) || design.rows()<=type_int(0)) return false;
			if(size_strat_0<=type_int(0) || size_strat_0>=n_component){
				size_strat_0=n_component/type_int(2);
			}
			type_int size_strat_1=n_component-size_strat_0; //the size of the other stratum;
			type_int constexpr n_type_tup_pos=4; //num_type_tuple_position;
			type_int i_col_1,i_tup_pos;
			Liuze::Type_ArrayTemp<type_int> n_rep_type_tup_pos(n_type_tup_pos);
			n_rep_type_tup_pos[0]=size_strat_0*(size_strat_0-type_int(1));
			n_rep_type_tup_pos[1]=size_strat_0*size_strat_1;
			n_rep_type_tup_pos[2]=n_rep_type_tup_pos[1];
			n_rep_type_tup_pos[3]=size_strat_1*(size_strat_1-type_int(1));
			i_col_1=n_component-type_int(1);
			for(i_tup_pos=0;i_tup_pos<n_type_tup_pos;++i_tup_pos){
				if(n_rep_type_tup_pos[i_tup_pos]%i_col_1==type_int(0)){
					n_rep_type_tup_pos[i_tup_pos]/=i_col_1;
				} else {
					return false;
				}
			}
			//checkings end.
			type_int i_row,i_col;
			std::function<type_int(type_int const &,type_int const &)> infun_type_tup_pos=
				[&size_strat_0](type_int const & x0,type_int const & x1)->type_int{
					return (x0<size_strat_0 ? type_int(0) : type_int(2))
						+ (x1<size_strat_0 ? type_int(0) : type_int(1));
				};
			Liuze::Type_ArrayTemp<type_int> n_rep_type_tup_pos_act(n_type_tup_pos);
			for(i_col=1;i_col<n_component;++i_col){
				for(i_col_1=0;i_col_1<i_col;++i_col_1){
					n_rep_type_tup_pos_act.assign(n_type_tup_pos,type_int(0));
					for(i_row=0;i_row<design.rows();++i_row){
						i_tup_pos=infun_type_tup_pos(design(i_row,i_col_1),design(i_row,i_col));
						++n_rep_type_tup_pos_act[i_tup_pos];
						if(n_rep_type_tup_pos_act[i_tup_pos]>n_rep_type_tup_pos[i_tup_pos]){
							return false;
						}
					}
				}
			}
			return true;
		} //fun: check;
}; //class Searcher_COA1_SO2;
/****************************************************************************************************/

#endif //#ifndef LZ_DEF_hpp_OofA_des_strat