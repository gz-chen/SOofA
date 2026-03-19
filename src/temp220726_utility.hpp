#ifndef LZ_DEF_hpp_temp220726_utility
#define LZ_DEF_hpp_temp220726_utility 202107L

#include<cstdarg>
#include<initializer_list>
#include<Eigen/SVD>
#include<LZ_H/math/function.hpp>

/****************************************************************************************************/
//convertions between Comb_enum and value vector{
template<typename ResType=void,typename UIntType=Liuze::Type_UInt,
	typename RetType=typename std::conditional<std::is_arithmetic<ResType>::value,ResType,UIntType>::type,
	typename=typename std::enable_if<std::is_arithmetic<RetType>::value &&
	std::is_integral<UIntType>::value>::type>
RetType Comb_enum_comb_index
(UIntType const & num_total,typename Liuze::math::Comb_enum_comb<UIntType>::type_val const & val){
	RetType index;
	Liuze::math::Comb_enum_comb<UIntType>::value_to_index(num_total,val,index);
	return index;
} //fun: Comb_enum_comb_index;
template<typename UIntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<UIntType>::value>::type>
typename Liuze::math::Comb_enum_comb<UIntType>::type_val
Comb_enum_comb_val
(UIntType const & num_total,UIntType const & num_select,UIntType const & index){
	typename Liuze::math::Comb_enum_comb<UIntType>::type_val val;
	Liuze::math::Comb_enum_comb<UIntType>::index_to_value(num_total,num_select,index,val);
	return val;
} //fun: Comb_enum_comb_val;

template<typename ResType=void,typename UIntType=Liuze::Type_UInt,
	typename RetType=typename std::conditional<std::is_arithmetic<ResType>::value,ResType,UIntType>::type,
	typename=typename std::enable_if<std::is_arithmetic<RetType>::value &&
	std::is_integral<UIntType>::value>::type>
RetType Comb_enum_perm_index
(typename Liuze::math::Comb_enum_perm<UIntType>::type_val const & val){
	RetType index;
	Liuze::math::Comb_enum_perm<UIntType>::value_to_index(val,index);
	return index;
} //fun: Comb_enum_perm_index;
template<typename UIntType=Liuze::Type_UInt,
	typename=typename std::enable_if<std::is_integral<UIntType>::value>::type>
typename Liuze::math::Comb_enum_perm<UIntType>::type_val
Comb_enum_perm_val
(UIntType const & num_order,UIntType const & index){
	typename Liuze::math::Comb_enum_perm<UIntType>::type_val val;
	Liuze::math::Comb_enum_perm<UIntType>::index_to_value(num_order,index,val);
	return val;
} //fun: Comb_enum_perm_val;
//}(convertions between Comb_enum and value vector)
/****************************************************************************************************/
//class Matrix_Inverse_MoorePenrose{
template<typename RealType=Liuze::Type_Real,
	typename=typename std::enable_if<std::is_arithmetic<RealType>::value>::type>
class Matrix_Inverse_MoorePenrose{
	public:
		typedef Matrix_Inverse_MoorePenrose<RealType> type_this;
		typedef Liuze::Type_Size type_size;
		typedef RealType type_real;
		typedef Liuze::Type_MatTemp<type_real> type_matrix;
		typedef Eigen::DiagonalMatrix<type_real,Eigen::Dynamic> type_matrix_diag;
		
		Matrix_Inverse_MoorePenrose()=default;
		template<typename EigenDerived,
			typename=typename std::enable_if<std::is_same<type_real,
			typename Eigen::DenseBase<EigenDerived>::Scalar>::value>::type>
		Matrix_Inverse_MoorePenrose(Eigen::DenseBase<EigenDerived> const & X){
			this->compute(X);
		}
		
		Eigen::BDCSVD<type_matrix> const & Get_SVD() const {
			return solver;
		}
		type_matrix Get_inverse() const {
			type_matrix res=
				solver.matrixV()*type_matrix_diag(diag_inv)*solver.matrixU().transpose();
			return tag_fat ? res.transpose() : res;
		}
		type_matrix Get_projM() const {
			type_matrix const & U= tag_fat ? solver.matrixV() : solver.matrixU();
			return U*type_matrix_diag(diag_proj)*U.transpose();
		}
		type_matrix Get_projMc() const {
			type_matrix const & U= tag_fat ? solver.matrixV() : solver.matrixU();
			return type_matrix::Identity(U.rows(),U.rows())
				-U*type_matrix_diag(diag_proj)*U.transpose();
		}
		type_matrix Get_singularValues() const {
			return solver.singularValues();
		}
		type_matrix Get_singularValues_nonzero() const {
			return vec_singularValues_nonzero;
		}
		type_size Get_dim_embed() const {
			return dim_embed;
		}
		type_size Get_dim() const {
			return dim;
		}
		type_real Get_det() const {
			return solver.singularValues().prod();
		}
		type_real Get_det_nonzero() const {
			return vec_singularValues_nonzero.prod();
		}
		
		template<typename EigenDerived,
			typename=typename std::enable_if<std::is_same<type_real,
			typename Eigen::DenseBase<EigenDerived>::Scalar>::value>::type>
		type_this & compute(Eigen::DenseBase<EigenDerived> const & X){
			tag_fat=(X.rows()<X.cols());
			solver.compute((tag_fat ? type_matrix(X).transpose() : type_matrix(X)),
				Eigen::ComputeThinU|Eigen::ComputeFullV);
			dim_embed=std::min(X.rows(),X.cols());
			dim=0;
			Eigen::VectorX<type_real> vec_singularValues=solver.singularValues();
			diag_proj.resize(dim_embed);
			diag_inv.resize(dim_embed);
			type_size i,i1;
			for(i=0;i<dim_embed;++i){
				if(Liuze::math::abs(vec_singularValues[i])
				>type_real(LZ_DEF_const_default_precision)){
					diag_inv[i]=type_real(1)/vec_singularValues[i];
					diag_proj[i]=type_real(1);
					++dim;
				} else {
					diag_inv[i]=diag_proj[i]=type_real(0);
				}
			}
			vec_singularValues_nonzero.resize(dim);
			for(i1=i=0;i1<dim && i<dim_embed;++i){
				if(diag_proj[i]!=type_real(0)){
					vec_singularValues_nonzero[i1]=vec_singularValues[i];
					++i1;
				}
			}
			return *this;
		}
		
		template<typename EigenDerived,
			typename=typename std::enable_if<std::is_same<type_real,
			typename Eigen::DenseBase<EigenDerived>::Scalar>::value>::type>
		type_matrix operator ()(Eigen::DenseBase<EigenDerived> const & X){
			return this->compute(X).Get_inverse();
		}
	protected:
		Eigen::BDCSVD<type_matrix> solver;
		Liuze::Type_Bool tag_fat;
		type_size dim_embed,dim;
		Eigen::VectorX<type_real> vec_singularValues_nonzero,diag_inv,diag_proj;
}; //class Matrix_Inverse_MoorePenrose;
//}(class Matrix_Inverse_MoorePenrose)
/****************************************************************************************************/
//class Calculator_norm_L_stream{
template<typename ValType=Liuze::Type_Real,
	typename RealType=Liuze::Type_Real,typename SizeType=Liuze::Type_Size,
	typename=typename std::enable_if<std::is_integral<SizeType>::value &&
	std::is_arithmetic<RealType>::value>::type>
class Calculator_norm_L_stream{
	public:
		typedef Calculator_norm_L_stream<ValType,RealType,SizeType,void> type_this;
		typedef ValType type_val;
		typedef RealType type_real;
		typedef SizeType type_size;
		typedef Liuze::Type_Bool type_bool;
		
		static constexpr type_real val_power_infinity=std::numeric_limits<type_real>::infinity();
		
		Calculator_norm_L_stream(type_real power=2,type_bool absolute_value=true):
		val_power(power>type_real(0) ? power :
		(power==type_real(0) ? type_this::val_power_infinity : type_real(2))),
		val_abs(absolute_value),val_len(0){
		}
		
		type_this & operator <<(type_val x){
			if(val_len==type_size(0)){
				val_norm=x-x;
			}
			if(val_abs) x=Liuze::math::abs(x);
			if(val_power==type_this::val_power_infinity){
				if(val_len==type_size(0) || x>val_norm) val_norm=x;
			} else {
				val_norm+=pow(x,val_power);
			}
			++val_len;
			return *this;
		}
		
		type_val Get_Norm() const {
			if(val_len==type_size(0)){
				type_val x;
				return x-x;
			}
			if(val_power==type_this::val_power_infinity){
				return val_norm;
			}
			if(val_power>type_real(1)){
				return pow(val_norm,type_real(1)/val_power);
			}
			return val_norm;
		}
		type_val Get_norm() const {
			if(val_len==type_size(0)){
				type_val x;
				return x-x;
			}
			if(val_power==type_this::val_power_infinity){
				return val_norm;
			}
			if(val_power>type_real(1)){
				return pow(val_norm/type_real(val_len),type_real(1)/val_power);
			}
			return val_norm/type_real(val_len);
		}
		type_val Get_sum() const {
			if(val_len==type_size(0)){
				type_val x;
				return x-x;
			}
			return val_norm;
		}
		type_val Get_mean() const {
			if(val_len==type_size(0)){
				type_val x;
				return x-x;
			}
			if(val_power==type_this::val_power_infinity){
				return val_norm;
			}
			return val_norm/type_real(val_len);
		}
		
		void Set_init(){
			val_len=type_size(0);
		}
		void Set_init(type_real const & power,type_bool absolute_value=true){
			val_power= power>type_real(0) ? power :
				(power==type_real(0) ? type_this::val_power_infinity : val_power);
			val_abs=absolute_value;
			this->Set_init();
		}
	protected:
		type_real val_power;
		type_bool val_abs;
		type_val val_norm;
		type_size val_len;
}; //class Calculator_norm_L_stream;
//}(Calculator_norm_L_stream)
/****************************************************************************************************/
//class Iterator_MatRowTranspose{
template<typename MatType,
	typename=typename std::enable_if<std::is_same<Liuze::Type_MatTemp<typename MatType::value_type>,
	typename std::remove_cv<MatType>::type>::value>::type>
class Iterator_MatRowTranspose{
	public:
		//types:
		typedef Iterator_MatRowTranspose<MatType> type_this;
		typedef typename std::conditional<std::is_const<MatType>::value,
				Eigen::Transpose<Eigen::Block<MatType,1,Eigen::Dynamic,false> const> const,
				Eigen::Transpose<Eigen::Block<MatType,1,Eigen::Dynamic,false> >
			>::type
			type_val;
		typedef typename MatType::Index type_size;
		//types: for iterator:
		typedef std::iterator<std::random_access_iterator_tag,type_val> iterator_base;
		typedef typename iterator_base::iterator_category iterator_category;
		typedef typename iterator_base::value_type value_type;
		typedef typename iterator_base::difference_type difference_type;
		typedef typename iterator_base::pointer pointer;
		typedef typename iterator_base::reference reference;
		
		//Constructor:
		Iterator_MatRowTranspose(void):ptr_mat(NULL),id_row(0){}
		explicit Iterator_MatRowTranspose(MatType & matrix):
		ptr_mat(std::addressof(matrix)),id_row(type_size(0)){
		}
		Iterator_MatRowTranspose(MatType & matrix,type_size const & id):
		ptr_mat(std::addressof(matrix)),id_row(id){
		}
		Iterator_MatRowTranspose(MatType && matrix)=delete;
		Iterator_MatRowTranspose(type_this const & iter):ptr_mat(iter.ptr_mat),id_row(iter.id_row){}
		
		//Destructor:
		~Iterator_MatRowTranspose(void)=default;
		
		//operator:
		type_this & operator =(type_this const & iter){
			if(&iter==this) return *this;
			ptr_mat=iter.ptr_mat;
			id_row=iter.id_row;
			return *this;
		}
		type_this & operator =(MatType & matrix){
			ptr_mat=std::addressof(matrix);
			id_row=type_size(0);
			return *this;
		}
		Liuze::Type_Bool operator ==(type_this const & iter) const {
			return (ptr_mat==iter.ptr_mat && id_row==iter.id_row);
		}
		Liuze::Type_Bool operator !=(type_this const & iter) const {
			return !(*this==iter);
		}
		Liuze::Type_Bool operator <(type_this const & iter) const {
			return ptr_mat==iter.ptr_mat ? id_row<iter.id_row : false;
		}
		Liuze::Type_Bool operator <=(type_this const & iter) const {
			return ptr_mat==iter.ptr_mat ? id_row<=iter.id_row : false;
		}
		Liuze::Type_Bool operator >(type_this const & iter) const {
			return ptr_mat==iter.ptr_mat ? id_row>iter.id_row : false;
		}
		Liuze::Type_Bool operator >=(type_this const & iter) const {
			return ptr_mat==iter.ptr_mat ? id_row>=iter.id_row : false;
		}
		type_val operator *() const {
			return (*ptr_mat).row(id_row).transpose();
		}
		type_val* operator ->() const {
			return std::addressof(**this);
		}
		type_this & operator ++(){
			if(ptr_mat) ++id_row;
			return *this;
		}
		type_this operator ++(int){
			type_this res=*this;
			++(*this);
			return res;
		}
		type_this & operator --(){
			if(ptr_mat) --id_row;
			return *this;
		}
		type_this operator --(int){
			type_this res=*this;
			--(*this);
			return res;
		}
		type_this operator +(difference_type const & offset) const {
			if(!ptr_mat) return *this;
			type_this res=*this;
			res.id_row+=offset;
			return res;
		}
		type_this operator -(difference_type const & offset) const {
			if(!ptr_mat) return *this;
			type_this res=*this;
			res.id_row-=offset;
			return res;
		}
		difference_type operator -(type_this const & iter) const {
			return ptr_mat==iter.ptr_mat ? difference_type(id_row-iter.id_row) : difference_type(0);
		}
		type_this & operator +=(difference_type const & offset){
			if(ptr_mat) id_row+=offset;
			return *this;
		}
		type_this & operator -=(difference_type const & offset){
			if(ptr_mat) id_row-=offset;
			return *this;
		}
		type_val operator [](difference_type const & offset) const {
			return (*ptr_mat).row(id_row+offset).transpose();
		}
		
		//Is:
		Liuze::Type_Bool Is_deref() const {
			return (ptr_mat && id_row>=type_size(0) && id_row<ptr_mat->rows());
		}
		
		//begin and end:
		type_this begin() const {
			return ptr_mat ? type_this(*ptr_mat,type_size(0)) : type_this();
		}
		type_this end() const {
			return ptr_mat ? type_this(*ptr_mat,ptr_mat->rows()) : type_this();
		}
		
		//Get:
		type_size size() const {
			return ptr_mat ? type_size(ptr_mat->rows()) : type_size(0);
		}
		
		//Set:
		void assign(MatType & matrix,type_size id=type_size(0)){
			ptr_mat=std::addressof(matrix);
			id_row=id;
		}
	protected:
		//member value:
		MatType * ptr_mat;
		type_size id_row;
}; //class Iterator_MatRowTranspose;

template<typename MatType>
inline Iterator_MatRowTranspose<MatType> operator +
(typename Iterator_MatRowTranspose<MatType>::difference_type const & offset,
Iterator_MatRowTranspose<MatType> const & iter){
	return iter+offset;
}
//}(Iterator_MatRowTranspose)
/****************************************************************************************************/
//ostream matters{
class fstream_null: public std::ostream{
	public:
		fstream_null(...){}
		~fstream_null()=default;
		fstream_null & open(...){return *this;}
		fstream_null & close(...){return *this;}
		fstream_null & operator <<(std::ostream&(*)(std::ostream&)){return *this;}
		template<typename Tp>
		fstream_null & operator <<(Tp && x){return *this;}
};
class fstream_std: public std::ostream{
	public:
		fstream_std(...){}
		~fstream_std()=default;
		fstream_std & open(...){return *this;}
		fstream_std & close(...){return *this;}
		fstream_std & operator <<(std::ostream&(*x)(std::ostream&)){
			::std::cout<<x;
			return *this;
		}
		template<typename Tp>
		fstream_std & operator <<(Tp && x){
			::std::cout<<x;
			return *this;
		}
};
template<Liuze::Type_Size size_buff=256>
class multi_ostream{
	public:
		//type:
		typedef multi_ostream type_this;
		typedef Liuze::Type_Size type_size;
		
		//Constructor:
		multi_ostream(void)=default;
		multi_ostream(Liuze::Type_ArrayTemp<std::ostream*> const & ptr_out_seq):
		ptr_out(ptr_out_seq){
		}
		multi_ostream(std::initializer_list<std::ostream*> const & ptr_out_init_list):
		ptr_out(ptr_out_init_list){
		}
		multi_ostream(type_this const & out):
		ptr_out(out.ptr_out){
		}
		multi_ostream(type_this && out):
		ptr_out(std::move(out.ptr_out)){
		}
		
		//Destructor:
		~multi_ostream(void)=default;
		
		//operator:
		type_this & operator =(type_this const & out){
			if(std::addressof(out)==this) return *this;
			ptr_out=out.ptr_out;
			return *this;
		}
		type_this & operator =(type_this && out){
			ptr_out=std::move(out.ptr_out);
			return *this;
		}
		type_this & operator =(Liuze::Type_ArrayTemp<std::ostream*> const & ptr_out_seq){
			ptr_out=ptr_out_seq;
			return *this;
		}
		type_this & operator =(std::initializer_list<std::ostream*> const & ptr_out_init_list){
			ptr_out=ptr_out_init_list;
			return *this;
		}
		template<typename Tp>
		type_this & operator <<(Tp const & obj){
			for(type_size i=0;i<ptr_out.size();++i){
				if(!(ptr_out[i])) continue;
				(*(ptr_out[i]))<<obj;
			}
			return *this;
		}
		type_this & operator <<(std::ios_base & (*pf)(std::ios_base &)){
			for(type_size i=0;i<ptr_out.size();++i){
				if(!(ptr_out[i])) continue;
				(*(ptr_out[i]))<<(*pf);
			}
			return *this;
		}
		type_this & operator <<(std::ios & (*pf)(std::ios &)){
			for(type_size i=0;i<ptr_out.size();++i){
				if(!(ptr_out[i])) continue;
				(*(ptr_out[i]))<<(*pf);
			}
			return *this;
		}
		type_this & operator <<(std::ostream & (*pf)(std::ostream &)){
			for(type_size i=0;i<ptr_out.size();++i){
				if(!(ptr_out[i])) continue;
				(*(ptr_out[i]))<<(*pf);
			}
			return *this;
		}
		std::ostream * & operator [](type_size const & id){
			return ptr_out[id];
		}
		
		//other function:
		std::ostream & at(type_size const & id){
			return (id>=type_size(0) && id<ptr_out.size()) ?
				(ptr_out[id] ? std::cout : *(ptr_out[id])) : std::cout;
		}
		type_this & flush() const {
			for(type_size i=0;i<ptr_out.size();++i){
				if(!(ptr_out[i])) continue;
				(*(ptr_out[i])).flush();
			}
			return *this;
		}
		type_this & printf(char const * fmt,...){
			va_list args;
			va_start(args,fmt);
			vsprintf(buff,fmt,args);
			va_end(args);
			(*this)<<buff;
			return *this;
		}
	protected:
		//member value:
		Liuze::Type_ArrayTemp<std::ostream*> ptr_out;
		char buff[size_buff];
}; //class multi_ostream;

template<typename OutType,typename ValType,
	typename=typename std::enable_if<std::is_base_of<std::ostream,OutType>::value>::type>
OutType & operator <<(OutType & out,Liuze::Type_ArrayTemp<ValType> const & arr){
	auto iter=arr.begin();
	auto iter_end=arr.end();
	if(iter==iter_end) return out;
	out<<*iter;
	for(++iter;iter!=iter_end;++iter) out<<" "<<*iter;
	return out;
}
//}(ostream matters)
/****************************************************************************************************/

#endif //#ifndef LZ_DEF_hpp_temp220726_utility