#if LZ_DEF_LZ_H_math_differentiation!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_diferentiation_CPP
#define LZ_DEF_LZ_H_math_src_diferentiation_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/
//class Derivative{
//public:
	
	//type: type_pre{
	//public:
		
		template<typename RealType>
		Derivative<RealType>::type_pre::type_pre
		(type_real length_initial,type_real error_max,type_size iteration_max,type_stat error_type):
		length_init(length_initial>type_real(0) ? length_initial : type_real(1)),
		err(error_max>type_real(0) ? error_max : type_real(LZ_DEF_const_default_precision)),
		n_iter_max((iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
		iteration_max : type_size(100)),
		err_type(error_type){
		}
		
	//}(type_pre)
	
//public:
	//Constructor:
	
	template<typename RealType>
	Derivative<RealType>::Derivative
	():
	pt_func(NULL),in_nrow(0),in_ncol(0),out_nrow(0),out_ncol(0),der_nrow(0),der_ncol(0),
	stat(type_this::val_stat_null),pre(){
	}
	
	template<typename RealType>
	template<typename FuncType,typename>
	Derivative<RealType>::Derivative
	(FuncType && func,type_size nrow_out,type_size ncol_out,type_size nrow_in,type_size ncol_in):
	pt_func(new type_func(std::forward<FuncType>(func))),pre(){
		if(pt_func){
			in_nrow= nrow_in>(type_size)0 ? nrow_in : (type_size)0;
			in_ncol= ncol_in>(type_size)0 ? ncol_in : (type_size)0;
			out_nrow= nrow_out>(type_size)0 ? nrow_out : (type_size)0;
			out_ncol= ncol_out>(type_size)0 ? ncol_out : (type_size)0;
			if(in_nrow==(type_size)0 || in_ncol==(type_size)0 ||
			out_nrow==(type_size)0 || out_ncol==(type_size)0){
				stat=type_this::val_stat_incomplete;
				der_nrow=der_ncol=(type_size)0;
			} else {
				stat=type_this::val_stat_complete;
				der_nrow=in_nrow*out_nrow;
				der_ncol=in_ncol*out_ncol;
			}
		} else {
			in_nrow=in_ncol=out_nrow=out_ncol=der_nrow=der_ncol=(type_size)0;
			stat=type_this::val_stat_null;
		}
	}
	
	//Destructor:
	
	template<typename RealType>
	Derivative<RealType>::~Derivative
	(){
		if(pt_func) delete pt_func;
	}
	
	//operator:
	
	//The following operator () calculates the derivative numerically.
	//The following comment introduces the calculating method.
	//Suppose a primitive function f: R -> R has 5-order continuous derivative,
	//then by the Talor formula at a real number x, for h>0,
	//f(x+h) = f(x) + f'(x)h + f''(x)/2*h^2 + f'''(x)/6*h^3 + f''''(x)/24*h^4 + O(h^5),
	//f(x-h) = f(x) - f'(x)h + f''(x)/2*h^2 - f'''(x)/6*h^3 + f''''(x)/24*h^4 + O(h^5),
	//so d(f,x,h) := (f(x+h)-f(x-h))/(2h) = f'(x) + f'''(x)/6*(h^2) + O(h^4),
	//and similarly, d(f,x,h/2) = f'(x) + f'''(x)/6*(h^2)/4 + O(h^4),
	//thus D(f,x,h/2) := (4d(f,x,h/2)-d(f,x,h))/3 = f'(x) + O(h^4).
	//The following algorithm bases on th above method.
	//If the primitive has (2K+1)-order continuous derivative,
	//we can use K iterations with concentrateing rate r to compute
	//d(f,x,h/r^0), ..., d(f,x,h/r^(K-1)),
	//and combine them using K coefficients
	//who are the entries of the first column of the inverse of the Vandermonde matrix
	//whose parameters are r^0, r^(-2), ..., r^(-2(K-1)),
	//i.e., the combination coefficient of d(f,x,h/r^i) is
	//prod{r^(-2j)/(r^(-2j)-r^(-2i)) : j in {0,...,K-1}\{i}};
	//doing so will improve the precision to O(h^(2K));
	//the following algorithm does not do this.
	template<typename RealType>
	typename Derivative<RealType>::type_matrix
	Derivative<RealType>::operator ()
	(type_matrix const & x){
		if(this->Cal_dim(x)==false) return type_matrix();
		type_matrix res(der_nrow,der_ncol);
		type_real const crr2=0.5;
		type_size ixr,ixc,iyr,iyc; //row and column indices for input, output;
		type_real tr0,tr1; //temp real variable;
		type_matrix yR,yL; //the left and right values of output;
		type_real h0,h1; //the step lengths;
		type_matrix d0,d1,D0,D1; //the calculated derivatives with respect to a 1 by 1 input component;
		type_matrix x1=x;
		type_size i_iter; //the index for iteration;
		for(ixr=0;ixr<in_nrow;++ixr){
			for(ixc=0;ixc<in_ncol;++ixc){
				//this block calculates (partial primitive)/(partial x(ixr,ixc)).
				tr0=Liuze::math::abs(x(ixr,ixc));
				h0=pre.length_init;
				h1=h0*crr2;
				x1(ixr,ixc)=x(ixr,ixc)+h1;
				yR=(*pt_func)(x1);
				x1(ixr,ixc)=x(ixr,ixc)-h1;
				yL=(*pt_func)(x1);
				D0=d0=(yR-yL)/h0;
				i_iter=0;
				while(pre.n_iter_max==type_pre::val_size_infinity || i_iter<pre.n_iter_max){
					h0=h1; //step length is multiplied by 0.5;
					h1*=crr2;
					x1(ixr,ixc)=x(ixr,ixc)+h1;
					yR=(*pt_func)(x1);
					x1(ixr,ixc)=x(ixr,ixc)-h1;
					yL=(*pt_func)(x1);
					d1=(yR-yL)/h0; //the new low-precision derivative;
					D1=(type_real(4)*d1-d0)/type_real(3); //the new high-precision derivative;
					++i_iter;
					tr1=tr0+h1;
					for(iyr=0;iyr<out_nrow;++iyr){
						for(iyc=0;iyc<out_ncol;++iyc){
							if(this->Cal_err(tr0,tr1,D0(iyr,iyc),D1(iyr,iyc))==false){
								goto GTS_err_large;
							}
						}
					}
					break;
					GTS_err_large:
					d0=d1;
					D0=D1;
				}
				//write (partial primitive)/(partial x(ixr,ixc)) to the result:
				for(iyr=0;iyr<out_nrow;++iyr){
					for(iyc=0;iyc<out_ncol;++iyc){
						res(iyr*in_nrow+ixr,iyc*in_ncol+ixc)=D1(iyr,iyc);
					}
				}
				x1(ixr,ixc)=x(ixr,ixc); //re-cover x1=x;
			}
		}
		return res;
	}
	
	//Get:
	
	template<typename RealType>
	inline
	typename Derivative<RealType>::type_matrix
	Derivative<RealType>::Get_function
	(type_matrix const & x){
		return stat!=type_this::val_stat_null ? (*pt_func)(x) : type_matrix();
	}
	
	//Set:
	
	template<typename RealType>
	inline
	void
	Derivative<RealType>::Set_function
	(){
		if(pt_func){
			delete pt_func;
			pt_func=NULL;
		}
		in_nrow=in_ncol=out_nrow=out_ncol=der_nrow=der_ncol=(type_size)0;
		stat=type_this::val_stat_null;
	}
	
	template<typename RealType>
	template<typename FuncType,typename>
	void
	Derivative<RealType>::Set_function
	(FuncType && func,type_size nrow_out,type_size ncol_out,type_size nrow_in,type_size ncol_in){
		if(pt_func){
			delete pt_func;
			pt_func=NULL;
		}
		pt_func=new type_func(std::forward<FuncType>(func));
		if(pt_func){
			in_nrow= nrow_in>(type_size)0 ? nrow_in : (type_size)0;
			in_ncol= ncol_in>(type_size)0 ? ncol_in : (type_size)0;
			out_nrow= nrow_out>(type_size)0 ? nrow_out : (type_size)0;
			out_ncol= ncol_out>(type_size)0 ? ncol_out : (type_size)0;
			if(in_nrow==(type_size)0 || in_ncol==(type_size)0 ||
			out_nrow==(type_size)0 || out_ncol==(type_size)0){
				stat=type_this::val_stat_incomplete;
				der_nrow=der_ncol=(type_size)0;
			} else {
				stat=type_this::val_stat_complete;
				der_nrow=in_nrow*out_nrow;
				der_ncol=in_ncol*out_ncol;
			}
		} else {
			in_nrow=in_ncol=out_nrow=out_ncol=der_nrow=der_ncol=(type_size)0;
			stat=type_this::val_stat_null;
		}
	}
	
	template<typename RealType>
	inline
	void
	Derivative<RealType>::Set_precision
	(type_real length_initial,type_real error_max,type_size iteration_max,type_stat error_type){
		pre.length_init= length_initial>type_real(0) ? length_initial : type_real(1);
		pre.err= error_max>type_real(0) ? error_max : type_real(LZ_DEF_const_default_precision);
		pre.n_iter_max= (iteration_max>type_size(0) || iteration_max==type_pre::val_size_infinity) ?
			iteration_max : type_size(100);
		pre.err_type=error_type;
	}
	
//protected:
	
	//Cal:
	
	template<typename RealType>
	Type_Bool
	Derivative<RealType>::Cal_dim
	(type_matrix const & x){
		if(stat==type_this::val_stat_null || x.rows()==(type_size)0 || x.cols()==(type_size)0){
			return false;
		}
		if(stat==type_this::val_stat_complete){
			return (x.rows()==in_nrow && x.cols()==in_ncol);
		}
		//the following deals with the case stat==type_this::val_stat_incomplete.
		type_size inr,inc,onr,onc; //in_nrow,in_ncol,out_nrow,out_ncol;
		if(in_nrow>(type_size)0){
			if(x.rows()!=in_nrow) return false;
			inr=in_nrow;
		} else {
			inr=x.rows();
		}
		if(in_ncol>(type_size)0){
			if(x.cols()!=in_ncol) return false;
			inc=in_ncol;
		} else {
			inc=x.cols();
		}
		type_matrix y=(*pt_func)(x);
		if(y.rows()==(type_size)0 || y.cols()==(type_size)0) return false;
		if(out_nrow>(type_size)0){
			if(y.rows()!=out_nrow) return false;
			onr=out_nrow;
		} else {
			onr=y.rows();
		}
		if(out_ncol>(type_size)0){
			if(y.cols()!=out_ncol) return false;
			onc=out_ncol;
		} else {
			onc=y.cols();
		}
		in_nrow=inr;
		in_ncol=inc;
		out_nrow=onr;
		out_ncol=onc;
		der_nrow=in_nrow*out_nrow;
		der_ncol=in_ncol*out_ncol;
		stat=type_this::val_stat_complete;
		return true;
	}
	
	template<typename RealType>
	Type_Bool
	Derivative<RealType>::Cal_err
	(type_real const & x0,type_real const & x1,type_real const & y0,type_real const & y1) const {
		switch(pre.err_type){
			case type_pre::val_err_in_abs:
				return Liuze::math::abs(x0-x1)<=pre.err;
			case type_pre::val_err_in_rel:
				if(x1!=(type_real)0) return Liuze::math::abs((x0-x1)/x1)<=pre.err;
				return x0==(type_real)0 ? true : (type_real)1<=pre.err;
			case type_pre::val_err_out_abs:
				return Liuze::math::abs(y0-y1)<=pre.err;
			case type_pre::val_err_out_rel:
				if(y1!=(type_real)0) return Liuze::math::abs((y0-y1)/y1)<=pre.err;
				return y0==(type_real)0 ? true : (type_real)1<=pre.err;
			default:
				return false;
		}
	}
	
//}(class Derivative)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_diferentiation_CPP
#endif //#if LZ_DEF_LZ_H_math_differentiation!=202105L