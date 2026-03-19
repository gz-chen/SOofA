#if LZ_DEF_LZ_H_math_prime!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_math_src_prime_CPP
#define LZ_DEF_LZ_H_math_src_prime_CPP 202105L

namespace Liuze{
namespace math{

/****************************************************************************************************/

//fun: Is_prime:
template<typename IntType,typename>
Type_Bool Is_prime(IntType number){
	typedef IntType type_int;
	if(number<type_int(0)) number=-number;
	if(number<=type_int(1)) return false;
	std::vector<bool> tag((Type_Size)number+(Type_Size)1,true); //true for prime and false for composite;
	type_int i,p=2; //`p' is a prime;
	while(true){
		i=p*p; //the next remained composite number;
		if(i>=type_int(tag.size())) return true;
		while(i<type_int(tag.size())){
			if(i==number) return false;
			tag[i]=false;
			i+=p;
		}
		do{
			++p;
			if(p==type_int(tag.size())) return true;
		}while(!(tag[p])); //to find the next prime number;
	}
}
//(fun: Is_prime)

//fun: Get_prime_all:
template<typename IntType,typename>
Type_ArrayTemp<IntType> Get_prime_all
(IntType const & upper){
	typedef IntType type_int;
	if(upper<=(type_int)1) return Type_ArrayTemp<type_int>(0);
	std::vector<bool> tag((Type_Size)upper+(Type_Size)1,true); //true for prime and false for composite;
	type_int i,p=2;
	Type_Size n=0; //number of prime numbers;
	while(true){
		++n; //p is a prime number;
		i=p*p; //the next remained composite number;
		if(i>=type_int(tag.size())) break;
		while(i<type_int(tag.size())){
			tag[i]=false;
			i+=p;
		}
		do{
			++p;
			if(p==type_int(tag.size())) goto GTS_end;
		}while(!(tag[p])); //to find the next prime number;
	}
	for(++p;p<type_int(tag.size());++p){
		if(tag[p]) ++n;
	}
	GTS_end:
	Type_ArrayTemp<type_int> res(n); //the list of prime numbers;
	i=0;
	for(p=2;i<(type_int)n;++p){
		if(tag[p]){
			res[i]=p;
			++i;
		}
	}
	return res;
}
//(fun: Get_prime_all)

//fun: Get_prime_all:
template<typename IntType,typename>
Type_ArrayTemp<IntType> Get_prime_all
(IntType const & upper,Type_ArrayTemp<IntType> const & prior_prime,IntType prior_upper){
	typedef IntType type_int;
	type_int n_prior=prior_prime.size(),n_post=0; //number of prime numbers;
	if(n_prior==(type_int)0) prior_upper=(type_int)1;
	else if(prior_upper<prior_prime[n_prior-1]) prior_upper=prior_prime[n_prior-1];
	if(upper<=prior_upper) return prior_prime;
	type_int i,i_p,p,i0=prior_upper+(type_int)1;
	//true for prime and false for composite, standing for {upper_prior+1,...,upper}:
	std::vector<bool> tag(upper-prior_upper,true);
	//screening using the prime numbers found before:
	for(i_p=0;i_p<n_prior;++i_p){
		p=prior_prime[i_p];
		i=(prior_upper/p)*p;
		while((i+=p)<=upper) tag[i-i0]=false;
	}
	//screening using the prime numbers newly found:
	p=i0;
	while(p<=upper){
		while(!(tag[p-i0])){
			++p;
			if(p>upper) goto GTS_end;
		}
		++n_post; //p is a new prime number;
		i=p;
		while((i+=p)<=upper) tag[i-i0]=false;
		++p;
	}
	//summarising the result:
	GTS_end:
	if(n_post==(type_int)0) return prior_prime;
	Type_ArrayTemp<type_int> res(n_prior+n_post);
	for(i_p=0;i_p<n_prior;++i_p) res[i_p]=prior_prime[i_p];
	for(i=0;i_p<res.size();++i){
		if(tag[i]){
			res[i_p]=i+i0;
			++i_p;
		}
	}
	return res;
}
//(fun: Get_prime_all)

//fun: Get_prime_atleast:
template<typename ForwardIterType,typename OutputIterType,typename IntType,typename>
Type_ArrayTemp<IntType> Get_prime_atleast
(ForwardIterType const & least_begin,ForwardIterType const & least_end,OutputIterType prime_begin,
Type_ArrayTemp<IntType> const prime_list,IntType prime_list_upper){
	typedef IntType type_int;
	if(least_begin==least_end) return Type_ArrayTemp<type_int>((Type_Size)0);
	type_int least_max=*(std::max_element(least_begin,least_end));
	Type_ArrayTemp<type_int> prime_list_full=Get_prime_all(least_max,prime_list,prime_list_upper);
	//least_max is composite:
	if(prime_list_full[prime_list_full.size()-1]!=least_max){
		type_int p=least_max,rp,i;
		Type_Bool tag=false;
		while(!tag){
			++p;
			rp=type_int(ceil(sqrt((Type_Real)p)));
			tag=true;
			for(i=0;tag && i<prime_list_full.size() && prime_list_full[i]<=rp;++i){
				if(p%prime_list_full[i]==(type_int)0) tag=false;
			}
		}
		prime_list_full.push_back(p); //p is a prime number;
	}
	Type_Size iL,iR,iM;
	for(ForwardIterType least=least_begin;least!=least_end;++least){
		iL=0;
		iR=prime_list_full.size()-(Type_Size)1;
		while(iL<iR){
			iM=(iL+iR)/(Type_Size)2;
			if((*least)>prime_list_full[iM]) iL=iM+1;
			else iR=iM;
		}
		(*prime_begin)=prime_list_full[iL];
		++prime_begin;
	}
	return prime_list_full;
}
//(fun: Get_prime_atleast)

//fun: integer_prime_decomposition:
template<typename IntType,typename>
Type_ArrayTemp<Type_ArrayTemp<IntType> > integer_prime_decomposition
(IntType const & number){
	typedef IntType type_int;
	typedef Type_ArrayTemp<type_int> type_arr_int;
	typedef Type_ArrayTemp<type_arr_int> type_res;
	type_int num= number>=type_int(0) ? number : -number;
	if(num<=type_int(1)) return type_res(0);
	type_res decomp(0);
	type_int prim=2,i;
	std::vector<bool> flag(num+type_int(1),true); //true for prime and false for composite;
	while(true){
		//`prim' is a prime number.
		if(num%prim==type_int(0)){
			i=decomp.size();
			decomp.push_back(type_arr_int({prim,type_int(1)}));
			while((num/=prim)%prim==type_int(0)) ++decomp[i][1];
			if(num==type_int(1)) return decomp;
		}
		i=prim*prim; //the next remained composite number;
		while(i<=num){
			flag[i]=false;
			i+=prim;
		}
		do{
			++prim;
		}while(!flag[prim]); //to find the next prime number;
	}
}
//(fun: integer_prime_decomposition)

/****************************************************************************************************/
//class: PrimePowerQuerier{
//public:
	
	//constructor:
	
	template<typename IntType,typename TraitsCheck>
	PrimePowerQuerier<IntType,TraitsCheck>::PrimePowerQuerier
	(type_int upper):
	m_upper(upper>=type_int(2) ? upper : type_int(0)),m_tab_asc(),m_tab_bas(){
		this->Cal_init();
	}
	
	template<typename IntType,typename TraitsCheck>
	PrimePowerQuerier<IntType,TraitsCheck>::PrimePowerQuerier
	(type_this const & querier):
	m_upper(querier.m_upper),m_tab_asc(querier.m_tab_asc),m_tab_bas(querier.m_tab_bas){
	}
	
	template<typename IntType,typename TraitsCheck>
	PrimePowerQuerier<IntType,TraitsCheck>::PrimePowerQuerier
	(type_this && querier):
	m_upper(std::move(querier.m_upper)),
	m_tab_asc(std::move(querier.m_tab_asc)),m_tab_bas(std::move(querier.m_tab_bas)){
	}
	
	//operator:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_this &
	PrimePowerQuerier<IntType,TraitsCheck>::operator =
	(type_this const & querier){
		if(std::addressof(querier)==this) return *this;
		m_upper=querier.m_upper;
		m_tab_asc=querier.m_tab_asc;
		m_tab_bas=querier.m_tab_bas;
		return *this;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_this &
	PrimePowerQuerier<IntType,TraitsCheck>::operator =
	(type_this && querier){
		m_upper=std::move(querier.m_upper);
		m_tab_asc=std::move(querier.m_tab_asc);
		m_tab_bas=std::move(querier.m_tab_bas);
		return *this;
	}
	
	//other function:
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_int const &
	PrimePowerQuerier<IntType,TraitsCheck>::get_upper
	() const {
		return m_upper;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_tab const &
	PrimePowerQuerier<IntType,TraitsCheck>::get_table_ascending
	() const {
		return m_tab_asc;
	}
	
	template<typename IntType,typename TraitsCheck>
	inline
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_tab const &
	PrimePowerQuerier<IntType,TraitsCheck>::get_table_base
	() const {
		return m_tab_bas;
	}
	
	template<typename IntType,typename TraitsCheck>
	typename PrimePowerQuerier<IntType,TraitsCheck>::type_int
	PrimePowerQuerier<IntType,TraitsCheck>::set_upper
	(type_int const & upper){
		if(upper<=type_int(1) || upper<=m_upper) return m_upper;
		Type_Size size_m_tab_asc_old=m_tab_asc.size();
		if(m_upper==type_int(0)) m_upper=type_int(1);
		type_int i_bas,p_bas,p_exp,n,n0=m_upper+type_int(1);
		//true for prime and false for composite, standing for {`m_upper'+1,...,`upper'}:
		std::vector<bool> flag(upper-m_upper,true);
		//screening using the prime numbers found before{
		for(i_bas=0;i_bas<type_int(m_tab_bas.size());++i_bas){
			p_exp=m_tab_bas[i_bas].size();
			p_bas=m_tab_bas[i_bas][0];
			//finding powers of previous primes{
			n=m_tab_bas[i_bas][p_exp-type_int(1)]*p_bas;
			while(n<=upper){
				++p_exp;
				m_tab_asc.push_back(type_subtab({n,m_tab_bas[i_bas][0],p_exp}));
				m_tab_bas[i_bas].push_back(n);
				n*=m_tab_bas[i_bas][0];
			}
			//}
			//finding multiples of previous primes{
			n=(m_upper/p_bas)*p_bas;
			while((n+=p_bas)<=upper) flag[n-n0]=false;
			//}
		}
		//}
		//screening using the prime numbers newly found{
		p_bas=m_upper;
		while(true){
			//finding the next prime number{
			do{
				++p_bas;
				if(p_bas>upper) goto GTL_end;
			}while(!flag[p_bas-n0]);
			//}
			//storing the prime number `p_bas'{
			p_exp=type_int(1);
			m_tab_asc.push_back(type_subtab({p_bas,p_bas,p_exp}));
			m_tab_bas.push_back(type_subtab(1,p_bas));
			//}
			//finding powers of `p_bas'{
			n=p_bas*p_bas;
			while(n<=upper){
				++p_exp;
				m_tab_asc.push_back(type_subtab({n,p_bas,p_exp}));
				m_tab_bas[m_tab_bas.size()-Type_Size(1)].push_back(n);
				n*=p_bas;
			}
			//}
			//finding multiples of `p_bas'{
			if(p_exp==type_int(1)) break;
			n=p_bas*p_bas; //the next remained composite number;
			do{
				flag[n-n0]=false;
				n+=p_bas;
			}while(n<=upper);
			//}
		}
		p_exp=type_int(1);
		for(++p_bas;p_bas<=upper;++p_bas){
			if(flag[p_bas-n0]){
				m_tab_asc.push_back(type_subtab({p_bas,p_bas,p_exp}));
				m_tab_bas.push_back(type_subtab(1,p_bas));
			}
		}
		//}
		GTL_end:
		std::sort(m_tab_asc.begin()+size_m_tab_asc_old,m_tab_asc.end(),
			[](type_subtab const & tab0,type_subtab const & tab1)->bool{return tab0[0]<tab1[0];});
		m_upper=upper;
		return m_upper;
	}
	
//protected:
	
	//other function:
	
	template<typename IntType,typename TraitsCheck>
	void
	PrimePowerQuerier<IntType,TraitsCheck>::Cal_init
	(){
		if(m_upper==type_int(0)) return;
		m_tab_asc.resize(0);
		m_tab_bas.resize(0);
		std::vector<bool> flag(m_upper+type_int(1),true); //true for prime and false for composite;
		type_int p_exp,n,p_bas=2;
		while(true){
			//storing the prime number `p_bas'{
			p_exp=type_int(1);
			m_tab_asc.push_back(type_subtab({p_bas,p_bas,p_exp}));
			m_tab_bas.push_back(type_subtab(1,p_bas));
			//}
			//finding powers of `p_bas'{
			n=p_bas*p_bas;
			while(n<=m_upper){
				++p_exp;
				m_tab_asc.push_back(type_subtab({n,p_bas,p_exp}));
				m_tab_bas[m_tab_bas.size()-Type_Size(1)].push_back(n);
				n*=p_bas;
			}
			//}
			//finding multiples of `p_bas'{
			if(p_exp==type_int(1)) break;
			n=p_bas*p_bas; //the next remained composite number;
			do{
				flag[n]=false;
				n+=p_bas;
			}while(n<=m_upper);
			//}
			//finding the next prime number{
			do{
				++p_bas;
				if(p_bas>m_upper) goto GTL_end;
			}while(!flag[p_bas]);
			//}
		}
		p_exp=type_int(1);
		for(++p_bas;p_bas<=m_upper;++p_bas){
			if(flag[p_bas]){
				m_tab_asc.push_back(type_subtab({p_bas,p_bas,p_exp}));
				m_tab_bas.push_back(type_subtab(1,p_bas));
			}
		}
		GTL_end:
		std::sort(m_tab_asc.begin(),m_tab_asc.end(),
			[](type_subtab const & tab0,type_subtab const & tab1)->bool{return tab0[0]<tab1[0];});
		return;
	}
	
//}(class: PrimePowerQuerier)
/****************************************************************************************************/

} //namespace math;
} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_math_src_prime_CPP
#endif //#if LZ_DEF_LZ_H_math_prime!=202105L