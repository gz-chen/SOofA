#if LZ_DEF_LZ_H_basic_src_basic_utility_time_recorder!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_basic_src_utility_time_recorder_CPP
#define LZ_DEF_LZ_H_basic_src_utility_time_recorder_CPP 202105L

namespace Liuze{

/****************************************************************************************************/
//A multi-way time recorder{
//class Time_recorder{
//public:
	
	//static functions:
	
	typename Time_recorder::type_str
	Time_recorder::ns_to_str
	(type_time time_ns){
		std::stringstream res;
		res<<time_ns/31104000000000000LL<<"/"; //year=12month;
		time_ns%=31104000000000000LL;
		res<<time_ns/2592000000000000LL<<"/"; //month=30day;
		time_ns%=2592000000000000LL;
		res<<time_ns/86400000000000LL<<"."; //day=24h;
		time_ns%=86400000000000LL;
		res<<time_ns/3600000000000LL<<":"; //h=60min;
		time_ns%=3600000000000LL;
		res<<time_ns/60000000000LL<<":"; //min=60s;
		time_ns%=60000000000LL;
		res<<time_ns/1000000000LL<<"."; //s=1000ms;
		time_ns%=1000000000LL;
		res<<time_ns/1000000<<"."; //ms=1000us;
		time_ns%=1000000;
		res<<time_ns/1000<<"."; //us=1000ns;
		time_ns%=1000;
		res<<time_ns; //ns;
		return res.str();
	}
	
	//Constructor:
	
	Time_recorder::Time_recorder
	():
	num(1),
	time_bg(1,std::chrono::steady_clock::now().time_since_epoch().count()),
	time_ed(1),
	time_last(1,type_time(0)),time_total(1,type_time(0)),
	func_ns_to_str(&type_this::ns_to_str){
	}
	
	Time_recorder::Time_recorder
	(type_size const & number):
	num(number),
	time_bg(number,std::chrono::steady_clock::now().time_since_epoch().count()),
	time_ed(number),
	time_last(number,type_time(0)),time_total(number,type_time(0)),
	func_ns_to_str(&type_this::ns_to_str){
	}
	
	//other functions:
	
	inline
	typename Time_recorder::type_time
	Time_recorder::last
	(type_size const & id) const {
		return id<num ? time_last[id] : (type_time)0;
	}
	
	inline
	typename Time_recorder::type_time
	Time_recorder::total
	(type_size const & id) const {
		return id<num ? time_total[id] : (type_time)0;
	}
	
	inline
	typename Time_recorder::type_str
	Time_recorder::last_str
	(type_size const & id) const {
		return func_ns_to_str(last(id));
	}
	
	inline
	typename Time_recorder::type_str
	Time_recorder::total_str
	(type_size const & id) const {
		return func_ns_to_str(total(id));
	}
	
	inline
	typename Time_recorder::type_time
	Time_recorder::start
	(type_size const & id){
		return id<num ?
			(time_bg[id]=std::chrono::steady_clock::now().time_since_epoch().count()) :
			type_time(0);
	}
	
	inline
	typename Time_recorder::type_time
	Time_recorder::pause
	(type_size const & id){
		if(id>=num) return type_time(0);
		time_ed[id]=std::chrono::steady_clock::now().time_since_epoch().count();
		time_last[id]=time_ed[id]-time_bg[id];
		time_total[id]+=time_last[id];
		return time_ed[id];
	}
	
	inline
	typename Time_recorder::type_time
	Time_recorder::paustart
	(type_size const & id){
		if(id>=num) return type_time(0);
		time_ed[id]=std::chrono::steady_clock::now().time_since_epoch().count();
		time_last[id]=time_ed[id]-time_bg[id];
		time_total[id]+=time_last[id];
		time_bg[id]=time_ed[id];
		return time_ed[id];
	}
	
	inline
	void
	Time_recorder::reset
	(){
		time_last.assign(num,type_time(0));
		time_total.assign(num,type_time(0));
		time_bg.assign(num,std::chrono::steady_clock::now().time_since_epoch().count());
	}
	
	inline
	void
	Time_recorder::reset
	(type_size const & id){
		if(id<num){
			time_last[id]=time_total[id]=type_time(0);
			time_bg[id]=std::chrono::steady_clock::now().time_since_epoch().count();
		} else if(id==num){
			this->reset();
		}
	}
	
	inline
	typename Time_recorder::type_size
	Time_recorder::size
	() const {
		return num;
	}
	
	template<typename FuncType,typename>
	void
	Time_recorder::Set_ns_to_str
	(FuncType && fun_ns_to_str){
		func_ns_to_str=fun_ns_to_str;
	}
	
//}(class Time_recorder)
//}(A multi-way time recorder)
/****************************************************************************************************/

} //namespace Liuze;

#endif //#ifndef LZ_DEF_LZ_H_basic_src_utility_time_recorder_CPP
#endif //#if LZ_DEF_LZ_H_basic_src_basic_utility_time_recorder!=202105L