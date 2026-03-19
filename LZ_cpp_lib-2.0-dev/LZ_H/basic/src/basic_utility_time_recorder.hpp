#if LZ_DEF_LZ_H_basic_basic!=202105L
	#error "This file should not be included alone!"
#else

#ifndef LZ_DEF_LZ_H_basic_src_basic_utility_time_recorder
#define LZ_DEF_LZ_H_basic_src_basic_utility_time_recorder 202105L

#include<chrono>
#include<string>
#include<sstream>

namespace Liuze{
/****************************************************************************************************/
//A multi-way time recorder{
class Time_recorder{
	public:
		//Type: this class:
		typedef Time_recorder type_this;
		//Type: the unsigned integral type for time in 10^(-9) second:
		typedef decltype(std::chrono::steady_clock::now().time_since_epoch().count()) type_time;
		typedef type_time value_type;
		//Type: the size:
		typedef Liuze::Type_Size type_size;
		//Type: the sequence of times:
		typedef Liuze::Type_ArrayTemp<type_time> type_time_seq;
		//Type: string:
		typedef std::string type_str;
		
		//static functions:
		static type_str ns_to_str(type_time time_ns);
		
		//Constructor:
		Time_recorder();
		explicit Time_recorder(type_size const & number);
		Time_recorder(type_this const & recorder)=delete;
		
		//Destructor:
		~Time_recorder()=default;
		
		//other functions:
		type_time last(type_size const & id) const;
		type_time total(type_size const & id) const;
		type_str last_str(type_size const & id) const;
		type_str total_str(type_size const & id) const;
		type_time start(type_size const & id);
		type_time pause(type_size const & id);
		type_time paustart(type_size const & id);
		void reset();
		void reset(type_size const & id);
		type_size size() const;
		
		template<typename FuncType,
			typename=typename std::enable_if<
			std::is_convertible<FuncType&&,std::function<type_str(type_time)> >::value>::type>
		void Set_ns_to_str(FuncType && fun_ns_to_str);
	protected:
		type_size num; //the number of recorders;
		type_time_seq time_bg,time_ed;
		type_time_seq time_last,time_total;
		std::function<type_str(type_time)> func_ns_to_str; //a function converting time in ns to a string;
}; //class Time_recorder;
//}(A multi-way time recorder)
/****************************************************************************************************/
} //namespace Liuze;

//implementation:
#include"LZ_H/basic/src/utility/time_recorder.cpp"

#endif //#ifndef LZ_DEF_LZ_H_basic_src_basic_utility_time_recorder
#endif //#if LZ_DEF_LZ_H_basic_basic!=202105L