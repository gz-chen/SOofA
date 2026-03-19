namespace Liuze{
namespace stats{

template<typename _OutType=_DFT_RealType,typename _InType=Eigen::MatrixX<_DFT_RealType>,
	typename _WeightType=_DFT_RealType>
class Data_Response{
	LZ_DEF_func_check_traits(std::is_floating_point<_WeightType>::value);
	public:
		//Type: this:
		typedef Data_Response<_OutType,_InType,_WeightType> type_this;
		//Type: output:
		typedef _OutType type_out;
		//Type: input:
		typedef _Input type_in;
		//Type: weight:
		typedef _WeightType type_weight;
		//Type: the iterators:
		template<typename _Tp>
		using type_iter=std::iterator<std::input_iterator_tag,_Tp const>;
		
		//Constructor:
		template<typename _ArrayOutType,typename _ArrayInType,typename _ArrayWeightType>
		Data_Response
		(_ArrayOutType const & output,_ArrayInType const & input,_ArrayWeightType const & weight):
		{
		}
	protected:
		//values:
		type_iter<type_out> 
};

}
}