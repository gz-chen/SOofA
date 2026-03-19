#ifndef LZ_DEF_LZ_H_basic_src_basic_Eigen
#define LZ_DEF_LZ_H_basic_src_basic_Eigen 202105L

#include<Eigen/Core>
/****************************************************************************************************/
//define the Eigen::MatrixX type{
namespace Eigen{
	#ifndef LZ_DEF_type_Eigen_MatrixX
	#define LZ_DEF_type_Eigen_MatrixX 202105L
	
	//typedef: dynamic matrix in Eigen.
	#if LZ_DEF_compile_ignore
	//`#if !EIGEN_VERSION_AT_LEAST(3,4,0)
	template<typename Scalar>
	using ArrayX=Array<Scalar,Dynamic,Dynamic>;
	#endif
	template<typename Scalar>
	using MatrixX=Matrix<Scalar,Dynamic,Dynamic>;
	template<typename Scalar>
	using VectorX=Matrix<Scalar,Dynamic,1>;
	template<typename Scalar>
	using RowVectorX=Matrix<Scalar,1,Dynamic>;
	
	#endif //#ifndef LZ_DEF_type_Eigen_MatrixX
} //namespace Eigen;
//}(define the Eigen::MatrixX type)
/****************************************************************************************************/

#endif //#ifndef LZ_DEF_LZ_H_basic_src_basic_Eigen