
#ifndef _CLUSTER_CLASS_
#define _CLUSTER_CLASS_

#include <Eigen/Dense>

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
 Cluster.
*/

using namespace Eigen;

//Row Vector typedefs, which are not defined in the default Eigen Library Typedefs

/** @brief Fixed size (N) Row Vector typedef (not defined in the default 
    Eigen Library Typedefs)
*/
template <typename T, int N>
using RowVector = Matrix<T,1,N>;
/** @brief Dynamic size Row Vector typedef (not defined in the default 
    Eigen Library Typedefs)
*/
template <typename T>
using RowVectorX = Matrix<T,1,Dynamic>;


/*---------------------------------------------------------------------------*/
/** @brief This class describes the basic structure of the Cluster object.
    These objects contain important information for each group (cluster) fitted
    by the Clustering algorithms (K-Means, GMM-EM, etc).
*/

//*****************************************************************************
//Cluster class declaration
//*****************************************************************************
template <typename T, int N>
class Cluster
{
public:
/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class Cluster. It initializes the 
    cluster centroid to the zero vector.
*/
   Cluster();

/*---------------------------------------------------------------------------*/
/** @brief Return a non-const reference to the Cluster centroid (read/write)
    @returns RowVector<T,N>& Reference to centroid vector.
*/
   
   RowVector<T,N>& centroid();
  
protected:
/*---------------------------------------------------------------------------*/
/** @brief Eigen RowVector corresponding to the Cluster centroid position.
*/
   RowVector<T,N> centroid_;
};

//*****************************************************************************
//Cluster class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
Cluster<T,N>::Cluster()
   : centroid_(RowVector<T,N>::Zero())
{}

template <typename T, int N>
RowVector<T,N>& Cluster<T,N>::centroid()
{
   return this->centroid_;
}

#endif
