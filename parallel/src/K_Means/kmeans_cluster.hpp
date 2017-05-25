
#ifndef _KMEANS_CLUSTER_CLASS_
#define _KMEANS_CLUSTER_CLASS_

#include "../cluster.hpp"

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
 KM_Cluster (K-Means Cluster).
*/

/*---------------------------------------------------------------------------*/
/** @brief This class extends the Cluster base class properties by adding an
    integer size for the cluster object, as required by the K-Means algorithm.
*/

//*****************************************************************************
//KM_Cluster class declaration
//*****************************************************************************
template <typename T, int N>
class KM_Cluster: public Cluster<T,N>
{
public:
/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class KM_Cluster. It calls the Cluster
    base class constructor and initializes the size member to zero.
*/
   KM_Cluster();

/*---------------------------------------------------------------------------*/
/** @brief Return a non-const reference to the Cluster size (read/write)
    @returns int& Reference to the cluster size.
*/
   int& size();
  
private:  
/*---------------------------------------------------------------------------*/
/** @brief Integer corresponding to the Cluster size.
 */
   int size_;
};

//*****************************************************************************
//KM_Cluster class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
KM_Cluster<T,N>::KM_Cluster()
   :Cluster<T,N>(), size_(0)
{}

template <typename T, int N>
int& KM_Cluster<T,N>::size()
{
   return this->size_;
}

#endif
