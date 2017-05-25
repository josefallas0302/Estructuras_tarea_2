
#ifndef _CLUSTERER_CLASS_
#define _CLUSTERER_CLASS_

#include <iostream>
#include <vector>
#include "cluster.hpp"

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
Clusterer.
 */

//TODO: Improve data_ptr
//TODO: Declare virtual log_likelihood 

/*---------------------------------------------------------------------------*/
/** @brief This class contains some basic functionality and data structures 
required by a set of clustering algorithms (K-Means, GMM-EM, etc). 

    The class manages the clustering data with a pointer to an std vector of 
    Eigen Library n-dimensional vectors (dataset). It also includes an Eigen 
    vector that corresponds to the output of hard-assign clustering algorithms.
    The functionality included is related to reassign the clustering data, 
    printing the dataset and cluster output, and a virtual declaration of the 
    cluster function, which executes the clustering algorithm in the 
    derived class. This class is templated with typename T corresponding to 
    the data type, and the integer N to each datapoint dimension.
*/

//*****************************************************************************
//Clusterer class declaration
//*****************************************************************************
template <typename T, int N> 
class Clusterer
{
public:  

/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class Clusterer. It receives a reference
    to the dataset, and uses it to initialize the data pointer.

    @param[in] n_dataset const std::vector<RowVector<T,N>>&. Reference to the 
    clustering dataset.
*/
   Clusterer(const std::vector<RowVector<T,N>>& n_dataset); 

/*---------------------------------------------------------------------------*/
/** @brief Virtual function corresponding to the clustering algorithm implemented
    in the derived classes.

    The output is a reference to an Eigen Library vector (cluster_out), 
    in which a given position corresponds to the index of the cluster that the 
    datapoint (same position in the dataset vector) was assigned to.
    
    @param[in] debug bool. Flag that when set to true indicates the cluster function 
    to print debugging information.
    @returns RowVectorXi& Reference to cluster_out vector.
*/
   virtual RowVectorXi& cluster(bool debug) = 0;

/*---------------------------------------------------------------------------*/
/** @brief Reassign the data pointer to a new dataset
    
    This function receives a reference to a new clustering dataset and 
    reconfigures the data pointer with this new data.

    @param[in] n_dataset const std::vector<RowVector<T,N>>&. Reference to the 
    clustering dataset.
*/
   void reassign_data(const std::vector<RowVector<T,N>>& n_dataset);

/*---------------------------------------------------------------------------*/
/** @brief Print the dataset points information for debugging purposes.
 */
   void print_data();

/*---------------------------------------------------------------------------*/
/** @brief Print the resulting datapoint/cluster hard assignments.

    @note It is expected to run this function after executing the corresponding
    derived class clustering algorithm (cluster function). Otherwise, the presented 
    values will be the default -1 for every position.
 */
   void print_cluster_out();
  
protected:

/*---------------------------------------------------------------------------*/
/** @brief Pointer to dataset 
 */
   const std::vector<RowVector<T,N>>* data_ptr;

/*---------------------------------------------------------------------------*/
/** @brief Dynamic size Eigen RowVector that contains the cluster index assignment
    for every point in the dataset (hard clustering output)
 */
   RowVectorXi cluster_out;
};

//*****************************************************************************
//Clusterer class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
Clusterer<T,N>::Clusterer(const std::vector<RowVector<T,N>>& n_dataset)
   : data_ptr(&(n_dataset)),
     cluster_out(RowVectorXi::Constant(n_dataset.size(), -1))
{}

template <typename T, int N>
void Clusterer<T,N>::reassign_data(const std::vector<RowVector<T,N>>& n_dataset)
{ 
   this->data_ptr = &(n_dataset);
   this->cluster_out.setZero(n_dataset.size());
}

template <typename T, int N>
void Clusterer<T,N>::print_data()
{
   int precision = 4;
   IOFormat fmt(precision, 0, ", ", "\n", "[", "]");
   for(auto const& it : *this->data_ptr) std::cout << it.format(fmt) << std::endl;
}
    
template <typename T, int N>
void Clusterer<T,N>::print_cluster_out()
{
   int precision = 4;
   IOFormat fmt(precision, 0, " ", "\n", "[", "]");
   std::cout << "Cluster out:\n" << std::endl;
   std::cout << cluster_out.format(fmt) << std::endl;    
}

//*****************************************************************************
//Distance functions namespace
//*****************************************************************************

/*---------------------------------------------------------------------------*/
/** @brief Namespace that includes distance functions (euclidean, etc) for
    several clustering algorithms.
*/
namespace distance
{
/*---------------------------------------------------------------------------*/
/** @brief Euclidean distance between two given points.
    @param[in] p1 const RowVector<T, N>&. First N-Dimensional Point.
    @param[in] p2 const RowVector<T, N>&. Second N-Dimensional Point.
    @returns double Euclidean distance value.
*/
   template<typename T, int N>
   double euclidean(const RowVector<T, N>& p1, const RowVector<T, N>& p2)
   {
      double eucl_dist=(p1-p2).norm(); 
      return eucl_dist;
   }
}


#endif
