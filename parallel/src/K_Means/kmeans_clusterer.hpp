
#ifndef _KMEANS_CLUSTERER_CLASS_
#define _KMEANS_CLUSTERER_CLASS_

#include <cstdlib>
#include "../clusterer.hpp"
#include "kmeans_cluster.hpp"
#include <thread>

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
 KM_Clusterer (K-Means Clusterer).
*/

//TODO: Implement configurable initialization method: Random, Kmeans++
//TODO: Implement reconfigure method
//TODO: Distance functor
//TODO: Add max_iterations to algorithm and/or treshold
//TODO: Check LLoyd's algorithm, several runs best fit (inertia)
//TODO: Implement cluster_vector from base class with pointers

/*---------------------------------------------------------------------------*/
/** @brief This class contains the implementation of the K-Means clustering 
    algorithm. It uses the KM_Cluster class to describe
    each cluster object (centroid and size).
*/

//*****************************************************************************
//KM_Clusterer class declaration
//*****************************************************************************
template <typename T, int N>
class KM_Clusterer: public Clusterer<T, N>
{
public:
/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class KM_Clusterer. It receives a 
    reference to the dataset (to initialize the data pointer) and the number of 
    clusters (K) to fit in the algorithm. 

    @param[in] n_dataset const std::vector<RowVector<T,N>>&. Reference to the 
    clustering dataset.
    @param[in] k_size int. Number of clusters K.
*/
   KM_Clusterer(const std::vector<RowVector<T,N>>& n_dataset, int k_size);

/*---------------------------------------------------------------------------*/
/** @brief Execute the K-Means clustering algorithm.

    This function iterates assigning each datapoint the its closest centroid 
    and then recalculates these centroids. The initial centroids locations
    are choosed randomly. The output is the datapoint-cluster assignment in the 
    cluster_out vector, indicating the closest cluster id as the value, and the 
    corresponding position as the datapoint (same position in the dataset vector).
    
    @param[in] print bool. Flag that when set to true indicates the cluster function 
    to print debugging information.
    @returns RowVectorXi& Reference to cluster_out vector.
*/
   RowVectorXi& cluster(bool print=false);

/*---------------------------------------------------------------------------*/
/** @brief Return a constant reference to the cluster vector, in which a certain 
    position contains a cluster fitted by the algorithm (this position corresponds
    to the cluster index). The size of this vector is K.

    @note It is expected to run this function after executing the clustering 
    algorithm (cluster function). Otherwise, the returned vector will contain the
    default clusters (zero vector as centroid, size 0).
    
    @returns const std::vector<KM_Cluster<double,N>>& Reference to the cluster vector.
*/   
   const std::vector<KM_Cluster<double,N>>& get_clusters();

/*---------------------------------------------------------------------------*/
/** @brief Print the information of all clusters in the cluster vector 
    (centroid and size).
*/      
   void print_cluster_info();

   void call_from_thread(int tid);
  
private:
/*---------------------------------------------------------------------------*/
/** @brief Vector of the K clusters fitted by the K-Means algorithm. A certain
    position number in this vector (starting from 0) is the corresponding 
    cluster index. 
*/
   std::vector<KM_Cluster<double,N>> cluster_vector;

/*---------------------------------------------------------------------------*/
/** @brief Number of clusters to be fitted by the K-Means algorithm.
*/
   int K;

/*---------------------------------------------------------------------------*/
/** @brief Assign the initial centroids positions according to an uniform random
    distribution, choosing K points from the Dataset. This is used inside the 
    cluster function (K-Means algorithm).
*/   
   void assign_random_clusters();

/*---------------------------------------------------------------------------*/
/** @brief Assign every datapoint to its closest cluster centroid, modifying the 
    cluster_out vector accordingly. This is used inside the cluster function 
    (K-Means algorithm).
*/ 
   bool assign_clusters();

/*---------------------------------------------------------------------------*/
/** @brief Get the nearest cluster id from a given datapoint. This is used inside 
    the assign_clusters function.

    @param[in] data_point const RowVector<T,N>&. Reference to the datapoint used to 
    find its nearest cluster index.
    @returns int Nearest cluster id.
    
*/  
   int nearest_centroid(const RowVector<T,N>& data_point);

/*---------------------------------------------------------------------------*/
/** @brief Recalculate each cluster centroid as the mean of the data points assigned
    to it (by assign_clusters). This is used inside the cluster function (K-Means 
    algorithm).
*/      
   void recalculate_centroids();  
};

//*****************************************************************************
//KM_Clusterer class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
KM_Clusterer<T,N>::KM_Clusterer(const std::vector<RowVector<T,N>>& n_dataset, int k_size)
   : Clusterer<T,N>(n_dataset),
   cluster_vector(k_size, KM_Cluster<double,N>()),
   K(k_size)
{}

template <typename T, int N>
RowVectorXi& KM_Clusterer<T,N>::cluster(bool debug)
{
   int iteration=1;
   bool continue_flag;
    
   if(debug) std::cout << "\n\n************** KMEANS CLUSTERING FUNCTION **************"
		       << std::endl;

   this->assign_random_clusters();
   this->recalculate_centroids();
    
   if(debug)
   {
      std::cout << "\n================= Initial cluster info =================\n" << std::endl;
      this->print_cluster_info(); std::cout<< std::endl; this->print_cluster_out();
      std::cout << "\n========================================================\n" << std::endl;
      std::cout << "\n================= Clustering algorithm =================" << std::endl;
   }
  
   do
   {
      continue_flag=this->assign_clusters();
      if(continue_flag)
      {
	 this->recalculate_centroids();	
	 if(debug)
	 {
	    std::cout << "\n-----------------------Iteration "<< iteration
		      <<"-----------------------\n\n\n";
	    this->print_cluster_info();
	    std::cout<< std::endl;
	    this->print_cluster_out();
	    std::cout << "\n---------------------------------------------------------\n\n";
	 }
	 ++iteration;
      }
   } while(continue_flag);

   if(debug)
   {
      std::cout << "=========================================================\n" << std::endl;
      std::cout << "\n=================== Final cluster info ==================\n" << std::endl;
      this->print_cluster_info(); std::cout<< std::endl; this->print_cluster_out();
      std::cout << "\n=========================================================\n" << std::endl;
      std::cout << "*********************************************************\n" << std::endl;
   }
   return this->cluster_out;
}

template <typename T, int N>
const std::vector<KM_Cluster<double,N>>& KM_Clusterer<T,N>::get_clusters()
{
   return this->cluster_vector;
}

template <typename T, int N>
void KM_Clusterer<T,N>::print_cluster_info()
{
   for(int i=0; i<(this->K); ++i)
   {
      std::cout <<"Centroid["<< i <<"]: "<< this->cluster_vector[i].centroid()
		<<"  Size="<< this->cluster_vector[i].size() << std::endl;
   }
}

//-------------------------------------------------------------------
//Private methods
//-------------------------------------------------------------------

template <typename T, int N>
void KM_Clusterer<T,N>::assign_random_clusters()
{
   int cluster_id;
   for(int i=0; i<this->cluster_out.size(); ++i)
   {
      cluster_id = rand()%(this->K);
      this->cluster_out[i] = cluster_id;
      this->cluster_vector[cluster_id].size() += 1;
   }
}

template <typename T, int N>
bool KM_Clusterer<T,N>::assign_clusters()
{    
   bool change_flag=false;
   int cluster_id=0;

   for(int i=0; i<(this->K); ++i) this->cluster_vector[i].size() = 0;

   int i=0;
   for(auto it=this->data_ptr->begin(); it!=this->data_ptr->end(); ++it)
   {
      cluster_id = nearest_centroid(*it);
      
      if(this->cluster_out[i] != cluster_id)
      {
	 change_flag=true;
	 this->cluster_out[i] = cluster_id;
      }
      this->cluster_vector[cluster_id].size()+=1;
      ++i;
   }
   return change_flag;
}

template <typename T, int N>
void KM_Clusterer<T,N>::call_from_thread(int tid) {
	std::cout << "Launched by thread " << tid << std::endl;
}


template <typename T, int N>
int KM_Clusterer<T,N>::nearest_centroid(const RowVector<T,N>& data_point)
{
   static const int num_threads = 3;
   int cluster_id=0;
   double curr_dist, curr_min_dist;
   std::thread t[num_threads];
   for(int i =0; i<num_threads; i++){
	t[i] = std::thread(call_from_thread, i);	
   }
   std::cout << "Launched from the main\n" << std::endl;
   for(int i = 0; i<num_threads; i++){
	t[i].join();
   }

   for(int i=0; i<K; ++i)
   {
      //For now centroid type is fixed to double => data_point casting
      curr_dist = distance::euclidean(RowVector<double,N>(data_point.template cast<double>()),
				      this->cluster_vector[i].centroid());   
      if(i==0) curr_min_dist = curr_dist;
      else if(curr_dist < curr_min_dist)
      {
	 curr_min_dist = curr_dist;
	 cluster_id = i;
      }
   }
   return cluster_id;
}

template <typename T, int N>
void KM_Clusterer<T,N>::recalculate_centroids()
{    
   for(int i=0; i<(this->K); ++i) this->cluster_vector[i].centroid().setZero();

   int centroid_id;
   int i=0;
   for(auto it=this->data_ptr->begin(); it!=this->data_ptr->end(); ++it)
   {
      centroid_id = this->cluster_out[i];
      this->cluster_vector[centroid_id].centroid() += (*it).template cast<double>(); 
      ++i;
   }
    
   for(int i=0; i<K; ++i)
   {
      this->cluster_vector[i].centroid() /= this->cluster_vector[i].size();
   }    
}

#endif
